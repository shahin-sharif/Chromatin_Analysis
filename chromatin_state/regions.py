"""Explicit half-open windows and auditable feature identity."""
import bisect
import re
from collections import defaultdict
from .common import require, text
from ChromHMMTools import merge


def read_sizes(path):
    sizes = {}
    with text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            chrom, size = line.split()
            require(chrom not in sizes and int(size)>0, 'Invalid/duplicate chromosome size')
            sizes[chrom] = int(size)
    require(sizes, 'Empty chromosome dictionary')
    return sizes


def models(path, level, selected=None):
    genes, transcripts = {}, {}
    with text(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            require(len(fields)==9, f'GTF row {n}: expected nine columns')
            chrom, _, kind, start, end, _, strand, _, attributes = fields
            if kind not in ('gene', 'transcript', 'exon'):
                continue
            start, end = int(start)-1, int(end)
            require(0<=start<end and strand in ('+','-'), f'GTF row {n}: invalid coordinates/strand')
            attr = dict(re.findall(r'(\S+)\s+"([^"]*)"\s*;?', attributes))
            gene, tx = attr.get('gene_id'), attr.get('transcript_id')
            require(gene, f'GTF row {n}: missing gene_id')
            for table, identity, declaration in ((genes,gene,'gene'), (transcripts,tx,'transcript')):
                if table is transcripts and kind == 'gene':
                    continue
                require(identity, f'GTF row {n}: missing transcript_id')
                entry = table.setdefault(identity, dict(entity_id=identity, gene_id=gene,
                    gene_name=attr.get('gene_name',''), chrom=chrom, strand=strand,
                    start=start, end=end, declaration=None, exons=[]))
                require((entry['chrom'],entry['strand'],entry['gene_id'])==(chrom,strand,gene), 'Conflicting GTF identity: '+identity)
                entry['start'] = min(entry['start'],start)
                entry['end'] = max(entry['end'],end)
                if kind == declaration:
                    require(entry['declaration'] is None, 'Duplicate declaration: '+identity)
                    entry['declaration'] = (start,end)
                if kind == 'exon':
                    entry['exons'].append((start,end))
    chosen = genes if level=='gene' else transcripts
    if selected:
        with text(selected) as handle:
            ids = [line.strip() for line in handle if line.strip() and not line.startswith('#')]
        require(ids and len(ids)==len(set(ids)) and set(ids)<=chosen.keys(), 'Selection has duplicate/unknown/empty entity IDs')
        chosen = {i:chosen[i] for i in ids}
    for entry in chosen.values():
        declared = entry.pop('declaration')
        require(declared is None or declared==(entry['start'],entry['end']), 'Features outside declared span: '+entry['entity_id'])
        entry['boundary_source'] = 'declared_'+level if declared else 'feature_span'
        entry['exons'] = merge(entry['exons'])
    require(chosen, 'No annotation entities')
    return chosen, genes


def build(config):
    sizes = read_sizes(config['chrom_sizes'])
    entities, all_genes = models(config['gtf'],config['entity_level'],config.get('entity_ids'))
    require(any(e['chrom'] in sizes for e in entities.values()), 'No GTF/reference chromosome names match')
    regions = []
    def add(entry, name, blocks, association):
        requested = sum(b-a for a,b in blocks)
        limit = sizes.get(entry['chrom'],0)
        blocks = merge([(max(0,a),min(limit,b)) for a,b in blocks if min(limit,b)>max(0,a)])
        effective = sum(b-a for a,b in blocks)
        row = {k:entry[k] for k in ('entity_id','gene_id','gene_name','chrom','strand','boundary_source')}
        row.update(region_id=f'r{len(regions)+1:08d}', window=name, blocks=blocks,
                   expression_entity_id=entry.get('expression_entity_id',entry['entity_id']),
                   requested_bp=requested,effective_bp=effective,clipped_bp=requested-effective,
                   association=association,status='ok' if effective else 'empty_or_missing_chromosome')
        regions.append(row)
    for entry in entities.values():
        plus = entry['strand']=='+'
        tss, tes = (entry['start'],entry['end']) if plus else (entry['end'],entry['start'])
        for window in config['windows']:
            if window['kind']=='body':
                blocks = entry['exons'] if window.get('mode','span')=='exons' else [(entry['start'],entry['end'])]
            else:
                anchor = tss if window['anchor']=='tss' else tes
                a,b = window['start'],window['end']
                blocks = [(anchor+a,anchor+b)] if plus else [(anchor-b,anchor-a)]
            add(entry,window['name'],blocks,'annotated_'+config['entity_level'])
    # Nearest-TSS assignment is optional for supplied distal BED regions. All ties
    # are reported as ambiguous rather than choosing an arbitrary target gene.
    anchors = defaultdict(list)
    for gene, e in all_genes.items():
        anchors[e['chrom']].append((e['start'] if e['strand']=='+' else e['end'],gene))
    for chrom in anchors:
        anchors[chrom].sort()
    anchor_positions={chrom:[p for p,g in items] for chrom,items in anchors.items()}
    for custom in config.get('custom_regions',[]):
        seen = set()
        with text(custom['path']) as handle:
            for n,line in enumerate(handle,1):
                if not line.strip() or line.startswith(('#','track ','browser ')): continue
                x = line.rstrip('\n').split('\t')
                require(len(x) in (4,5,6), 'Custom BED requires 4–6 columns (BED12 is not supported here)')
                a,b = int(x[1]),int(x[2]); name=x[3]
                require(0<=a<b and name and name not in seen,'Invalid/duplicate custom BED row')
                seen.add(name)
                gene = ''; association='unassigned_custom'
                if custom.get('nearest_tss',False):
                    items = anchors[x[0]]; positions=anchor_positions.get(x[0],[])
                    # Distance from interval to TSS boundary; overlapping boundary has distance zero.
                    lo,hi=bisect.bisect_left(positions,a),bisect.bisect_left(positions,b)
                    candidates=items[lo:hi]
                    if not candidates:
                        candidates=items[max(0,lo-1):min(len(items),lo+1)]
                        if candidates:
                            distance=min(a-p if p<a else p-b for p,g in candidates)
                            positions_to_keep={p for p,g in candidates if (a-p if p<a else p-b)==distance}
                            candidates=[item for pos in positions_to_keep for item in items[bisect.bisect_left(positions,pos):bisect.bisect_right(positions,pos)]]
                    ids=sorted({g for p,g in candidates})
                    if len(ids)==1:
                        gene=ids[0]; association='nearest_tss_proximity'
                    elif ids:
                        association='ambiguous_nearest_tss'
                entry=dict(entity_id='custom:'+custom['name']+':'+name,gene_id=gene,gene_name=all_genes[gene]['gene_name'] if gene else '',
                           chrom=x[0],strand=x[5] if len(x)==6 else '.',boundary_source='custom_BED',
                           expression_entity_id=gene if custom.get('join_nearest_expression') and gene else '')
                require(entry['strand'] in ('+','-','.'),'Invalid custom BED strand')
                add(entry,custom['name'],[(a,b)],association)
    return regions,sizes


class TileIndex:
    def __init__(self,tiles):
        self.rows=defaultdict(list)
        for chrom,a,b,w,m in tiles:
            self.rows[chrom].append((a,b,w,m))
        self.ends={c:[r[1] for r in rows] for c,rows in self.rows.items()}

    def intersect(self,region):
        rows=self.rows[region['chrom']]; ends=self.ends.get(region['chrom'],[])
        for a,b in region['blocks']:
            i=bisect.bisect_right(ends,a)
            while i<len(rows) and rows[i][0]<b:
                c,d,w,m=rows[i]
                yield max(a,c),min(b,d),w,m
                i+=1
