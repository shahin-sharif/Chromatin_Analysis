"""Coverage-aware entity summaries and local, reproducible enrichment."""
import math
from collections import Counter,defaultdict
from .common import require, read_tsv
from .regions import TileIndex


def summarize(regions,comparison,expression,config,has_expression):
    index=TileIndex(comparison['tiles']);metrics=[];composition=[];aggregates=[];selected=[]
    counts_by_region={};groups=defaultdict(list)
    for region in regions:
        counts=Counter()
        for a,b,w,m in index.intersect(region):counts[(w,m)]+=b-a
        covered=sum(counts.values());effective=region['effective_bp']
        fraction=covered/effective if effective else None
        eligible=covered>0 and fraction>=config['min_coverage']
        ex=expression.get(region.get('expression_entity_id',region['entity_id']),dict(log2FoldChange=None,padj=None,abundance=None,
            expression_group='unmatched' if has_expression else 'not_supplied',abundance_group='abundance_missing'))
        # Proximity expression joins require an explicit per-custom-region option.
        changed=None if comparison['info']['state_mode']=='distinct' else sum(n for (w,m),n in counts.items() if w!=m)
        row={k:region[k] for k in ('region_id','entity_id','gene_id','window','association','effective_bp','requested_bp','clipped_bp')}
        row.update(ex,expression_entity_id=region.get('expression_entity_id',region['entity_id']),shared_bp=covered,coverage_fraction=fraction,eligible=eligible,
                   changed_bp=changed,changed_fraction=changed/covered if changed is not None and covered else None,
                   status='ok' if eligible else 'no_shared_coverage' if not covered else 'insufficient_coverage')
        top=max(counts.values(),default=0)
        winners=sorted(f'{w}->{m}' for (w,m),n in counts.items() if n==top)
        row['dominant_pair']=winners[0] if len(winners)==1 else 'tie' if winners else None
        metrics.append(row);counts_by_region[region['region_id']]=counts
        for (w,m),n in sorted(counts.items()):
            composition.append(dict(region_id=region['region_id'],entity_id=region['entity_id'],window=region['window'],
                                    wt_state=w,mt_state=m,bp=n,fraction=n/covered,eligible=eligible))
        if eligible:
            for kind,label in [('all','all'),('expression',ex['expression_group']),('abundance',ex['abundance_group'])]:
                groups[(region['window'],kind,label)].append((row,counts))
    for (window,kind,label),members in sorted(groups.items()):
        total=sum(r['shared_bp'] for r,c in members)
        for w in comparison['wt_states']:
            for m in comparison['mt_states']:
                bp=sum(c[(w,m)] for r,c in members)
                equal=sum(c[(w,m)]/r['shared_bp'] for r,c in members)/len(members)
                wt_support=sum(n for r,c in members for (a,b),n in c.items() if a==w)
                aggregates.append(dict(window=window,group_kind=kind,group=label,wt_state=w,mt_state=m,
                    bp=bp,attributed_shared_bp=total,eligible_entities=len(members),bp_fraction=bp/total,
                    entity_mean_fraction=equal,p_mt_given_wt=bp/wt_support if wt_support else None))
    for selection in config['selections']:
        pairs={tuple(p) for p in selection['pairs']}
        require(all(w in comparison['wt_states'] and m in comparison['mt_states'] for w,m in pairs),
                'Selection state absent from comparison axes: '+selection['name'])
        for row in metrics:
            if not row['eligible'] or (selection.get('windows') and row['window'] not in selection['windows']):continue
            if selection.get('expression_groups') and row['expression_group'] not in selection['expression_groups']:continue
            bp=sum(counts_by_region[row['region_id']][p] for p in pairs)
            if bp>=selection['min_bp'] and bp/row['shared_bp']>=selection['min_fraction']:
                selected.append(dict(**row,selection=selection['name'],selected_bp=bp,selected_fraction=bp/row['shared_bp']))
    state_rows=[]
    for condition,field,states in [('WT','wt_state',comparison['wt_states']),('MT','mt_state',comparison['mt_states'])]:
        for key,members in sorted(groups.items()):
            for state in states:
                values=[a for a in aggregates if (a['window'],a['group_kind'],a['group'])==key and a[field]==state]
                state_rows.append(dict(window=key[0],group_kind=key[1],group=key[2],condition=condition,state=state,
                    bp=sum(v['bp'] for v in values),bp_fraction=sum(v['bp_fraction'] for v in values),
                    entity_mean_fraction=sum(v['entity_mean_fraction'] for v in values),eligible_entities=len(members)))
    return metrics,composition,aggregates,state_rows,selected


def hypergeom_tail(k,N,K,n):
    """P[X >= k] for sampling n genes without replacement; stable log arithmetic."""
    if k<=max(0,n-(N-K)):return 1.0
    if k>min(K,n):return 0.0
    def logchoose(a,b):
        return math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
    logs=[logchoose(K,i)+logchoose(N-K,n-i)-logchoose(N,n) for i in range(k,min(K,n)+1) if 0<=n-i<=N-K]
    peak=max(logs)
    return min(1.,math.exp(peak)*sum(math.exp(x-peak) for x in logs))


def bh(values):
    result=[0.]*len(values);previous=1.
    for rank,i in reversed(list(enumerate(sorted(range(len(values)),key=values.__getitem__),1))):
        previous=min(previous,values[i]*len(values)/rank);result[i]=previous
    return result


def gene_sets(config):
    result={};descriptions={}
    for row in read_tsv(config['gene_sets'],['term_id','term_name','gene_id']):
        require(all(row[k] for k in ('term_id','term_name','gene_id')),'Empty gene-set field')
        term=row['term_id']
        require(term not in descriptions or descriptions[term]==row['term_name'],'Conflicting term names')
        descriptions[term]=row['term_name'];result.setdefault(term,set()).add(row['gene_id'])
    require(result,'Empty gene-set membership file')
    return result,descriptions


def enrichment(metrics,selected,config,sets,names,has_expression):
    results=[];backgrounds=[];universes=[]
    annotated=set().union(*sets.values())
    for selection in config['selections']:
        windows=selection.get('windows') or sorted({r['window'] for r in metrics})
        for window in windows:
            candidates=[r for r in metrics if r['window']==window and r['eligible'] and r['gene_id']
                        and (not selection.get('expression_groups') or r['expression_group'] in selection['expression_groups'])
                        and (not has_expression or r['expression_group'] in ('up','down','not_significant'))]
            eligible={r['gene_id'] for r in candidates}
            universe=eligible & annotated
            chosen={r['gene_id'] for r in selected if r['selection']==selection['name'] and r['window']==window}&universe
            audit=dict(selection=selection['name'],window=window,eligible_genes=len(eligible),annotated_background=len(universe),
                       selected_annotated_genes=len(chosen),unmapped_background=len(eligible-universe))
            backgrounds.append(audit)
            universes.extend(dict(selection=selection['name'],window=window,gene_id=g,
                                  annotated=g in annotated,selected=g in chosen) for g in sorted(eligible))
            if not universe or not chosen:continue
            family=[]
            for term,members in sorted(sets.items()):
                K=len(members&universe)
                if not config['enrichment']['min_size']<=K<=config['enrichment']['max_size']:continue
                overlap=members&chosen;k=len(overlap);N=len(universe);n=len(chosen)
                # Include zero-hit eligible terms in this BH family.
                family.append(dict(selection=selection['name'],window=window,term_id=term,term_name=names[term],
                    overlap_genes=k,selected_genes=n,term_background=K,background_genes=N,
                    fold_enrichment=(k/n)/(K/N),pvalue=hypergeom_tail(k,N,K,n),genes=';'.join(sorted(overlap))))
            for row,padj in zip(family,bh([r['pvalue'] for r in family])):
                row['padj']=padj;results.append(row)
    return results,backgrounds,universes


def genomic_context(tiles,config):
    """Disjoint base-pair feature annotation, separate from entity attribution."""
    from ChromHMMTools import annotation_partition,split_annotation
    from .regions import models,read_sizes
    genes,_=models(config['gtf'],'gene')
    sizes=read_sizes(config['chrom_sizes']);features=defaultdict(list)
    left,right=config['annotation_promoter']
    for g in genes.values():
        chrom=g['chrom'];limit=sizes.get(chrom,0)
        if not limit:continue
        a,b=g['start'],min(g['end'],limit)
        if a<b:features[chrom].append((a,b,'intron' if g['exons'] else 'genic_unresolved'))
        for a,b in g['exons']:
            if a<min(b,limit):features[chrom].append((a,min(b,limit),'exon'))
        tss=g['start'] if g['strand']=='+' else g['end']
        a,b=(tss+left,tss+right) if g['strand']=='+' else (tss-right,tss-left)
        if max(0,a)<min(limit,b):features[chrom].append((max(0,a),min(limit,b),'promoter'))
    part=annotation_partition(features,['promoter','exon','intron','genic_unresolved','intergenic'])
    partition={c:(r,[x[1] for x in r]) for c,r in part.items()}
    counters={'all_shared':Counter()}
    selections={s['name']:{tuple(p) for p in s['pairs']} for s in config['selections']}
    counters.update({name:Counter() for name in selections})
    for chrom,a,b,w,m in tiles:
        for l,r,category in split_annotation(chrom,a,b,partition):
            counters['all_shared'][category]+=r-l
            for name,pairs in selections.items():
                if (w,m) in pairs:counters[name][category]+=r-l
    return [dict(selection=name,category=category,bp=counter[category],total_bp=sum(counter.values()),
                 fraction=counter[category]/sum(counter.values()) if sum(counter.values()) else None)
            for name,counter in counters.items() for category in ('promoter','exon','intron','genic_unresolved','intergenic')]
