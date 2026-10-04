"""Configuration validation and verified reuse of ChromHMMTools comparisons."""
from pathlib import Path
from collections import Counter
import math
import numpy as np
import ChromHMMTools as core
from .common import require, keys, safe_name, load, sha, read_tsv, number

DEFAULT_WINDOWS = [
    dict(name='tss_core',kind='window',anchor='tss',start=-25,end=25),
    dict(name='upstream',kind='window',anchor='tss',start=-1000,end=-25),
    dict(name='downstream',kind='window',anchor='tss',start=25,end=1000),
    dict(name='broad_promoter',kind='window',anchor='tss',start=-3000,end=3000),
    dict(name='body',kind='body',mode='span'),
    dict(name='tes',kind='window',anchor='tes',start=-500,end=500)]


def configuration(path):
    path=Path(path).resolve(); cfg=load(path)
    keys(cfg,['schema_version','outdir','gtf','chrom_sizes','entity_level','entity_ids','windows',
              'custom_regions','comparisons','models','expression','min_coverage','selections',
              'enrichment','plots','abundance_breaks','annotation_promoter'],['schema_version','outdir','gtf','chrom_sizes','comparisons','models'])
    require(cfg['schema_version']==1,'Unsupported configuration version')
    base=path.parent
    def resolve(value):
        return str((base/Path(value).expanduser()).resolve())
    cfg['outdir']=resolve(cfg['outdir'])
    for key in ('gtf','chrom_sizes','entity_ids'):
        if cfg.get(key): cfg[key]=resolve(cfg[key])
    cfg.setdefault('entity_level','gene'); require(cfg['entity_level'] in ('gene','transcript'),'Invalid entity_level')
    cfg.setdefault('windows',DEFAULT_WINDOWS)
    require(cfg['windows'],'At least one window required')
    names=[]
    for w in cfg['windows']:
        keys(w,['name','kind','anchor','start','end','mode'],['name','kind'])
        names.append(safe_name(w['name']))
        if w['kind']=='window':
            require(w.get('anchor') in ('tss','tes') and all(type(w.get(k)) is int for k in ('start','end'))
                    and w['start']<w['end'] and 'mode' not in w,'Invalid anchor window')
        else:
            require(w['kind']=='body' and w.get('mode','span') in ('span','exons') and not ({'anchor','start','end'}&w.keys()),'Invalid body window')
    for custom in cfg.get('custom_regions',[]):
        keys(custom,['name','path','nearest_tss','join_nearest_expression'],['name','path'])
        names.append(safe_name(custom['name']));custom['path']=resolve(custom['path'])
        require(type(custom.get('nearest_tss',False)) is bool and type(custom.get('join_nearest_expression',False)) is bool,'Custom association flags must be boolean')
        require(not custom.get('join_nearest_expression') or (custom.get('nearest_tss') and cfg['entity_level']=='gene'),'Nearest expression join requires nearest_tss and gene-level expression')
    require(len(names)==len(set(names)),'Window names must be unique')
    cfg.setdefault('annotation_promoter',[-1000,200])
    require(isinstance(cfg['annotation_promoter'],list) and len(cfg['annotation_promoter'])==2 and all(type(x) is int for x in cfg['annotation_promoter']) and cfg['annotation_promoter'][0]<cfg['annotation_promoter'][1], 'Invalid annotation_promoter offsets')
    cfg.setdefault('min_coverage',0.8)
    require(isinstance(cfg['min_coverage'],(float,int)) and math.isfinite(cfg['min_coverage']) and 0<=cfg['min_coverage']<=1,'Invalid min_coverage')
    cfg.setdefault('plots',True);require(type(cfg['plots']) is bool,'plots must be boolean')
    cfg.setdefault('abundance_breaks',[100,500])
    breaks=cfg['abundance_breaks']
    require(isinstance(breaks,list) and all(type(v) in (float,int) and math.isfinite(v) and v>=0 for v in breaks)
            and sorted(set(breaks))==breaks,'Abundance breaks must be finite, nonnegative and increasing')
    cfg.setdefault('expression',{})
    ex=cfg['expression']
    keys(ex,['id_column','gene_column','log2fc_column','padj_column','abundance_column','fdr','min_abs_log2fc','id_map'])
    defaults=dict(id_column='entity_id',log2fc_column='log2FoldChange',padj_column='padj',abundance_column='baseMean',fdr=.05,min_abs_log2fc=0.)
    for k,v in defaults.items(): ex.setdefault(k,v)
    require(math.isfinite(ex['fdr']) and 0<ex['fdr']<=1 and math.isfinite(ex['min_abs_log2fc']) and ex['min_abs_log2fc']>=0,'Invalid expression thresholds')
    if ex.get('id_map'): ex['id_map']=resolve(ex['id_map'])
    require(cfg['models'] and isinstance(cfg['models'],dict),'Models must be a nonempty object')
    for name,model in cfg['models'].items():
        safe_name(name)
        keys(model,['emissions','state_annotations','state_id_map','context_tables'],['emissions'])
        model['emissions']=resolve(model['emissions'])
        for key in ('state_annotations','state_id_map'):
            if model.get(key): model[key]=resolve(model[key])
        for context in model.get('context_tables',[]):
            keys(context,['name','path','kind','value_unit','columns'],['name','path','kind','value_unit'])
            safe_name(context['name']);context['path']=resolve(context['path'])
            require(context['kind'] in ('overlap','neighborhood') and context['value_unit'] in ('fold_enrichment','fraction','percent','count'),'Invalid context table kind/units')
    names=[]
    require(isinstance(cfg['comparisons'],list) and cfg['comparisons'],'At least one comparison required')
    for c in cfg['comparisons']:
        keys(c,['name','directory','reference_model','target_model','expression','input_overrides'],['name','directory','reference_model','target_model'])
        names.append(safe_name(c['name']));c['directory']=resolve(c['directory'])
        require(c['reference_model'] in cfg['models'] and c['target_model'] in cfg['models'],'Unknown model')
        if c.get('expression'):c['expression']=resolve(c['expression'])
        for key,value in c.get('input_overrides',{}).items():c['input_overrides'][key]=resolve(value)
    require(len(names)==len(set(names)),'Duplicate comparison names')
    names=[]
    cfg.setdefault('selections',[])
    for sel in cfg['selections']:
        keys(sel,['name','pairs','min_fraction','min_bp','windows','expression_groups'],['name','pairs'])
        names.append(safe_name(sel['name']))
        require(sel['name']!='all_shared','Selection name all_shared is reserved')
        require(sel['pairs'] and all(isinstance(p,list) and len(p)==2 and all(isinstance(x,str) for x in p) for p in sel['pairs']),'Selections require explicit state pairs')
        require(len({tuple(p) for p in sel['pairs']})==len(sel['pairs']),'Duplicate selection pair')
        sel.setdefault('min_fraction',0.);sel.setdefault('min_bp',1)
        require(math.isfinite(sel['min_fraction']) and 0<=sel['min_fraction']<=1 and type(sel['min_bp']) is int and sel['min_bp']>0,'Invalid selection thresholds')
        require(set(sel.get('windows',[]))<=set(w['name'] for w in cfg['windows'])|set(w['name'] for w in cfg.get('custom_regions',[])),'Unknown selection window')
        require(set(sel.get('expression_groups',[]))<= {'up','down','not_significant','not_tested','unmatched','not_supplied'},'Invalid expression group')
    require(len(names)==len(set(names)),'Duplicate selections')
    if cfg.get('enrichment'):
        e=cfg['enrichment'];keys(e,['gene_sets','source','release','membership','min_size','max_size'],['gene_sets','source','release','membership'])
        e['gene_sets']=resolve(e['gene_sets'])
        require(e['source'] and e['release'] and e['membership']=='explicit_expanded','Supply source/release and membership=explicit_expanded (GO ancestors already expanded)')
        e.setdefault('min_size',5);e.setdefault('max_size',500)
        require(type(e['min_size']) is int and type(e['max_size']) is int and 0<e['min_size']<=e['max_size'],'Invalid term sizes')
        require(cfg['selections'],'Enrichment requires selections')
    return cfg


def model_data(model):
    matrix,states,marks=core.read_matrix(model['emissions'])
    require((matrix<=1).all(),'Emission probabilities must be in [0,1]')
    mapping={}
    if model.get('state_id_map'):
        for r in read_tsv(model['state_id_map'],['emission_state','segment_state']):
            require(r['emission_state'] not in mapping,'Duplicate emission-state mapping')
            mapping[r['emission_state']]=r['segment_state']
        require(set(mapping)==set(states) and len(set(mapping.values()))==len(states),'Emission state mapping must be complete and one-to-one')
        states=[mapping[s] for s in states]
    require(all(states),'Empty emission state label')
    labels={s:dict(state=s,label=s,color='#808080',order=i) for i,s in enumerate(states)}
    if model.get('state_annotations'):
        rows=read_tsv(model['state_annotations'],['state','label','color','order'])
        require(len(rows)==len(states) and {r['state'] for r in rows}==set(states),'State annotations must cover model exactly once')
        import re
        for r in rows:
            require(re.fullmatch(r'#[0-9a-fA-F]{6}',r['color']) and r['label'],'Invalid color/label')
            r['order']=int(r['order']);labels[r['state']]=r
        require(len({r['order'] for r in rows})==len(rows),'State order must be unique')
    contexts=[]
    names=set()
    for ctx in model.get('context_tables',[]):
        require(ctx['name'] not in names,'Duplicate model context name');names.add(ctx['name'])
        values,rs,cs=core.read_matrix(ctx['path'])
        rs=[mapping.get(s,s) for s in rs]
        if ctx.get('columns'):
            require(len(set(ctx['columns']))==len(ctx['columns']) and set(ctx['columns'])<=set(cs),'Invalid context column selection')
            values=values[:,[cs.index(s) for s in ctx['columns']]];cs=ctx['columns']
        require(ctx['value_unit']!='fold_enrichment' or not any('genome' in s.lower() and '%' in s for s in cs), 'Mixed context units: use columns to exclude Genome % from fold-enrichment panels')
        if ctx['value_unit']=='fraction': require((values<=1).all(),'Context fraction exceeds 1')
        if ctx['value_unit']=='percent': require((values<=100).all(),'Context percent exceeds 100')
        # ChromHMM overlap tables may contain a genome/background summary row.
        background=[s for s in rs if s not in states]
        require(set(states)<=set(rs) and all(s.lower() in ('base','genome','genome %','base (genome %)') for s in background),'Context table has missing/unknown state rows')
        contexts.append(dict(**ctx,states=rs,columns=cs,values=values.tolist(),background_rows=background))
    return dict(states=states,marks=marks,values=matrix.tolist(),annotations=labels,contexts=contexts)


def comparison_data(c,models):
    folder=Path(c['directory']);info=load(folder/'RUNINFO.json')
    require(info.get('version')=='1.0.0','Unsupported ChromHMMTools output version')
    require(info.get('state_mode') in ('shared','mapped','distinct'),'Invalid upstream state mode')
    paths={}
    overrides=c.get('input_overrides',{})
    require(set(overrides)<=info['inputs'].keys(),'Unknown input override')
    for name,entry in info['inputs'].items():
        path=Path(overrides.get(name,entry['path'])).expanduser()
        # Upstream relative paths depend on its original working directory; never guess.
        require(path.is_absolute(),f'Upstream {name} path is relative; supply input_overrides.{name}')
        require(path.is_file() and sha(path)==entry['sha256'],f'Upstream input missing/changed: {name}; restore exact input or rerun ChromHMMTools')
        paths[name]=str(path)
    wt=core.segmentation(paths['wt']); mt=core.segmentation(paths['mt'])
    require({r[2] for v in wt.values() for r in v}<=set(models[c['reference_model']]['states']),'WT states missing from emissions; supply explicit state_id_map')
    require({r[2] for v in mt.values() for r in v}<=set(models[c['target_model']]['states']),'KO states missing from emissions; supply explicit state_id_map')
    if info['state_mode']=='shared':
        require(c['reference_model']==c['target_model'],'Shared state mode requires one declared model; otherwise use mapped/distinct upstream')
    mapping={}
    if info['state_mode']=='mapped':
        with open(paths['mapping']) as handle:
            for line in handle:
                if line.strip() and not line.startswith('#'):
                    a,b=line.split(); require(a not in mapping,'Duplicate state mapping');mapping[a]=b
        require(set(mapping)=={r[2] for v in mt.values() for r in v} and len(set(mapping.values()))==len(mapping)
                and set(mapping.values())<={r[2] for v in wt.values() for r in v},'Invalid upstream state correspondence')
        mt={c:[(a,b,mapping[s]) for a,b,s in rows] for c,rows in mt.items()}
    settings=info['settings']
    if settings.get('autosomes_only'):
        keep={str(i) for i in range(1,23)}|{f'chr{i}' for i in range(1,23)}
        wt={c:r for c,r in wt.items() if c in keep};mt={c:r for c,r in mt.items() if c in keep}
    inc=core.read_bed(paths['include'],False) if 'include' in paths else None
    exc=core.read_bed(paths['exclude'],False) if 'exclude' in paths else None
    wt=core.clip_regions(wt,inc,exc);mt=core.clip_regions(mt,inc,exc)
    tiles=list(core.overlap_tiles(wt,mt))
    matrix,ws,ms=core.read_matrix(folder/'overlap.bp.tsv')
    counts=Counter()
    for chrom,a,b,w,m in tiles:counts[(w,m)]+=b-a
    require(set(w for w,m in counts)<=set(ws) and set(m for w,m in counts)<=set(ms),'Upstream matrix state axes differ from inputs')
    actual=np.array([[counts[(w,m)] for m in ms] for w in ws],dtype=np.int64)
    require(np.array_equal(actual,matrix) and actual.sum()>0 and int(actual.sum())==info['shared_bp'],
            'Reconstructed full overlaps disagree with saved comparison counts')
    # Retain every model state, including unobserved states, without positional padding.
    if info['state_mode'] in ('shared','mapped'):
        ws=ms=list(models[c['reference_model']]['states'])
    else:
        ws=list(models[c['reference_model']]['states']);ms=list(models[c['target_model']]['states'])
    actual=np.array([[counts[(w,m)] for m in ms] for w in ws],dtype=np.int64)
    return dict(info=info,paths=paths,tiles=tiles,wt_states=ws,mt_states=ms,matrix=actual.tolist(),
                mapping=mapping,wt_coverage={s:sum(b-a for rows in wt.values() for a,b,t in rows if s==t) for s in ws},
                mt_coverage={s:sum(b-a for rows in mt.values() for a,b,t in rows if s==t) for s in ms})


def expression_data(path,config,entities):
    if not path:return {},[]
    ex=config['expression'];columns=[ex['id_column'],ex['log2fc_column'],ex['padj_column']]
    if ex.get('abundance_column'):columns.append(ex['abundance_column'])
    if ex.get('gene_column'):columns.append(ex['gene_column'])
    mapping={}
    if ex.get('id_map'):
        for r in read_tsv(ex['id_map'],['expression_id','entity_id']):
            require(r['expression_id'] not in mapping,'Duplicate expression ID mapping')
            mapping[r['expression_id']]=r['entity_id']
        require(len(set(mapping.values()))==len(mapping),'Ambiguous many-to-one expression mapping; supply an explicit gene-level table')
    result={};audit=[];seen=set()
    for r in read_tsv(path,columns):
        original=r[ex['id_column']]
        require(original and original not in seen,'Empty/duplicate expression ID: '+original);seen.add(original)
        identity=mapping.get(original) if mapping else original
        matched=identity in entities
        audit.append(dict(expression_id=original,entity_id=identity,status='matched' if matched else 'unmatched'))
        lfc=number(r[ex['log2fc_column']]);padj=number(r[ex['padj_column']])
        abundance=number(r[ex['abundance_column']]) if ex.get('abundance_column') else None
        require(padj is None or 0<=padj<=1,'Adjusted p-value outside [0,1]')
        require(abundance is None or abundance>=0,'Negative expression abundance')
        if not matched:continue
        require(identity not in result,'Multiple expression rows map to one entity')
        if ex.get('gene_column'):require(r[ex['gene_column']]==entities[identity],'Expression gene ID disagrees with GTF')
        group='not_tested' if lfc is None or padj is None else 'not_significant'
        if group!='not_tested' and padj<=ex['fdr'] and abs(lfc)>=ex['min_abs_log2fc'] and lfc!=0:
            group='up' if lfc>0 else 'down'
        breaks=config['abundance_breaks']
        if abundance is None:abundance_group='abundance_missing'
        else:
            i=sum(abundance>=v for v in breaks)
            left=breaks[i-1] if i else 0
            right=breaks[i] if i<len(breaks) else 'inf'
            abundance_group=f'abundance_[{left},{right})'
        result[identity]=dict(log2FoldChange=lfc,padj=padj,abundance=abundance,
                              expression_group=group,abundance_group=abundance_group)
    require(result,'Expression table has no matching annotation IDs')
    for entity in entities.keys()-result.keys():
        audit.append(dict(expression_id=None,entity_id=entity,status='annotation_without_expression'))
    return result,audit
