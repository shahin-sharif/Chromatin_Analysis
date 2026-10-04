import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from . import VERSION
from .common import require,load,dump,sha,write_tsv
from .inputs import configuration,model_data,comparison_data,expression_data
from .regions import build
from .analysis import summarize,gene_sets,enrichment,genomic_context
from .report import report


def prepare(config_path):
    cfg=configuration(config_path)
    models={name:model_data(m) for name,m in cfg['models'].items()}
    regions,sizes=build(cfg)
    entities={r['entity_id']:r['gene_id'] for r in regions if not r['entity_id'].startswith('custom:')}
    comparisons={};expressions={};audits={}
    paths={str(Path(config_path).resolve()),cfg['gtf'],cfg['chrom_sizes']}
    for key in ('entity_ids',):
        if cfg.get(key):paths.add(cfg[key])
    if cfg['expression'].get('id_map'):paths.add(cfg['expression']['id_map'])
    for custom in cfg.get('custom_regions',[]):paths.add(custom['path'])
    for model in cfg['models'].values():
        for key in ('emissions','state_annotations','state_id_map'):
            if model.get(key):paths.add(model[key])
        paths.update(x['path'] for x in model.get('context_tables',[]))
    for c in cfg['comparisons']:
        comp=comparison_data(c,models)
        require(all(chrom in sizes and b<=sizes[chrom] for chrom,a,b,w,m in comp['tiles']), 'Shared segmentation territory exceeds/misses reference dictionary')
        comparisons[c['name']]=comp
        expressions[c['name']],audits[c['name']]=expression_data(c.get('expression'),cfg,entities)
        paths.update(comp['paths'].values())
        paths.update(str(p.resolve()) for p in Path(c['directory']).iterdir() if p.is_file() and (p.name=='RUNINFO.json' or p.suffix=='.tsv'))
        if c.get('expression'):paths.add(c['expression'])
        for selection in cfg['selections']:
            require(all(w in comp['wt_states'] and m in comp['mt_states'] for w,m in selection['pairs']),'Selected state absent from '+c['name'])
    sets=names=None
    if cfg.get('enrichment'):
        sets,names=gene_sets(cfg['enrichment']);paths.add(cfg['enrichment']['gene_sets'])
        require(set(entities.values()) & set().union(*sets.values()),'Gene-set IDs do not match any annotation gene IDs; supply matched memberships')
    if cfg['plots']:
        import matplotlib
    # Source code is part of the fingerprint so resume never mixes software versions.
    paths.update(str(p.resolve()) for p in Path(__file__).parent.glob('*.py'))
    paths.add(str((Path(__file__).parent.parent/'ChromHMMTools.py').resolve()))
    paths.add(str((Path(__file__).parent.parent/'ChromatinStateAnalysis.py').resolve()))
    inputs={p:sha(p) for p in sorted(paths)}
    digest=__import__('hashlib').sha256(json.dumps(dict(config=cfg,inputs=inputs),sort_keys=True).encode()).hexdigest()
    output=Path(cfg['outdir'])
    require(all(output!=Path(p) and output not in Path(p).parents for p in inputs),'Output contains/aliases an input')
    return cfg,models,regions,comparisons,expressions,audits,sets,names,inputs,digest


def stage(root,name,writer,resume=False):
    dest=root/name
    if dest.exists():
        require(resume,'Stage already exists: '+name)
        manifest=load(dest/'stage.json')
        require(manifest.get('complete') and all((dest/p).is_file() and sha(dest/p)==h for p,h in manifest['files'].items()),'Modified/incomplete stage '+name+'; preserve results and choose a new output directory')
        print('Reusing verified stage '+name,flush=True)
        return
    temp=Path(tempfile.mkdtemp(prefix='.'+name+'-',dir=root))
    try:
        writer(temp)
        files={str(p.relative_to(temp)):sha(p) for p in temp.rglob('*') if p.is_file()}
        dump(temp/'stage.json',dict(complete=True,files=files))
        temp.rename(dest)
    finally:
        if temp.exists():shutil.rmtree(temp)
    print('Completed stage '+name,flush=True)


def run(config_path,resume=False):
    prepared=prepare(config_path)
    cfg,models,regions,comparisons,expressions,audits,sets,names,inputs,digest=prepared
    root=Path(cfg['outdir'])
    if root.exists():
        require(resume and (root/'manifest.json').is_file(),'Output exists; use a new directory or --resume on a matching run')
        require(load(root/'manifest.json')['fingerprint']==digest,'Inputs/configuration/software changed; use a new output directory')
    else:
        root.mkdir(parents=True)
        import numpy
        versions=dict(python=sys.version,numpy=numpy.__version__)
        if cfg['plots']:
            import matplotlib
            versions['matplotlib']=matplotlib.__version__
        dump(root/'manifest.json',dict(tool='ChromatinStateAnalysis',version=VERSION,fingerprint=digest,
             inputs=inputs,configuration=cfg,software=versions,
             conventions='BED half-open; TSS boundary windows; shared coverage denominator; no differential chromatin significance testing'))
    # Prevent concurrent writers. A stale lock after an abrupt kill must be reviewed.
    lock=root/'.run.lock'
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        (root/'INCOMPLETE.txt').write_text('Requested stages have not all completed.\n')
        if (root/'COMPLETE.txt').exists():(root/'COMPLETE.txt').unlink()
        def region_writer(out):
            dump(out/'regions.json',regions)
            write_tsv(out/'region_qc.tsv',regions,[k for k in regions[0] if k!='blocks'])
            with (out/'regions.bed').open('w') as f:
                for r in regions:
                    for i,(a,b) in enumerate(r['blocks'],1):
                        f.write(f'{r["chrom"]}\t{a}\t{b}\t{r["region_id"]}.block{i}\t0\t{r["strand"]}\n')
        stage(root,'regions',region_writer,resume)
        stage(root,'models',lambda out:dump(out/'models.json',models),resume)
        for c in cfg['comparisons']:
            name=c['name'];comp=comparisons[name]
            def comparison_writer(out):
                metrics,composition,aggregates,states,selected=summarize(regions,comp,expressions[name],cfg,bool(c.get('expression')))
                context=genomic_context(comp['tiles'],cfg)
                from collections import Counter
                data=dict(expression_id_counts=dict(Counter(r['status'] for r in audits[name])),genomic_context=context,metrics=metrics,composition=composition,aggregates=aggregates,state_rows=states,selected=selected,
                          upstream=comp['info'],wt_states=comp['wt_states'],mt_states=comp['mt_states'],global_matrix=comp['matrix'])
                dump(out/'analysis.json',data)
                with (out/'all_overlap_tiles.bed').open('w') as f:
                    for chrom,a,b,w,m in comp['tiles']:f.write(f'{chrom}\t{a}\t{b}\t{w}\t{m}\n')
                with (out/'selected_transitions.bed').open('w') as f:
                    for sel in cfg['selections']:
                        pairs={tuple(p) for p in sel['pairs']}
                        for chrom,a,b,w,m in comp['tiles']:
                            if (w,m) in pairs:f.write(f'{chrom}\t{a}\t{b}\t{sel["name"]}:{w}->{m}\t0\t.\n')
                specifications=[('entity_metrics.tsv',metrics,list(metrics[0])),
                    ('entity_state_pairs.tsv',composition,['region_id','entity_id','window','wt_state','mt_state','bp','fraction','eligible']),
                    ('group_state_pairs.tsv',aggregates,['window','group_kind','group','wt_state','mt_state','bp','attributed_shared_bp','eligible_entities','bp_fraction','entity_mean_fraction','p_mt_given_wt']),
                    ('group_states.tsv',states,['window','group_kind','group','condition','state','bp','bp_fraction','entity_mean_fraction','eligible_entities']),
                    ('selected_entities.tsv',selected,list(metrics[0])+['selection','selected_bp','selected_fraction']),
                    ('expression_id_audit.tsv',audits[name],['expression_id','entity_id','status'])]
                for filename,rows,columns in specifications:write_tsv(out/filename,rows,columns)
                write_tsv(out/'genomic_context.tsv',context,['selection','category','bp','total_bp','fraction'])
                # Preserve global state coverage with explicitly different denominators.
                coverage=[]
                for condition,key in [('WT','wt_coverage'),('MT','mt_coverage')]:
                    total=sum(comp[key].values())
                    for state,bp in comp[key].items():coverage.append(dict(condition=condition,state=state,bp=bp,condition_covered_bp=total,fraction=bp/total))
                write_tsv(out/'global_state_coverage.tsv',coverage,['condition','state','bp','condition_covered_bp','fraction'])
            stage(root,'comparison_'+name,comparison_writer,resume)
            if cfg.get('enrichment'):
                def enrichment_writer(out):
                    data=load(root/('comparison_'+name)/'analysis.json')
                    results,backgrounds,universes=enrichment(data['metrics'],data['selected'],cfg,sets,names,bool(c.get('expression')))
                    dump(out/'enrichment.json',dict(results=results,backgrounds=backgrounds,resource=cfg['enrichment']))
                    write_tsv(out/'go_enrichment.tsv',results,['selection','window','term_id','term_name','overlap_genes','selected_genes','term_background','background_genes','fold_enrichment','pvalue','padj','genes'])
                    write_tsv(out/'background_audit.tsv',backgrounds,['selection','window','eligible_genes','annotated_background','selected_annotated_genes','unmapped_background'])
                    write_tsv(out/'gene_universes.tsv',universes,['selection','window','gene_id','annotated','selected'])
                stage(root,'enrichment_'+name,enrichment_writer,resume)
        def concordance_writer(out):
            all_rows={}
            for c in cfg['comparisons']:
                for r in load(root/('comparison_'+c['name'])/'analysis.json')['metrics']:
                    key=(r['entity_id'],r['window']);row=all_rows.setdefault(key,dict(entity_id=key[0],window=key[1]))
                    for field in ('changed_fraction','dominant_pair','expression_group','eligible'):
                        row[c['name']+':'+field]=r[field]
            rows=list(all_rows.values())
            for row in rows:
                groups=[row[c['name']+':expression_group'] for c in cfg['comparisons']]
                row['consistent_expression_direction']=len(groups)>1 and len(set(groups))==1 and groups[0] in ('up','down')
                comparable=len(cfg['comparisons'])>1 and len({c['reference_model'] for c in cfg['comparisons']})==1 and all(comparisons[c['name']]['info']['state_mode'] in ('shared','mapped') for c in cfg['comparisons'])
                pairs=[row[c['name']+':dominant_pair'] for c in cfg['comparisons']]
                supported=comparable and all(row[c['name']+':eligible'] for c in cfg['comparisons']) and all(p not in (None,'tie') for p in pairs)
                row['consistent_dominant_pair']=len(set(pairs))==1 if supported else None
            dump(out/'concordance.json',rows)
            write_tsv(out/'concordance.tsv',rows,list(rows[0]))
        stage(root,'concordance',concordance_writer,resume)
        stage(root,'report',lambda out:report(root,out,cfg),resume)
        (root/'COMPLETE.txt').write_text('All requested post-segmentation stages completed. Review QC and descriptive-analysis limits.\n')
        (root/'INCOMPLETE.txt').unlink()
        print('Completed. Report: '+str(root/'report'/'report.html'))
    finally:
        lock.unlink(missing_ok=True)


def main(argv=None):
    parser=argparse.ArgumentParser(description='Post-segmentation analysis of ChromHMM models and verified ChromHMMTools outputs')
    parser.add_argument('--version',action='version',version=VERSION)
    sub=parser.add_subparsers(dest='command',required=True)
    for cmd in ('check','run'):
        p=sub.add_parser(cmd);p.add_argument('--config',required=True)
        if cmd=='run':p.add_argument('--resume',action='store_true')
    p=sub.add_parser('report');p.add_argument('--results',required=True);p.add_argument('--outdir',required=True)
    args=parser.parse_args(argv)
    try:
        if args.command=='check':
            cfg,models,regions,comparisons,*_=prepare(args.config)
            print(f'Preflight passed: {len(models)} models, {len(comparisons)} comparisons, {len(regions)} windows. No outputs created.')
        elif args.command=='run':run(args.config,args.resume)
        else:
            root=Path(args.results).resolve();dest=Path(args.outdir).resolve()
            require(not dest.exists(),'Report destination exists; choose a new directory')
            manifest=load(root/'manifest.json')
            require(manifest.get('tool')=='ChromatinStateAnalysis','Not a pipeline result directory')
            # Verify every persisted numerical stage before rebuilding the report.
            for folder in root.iterdir():
                if folder.is_dir() and not folder.name.startswith('.') and folder.name!='report':
                    marker=load(folder/'stage.json')
                    require(marker.get('complete') and all((folder/p).is_file() and sha(folder/p)==h for p,h in marker['files'].items()),'Changed/incomplete stage: '+folder.name)
            dest.parent.mkdir(parents=True,exist_ok=True)
            temp=Path(tempfile.mkdtemp(prefix='.report-',dir=dest.parent))
            try:
                report(root,temp,manifest['configuration']);temp.rename(dest)
            finally:
                if temp.exists():shutil.rmtree(temp)
            print('Report rebuilt: '+str(dest/'report.html'))
        return 0
    except (OSError,ValueError,KeyError,TypeError,ImportError,OverflowError,RuntimeError) as error:
        print('ERROR: '+str(error),file=sys.stderr);return 2
