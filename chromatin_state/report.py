"""Portable report generation from saved numerical results, without refitting."""
import base64
import html
import io
from pathlib import Path
from .common import load


def table(rows, columns, limit=25):
    esc=lambda x:html.escape('NA' if x is None else str(x))
    return '<table><tr>'+''.join('<th>'+esc(c)+'</th>' for c in columns)+'</tr>'+''.join(
        '<tr>'+''.join('<td>'+esc(r.get(c))+'</td>' for c in columns)+'</tr>' for r in rows[:limit])+'</table>'


def report(root,out,config):
    models=load(root/'models'/'models.json')
    pieces=['<!doctype html><html lang="en"><meta charset="utf-8"><title>Chromatin state analysis</title>',
            '<style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 20px;color:#182838}table{border-collapse:collapse;font-size:13px;display:block;overflow:auto}td,th{padding:7px;border:1px solid #ccc}img{max-width:100%}section{margin:40px 0}code{background:#eef2f5}h1,h2{color:#18476b}</style>',
            '<h1>Chromatin state analysis</h1><p>Post-segmentation analysis using verified ChromHMMTools comparisons.</p>',
            '<p>WT→KO differences are between-condition overlaps. They are distinct from spatial adjacency within one segmentation and from fitted HMM transition parameters. Summaries are descriptive; bins are not biological replicates.</p>',
            '<p>Tables below are previews. Full numerical tables are saved beside this report in the analysis directory. Images are embedded so this HTML remains viewable when copied.</p>']
    count=0
    if config['plots']:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    def figure(fig,name):
        fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=130)
        buf=io.BytesIO();fig.savefig(buf,format='png',dpi=110);plt.close(fig)
        return '<img alt="'+html.escape(name)+'" src="data:image/png;base64,'+base64.b64encode(buf.getvalue()).decode()+'">'
    def heat(matrix,rows,cols,title,filename,fixed=False):
        if not config['plots']:return ''
        fig,ax=plt.subplots(figsize=(max(5,min(14,len(cols)*.5)),max(3,min(14,len(rows)*.38))))
        data=np.array([[float('nan') if v is None else v for v in r] for r in matrix],dtype=float)
        image=ax.imshow(data,aspect='auto',vmin=0,vmax=1 if fixed else None,cmap='viridis')
        ax.set_xticks(range(len(cols)));ax.set_xticklabels(cols,rotation=70,ha='right',fontsize=8)
        ax.set_yticks(range(len(rows)));ax.set_yticklabels(rows,fontsize=8);ax.set_title(title)
        fig.colorbar(image,ax=ax)
        return figure(fig,filename)
    for name,model in models.items():
        states=sorted(model['states'],key=lambda s:model['annotations'][s]['order'])
        values=[model['values'][model['states'].index(s)] for s in states]
        labels=[s+': '+model['annotations'][s]['label'] for s in states]
        pieces.append('<h2>Model '+html.escape(name)+'</h2>')
        pieces.append(table([model['annotations'][s] for s in states],['state','label','order'],1000))
        pieces.append(heat(values,labels,model['marks'],'Emission probability (0–1)',name+'_emissions',True))
        for i,ctx in enumerate(model['contexts']):
            pieces.append('<h3>'+html.escape(ctx['name'])+'</h3><p>Imported '+html.escape(ctx['kind'])+' table; unit: '+html.escape(ctx['value_unit'])+'. Background summary rows: '+html.escape(','.join(ctx['background_rows']))+'</p>')
            pieces.append(heat(ctx['values'],ctx['states'],ctx['columns'],ctx['name']+' ('+ctx['value_unit']+')',name+'_context_'+str(i)))
    concordance=load(root/'concordance'/'concordance.json')
    if len(config['comparisons'])>1:
        pieces.append('<h2>Across-comparison overview</h2><p>Agreement is descriptive. Comparisons may share the same reference; they are not independent replicates.</p>')
        pieces.append(table(concordance,list(concordance[0]) if concordance else ['entity_id','window']))
    for c in config['comparisons']:
        data=load(root/('comparison_'+c['name'])/'analysis.json')
        metrics=data['metrics'];agg=data['aggregates'];upstream=data['upstream']
        pieces.append('<section><h2>'+html.escape(c['name'])+'</h2>')
        pieces.append('<p>'+str(sum(r['eligible'] for r in metrics))+' / '+str(len(metrics))+' entity windows meet shared-coverage requirements. Model correspondence: '+html.escape(upstream['state_mode'])+'.</p>')
        pieces.append('<p>Window bp totals are attributed to entities: overlapping promoters/transcripts may reuse genomic bases. Entity means give each eligible entity-window equal weight. Genomic unique coverage is reported by the upstream comparison.</p>')
        pieces.append('<p>Expression ID audit: '+html.escape(str(data['expression_id_counts']))+'</p>')
        pieces.append(table([dict(condition=k,**v) for k,v in upstream['conditions'].items()],['condition','covered_bp','unshared_bp','full_bins','transitions']))
        pieces.append(table(metrics,['entity_id','window','effective_bp','shared_bp','coverage_fraction','expression_group','changed_fraction']))
        ws,ms=data['wt_states'],data['mt_states']
        pieces.append('<p>For comparison heatmaps, WT states are rows and KO states are columns.</p>')
        pieces.append(heat(data['global_matrix'],ws,ms,'Global shared overlap (bp)',c['name']+'_global'))
        # One panel for every requested window/group, including non-significant and missing expression.
        keys=sorted({(r['window'],r['group_kind'],r['group']) for r in agg})
        for window,kind,group in keys:
            rows=[r for r in agg if (r['window'],r['group_kind'],r['group'])==(window,kind,group)]
            lookup={(r['wt_state'],r['mt_state']):r for r in rows};count+=1
            title=f'{window} / {kind}: {group}'
            pieces.append('<h3>'+html.escape(title)+'</h3>')
            for metric,label in [('bp_fraction','Attributed bp fraction'),('entity_mean_fraction','Equal-entity mean fraction'),('p_mt_given_wt','P(KO | WT)')]:
                pieces.append(heat([[lookup[(w,m)][metric] for m in ms] for w in ws],ws,ms,title+' — '+label,f'panel_{count}_{metric}',True))
        if config['plots']:
            for wi,window in enumerate(sorted({r['window'] for r in data['state_rows']})):
                for kind in ('all','expression','abundance'):
                    rows=[r for r in data['state_rows'] if r['window']==window and r['group_kind']==kind]
                    if not rows:continue
                    groups=sorted({r['group'] for r in rows});states=sorted({r['state'] for r in rows})
                    for metric in ('bp_fraction','entity_mean_fraction'):
                        fig,ax=plt.subplots(figsize=(max(7,len(groups)*1.8),4))
                        positions=np.arange(len(groups)*2);bottom=np.zeros(len(positions))
                        tracks=[(condition,state) for condition in ('WT','MT') for state in states] if upstream['state_mode']=='distinct' else [(None,state) for state in states]
                        for only_condition,state in tracks:
                            values=[next((r[metric] for r in rows if r['group']==g and r['condition']==condition and r['state']==state and (only_condition is None or condition==only_condition)),0)
                                    for g in groups for condition in ('WT','MT')]
                            model_name=c['target_model'] if only_condition=='MT' else c['reference_model']
                            color=models[model_name]['annotations'].get(state,{}).get('color','#808080')
                            ax.bar(positions,values,bottom=bottom,label=(only_condition+':' if only_condition else '')+state,color=color);bottom+=np.array(values)
                        ax.set_xticks(positions);ax.set_xticklabels([g+' '+cnd for g in groups for cnd in ('WT','MT')],rotation=55,ha='right',fontsize=8)
                        ax.set_ylim(0,1);ax.set_ylabel(metric);ax.set_title(window+' / '+kind);ax.legend(fontsize=7,bbox_to_anchor=(1.01,1))
                        pieces.append(figure(fig,f'{c["name"]}_states_{wi}_{kind}_{metric}'))
        pieces.append('<h3>Genomic context</h3><p>Disjoint bp annotation with promoter &gt; exon &gt; intron &gt; intergenic priority. Selection context includes all selected state-pair tiles, before expression/window filtering.</p>')
        pieces.append(table(data['genomic_context'],['selection','category','bp','fraction'],1000))
        if config['plots']:
            for si,selection in enumerate(sorted({r['selection'] for r in data['genomic_context']})):
                rows=[r for r in data['genomic_context'] if r['selection']==selection and r['bp']>0]
                if rows:
                    fig,ax=plt.subplots(figsize=(5,4));ax.pie([r['bp'] for r in rows],labels=[r['category'] for r in rows],autopct='%1.1f%%');ax.set_title(selection+' — covered bp')
                    pieces.append(figure(fig,c['name']+'_context_'+str(si)))
        pieces.append('<h3>Selected transitions</h3>')
        pieces.append(table(data['selected'],['selection','entity_id','gene_id','window','association','selected_bp','selected_fraction','expression_group']))
        eroot=root/('enrichment_'+c['name'])
        if eroot.exists():
            en=load(eroot/'enrichment.json')
            pieces.append('<h3>Gene-set enrichment</h3><p>One-sided hypergeometric over-representation, BH correction within comparison × selection × window. Backgrounds are eligible genes with membership in the supplied gene-set resource. Association is not evidence of causation.</p>')
            pieces.append(table(en['backgrounds'],['selection','window','eligible_genes','annotated_background','selected_annotated_genes','unmapped_background']))
            ranked=sorted(en['results'],key=lambda r:r['padj'])
            pieces.append('<p>Ranked tested terms are shown even when not significant; inspect adjusted p-values. Dot size represents overlapping genes.</p>')
            pieces.append(table(ranked,['selection','window','term_id','term_name','overlap_genes','fold_enrichment','padj']))
            if config['plots'] and ranked:
                top=ranked[:15];fig,ax=plt.subplots(figsize=(9,max(3,min(8,.35*len(top)+1.5))))
                points=ax.scatter([r['fold_enrichment'] for r in top],range(len(top)),s=[20+15*r['overlap_genes'] for r in top],
                                  c=[-math.log10(max(r['padj'],1e-300)) for r in top],cmap='viridis')
                ax.set_yticks(range(len(top)));ax.set_yticklabels([r['selection']+'/'+r['window']+': '+r['term_name'] for r in top],fontsize=8)
                ax.set_xlabel('Fold enrichment');fig.colorbar(points,ax=ax,label='−log10(BH adjusted p)')
                pieces.append(figure(fig,c['name']+'_enrichment'))
        pieces.append('</section>')
    pieces.append('<h2>Configuration and provenance</h2><pre>'+html.escape(__import__('json').dumps(load(root/'manifest.json'),indent=2))+'</pre></html>')
    (out/'report.html').write_text('\n'.join(pieces))


import math
