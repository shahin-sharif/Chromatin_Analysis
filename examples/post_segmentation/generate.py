#!/usr/bin/env python3
"""Create an artificial complete post-segmentation example in a new directory."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[2]


def generate(out,plots=True):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'chrom.sizes').write_text('chrToy\t1000\n')
    (out/'wt.bed').write_text('chrToy\t0\t400\tE1\nchrToy\t400\t1000\tE2\n')
    (out/'ko.bed').write_text('chrToy\t0\t200\tE1\nchrToy\t200\t1000\tE2\n')
    (out/'emissions.tsv').write_text('state\tmarkA\tmarkB\n1\t0.9\t0.1\n2\t0.1\t0.9\n')
    (out/'state_ids.tsv').write_text('emission_state\tsegment_state\n1\tE1\n2\tE2\n')
    (out/'state_annotations.tsv').write_text('state\tlabel\tcolor\torder\nE1\tSynthetic promoter\t#bb4444\t1\nE2\tSynthetic repressed\t#5555aa\t2\n')
    (out/'neighborhood.tsv').write_text('state\t-200\t0\t200\n1\t1\t2\t1\n2\t1\t0.5\t1\n')
    (out/'overlap.tsv').write_text('state\tpromoter\texon\n1\t2\t1\n2\t0.5\t1\n')
    gtf=[]
    for i,(a,b,strand) in enumerate([(100,300,'+'),(300,500,'-'),(500,700,'+'),(700,900,'-')],1):
        attr=f'gene_id "g{i}.1"; gene_name "G{i}";'
        gtf.append(f'chrToy\tsynthetic\tgene\t{a+1}\t{b}\t.\t{strand}\t.\t{attr}\n')
        gtf.append(f'chrToy\tsynthetic\ttranscript\t{a+1}\t{b}\t.\t{strand}\t.\t{attr} transcript_id "t{i}.1";\n')
        for l,r in [(a,a+50),(b-50,b)]:gtf.append(f'chrToy\tsynthetic\texon\t{l+1}\t{r}\t.\t{strand}\t.\t{attr} transcript_id "t{i}.1";\n')
    (out/'genes.gtf').write_text(''.join(gtf))
    (out/'expression.tsv').write_text('entity_id\tlog2FoldChange\tpadj\tbaseMean\ng1.1\t2\t0.001\t600\ng2.1\t-1\t0.01\t150\ng3.1\t0.2\t0.5\t20\ng4.1\tNA\tNA\t0\n')
    (out/'terms.tsv').write_text('term_id\tterm_name\tgene_id\nTOY:1\tToy process A\tg1.1\nTOY:1\tToy process A\tg2.1\nTOY:2\tToy process B\tg3.1\nTOY:2\tToy process B\tg4.1\n')
    command=[sys.executable,str(ROOT/'ChromHMMTools.py'),'compare','--wt',str(out/'wt.bed'),'--mt',str(out/'ko.bed'),
             '--state-mode','shared','--bin-size','100','--write-change-bed','--outdir',str(out/'comparison')]
    subprocess.run(command,check=True)
    cfg=dict(schema_version=1,outdir='results',gtf='genes.gtf',chrom_sizes='chrom.sizes',entity_level='gene',plots=plots,
       windows=[dict(name='body',kind='body',mode='span'),dict(name='promoter',kind='window',anchor='tss',start=-100,end=100)],
       models={'shared':dict(emissions='emissions.tsv',state_id_map='state_ids.tsv',state_annotations='state_annotations.tsv',
            context_tables=[dict(name='tss_neighborhood',path='neighborhood.tsv',kind='neighborhood',value_unit='fold_enrichment'),
                            dict(name='annotation_overlap',path='overlap.tsv',kind='overlap',value_unit='fold_enrichment')])},
       comparisons=[dict(name='KO_vs_WT',directory='comparison',reference_model='shared',target_model='shared',expression='expression.tsv')],
       selections=[dict(name='E1_to_E2',pairs=[['E1','E2']],windows=['body'])],
       enrichment=dict(gene_sets='terms.tsv',source='synthetic known-answer memberships (not real GO)',release='1',membership='explicit_expanded',min_size=1,max_size=100))
    (out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    return out/'config.json'


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--no-plots',action='store_true')
    a=p.parse_args();print(generate(a.outdir,not a.no_plots))
