#!/usr/bin/env python3
"""Exercise modeling inputs -> ChromHMMTools -> post-segmentation analysis."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from generate import generate,ROOT

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--outdir',required=True)
p.add_argument('--no-plots',action='store_true')
a=p.parse_args()
config=generate(a.outdir,not a.no_plots)
for command in ('check','run'):
    subprocess.run([sys.executable,str(ROOT/'ChromatinStateAnalysis.py'),command,'--config',str(config)],check=True)
root=config.parent/'results'
data=json.loads((root/'comparison_KO_vs_WT/analysis.json').read_text())
rows=[r for r in data['metrics'] if r['window']=='body']
assert [r['changed_fraction'] for r in rows]==[.5,.5,0,0]
assert {r['gene_id'] for r in data['selected']}=={'g1.1','g2.1'}
go=json.loads((root/'enrichment_KO_vs_WT/enrichment.json').read_text())['results']
assert abs(go[0]['pvalue']-1/3)<1e-12 and abs(go[0]['padj']-2/3)<1e-12
assert (root/'COMPLETE.txt').is_file()
if not a.no_plots:
    assert 'data:image/png;base64,' in (root/'report/report.html').read_text()
(root/'example-validation.json').write_text(json.dumps(dict(status='PASS',synthetic=True,
    changed_body_fractions=[.5,.5,0,0],selected_genes=['g1.1','g2.1'],enrichment_pvalue=go[0]['pvalue'],
    enrichment_padj=go[0]['padj'],plots=not a.no_plots),indent=2)+'\n')
print('PASS: windows, unchanged territory, expression groups, selected genes, explicit GO background and complete report')
