import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from chromatin_state import pipeline
from chromatin_state.inputs import configuration,model_data,expression_data
from chromatin_state.regions import build,TileIndex
from chromatin_state.analysis import summarize,hypergeom_tail,bh
from chromatin_state.common import load
spec=importlib.util.spec_from_file_location('post_generator',ROOT/'examples/post_segmentation/generate.py')
gen=importlib.util.module_from_spec(spec);spec.loader.exec_module(gen)


class PostSegmentation(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)/'case';self.config=gen.generate(self.base,False)

    def edit(self,change):
        cfg=load(self.config);change(cfg);self.config.write_text(json.dumps(cfg));return cfg

    def run_cli(self,*args,ok=True):
        p=subprocess.run([sys.executable,str(ROOT/'ChromatinStateAnalysis.py'),*map(str,args)],capture_output=True,text=True)
        self.assertEqual(p.returncode,0 if ok else 2,p.stdout+'\n'+p.stderr)
        return p

    def test_end_to_end_known_answers(self):
        self.run_cli('run','--config',self.config)
        root=self.base/'results';data=load(root/'comparison_KO_vs_WT/analysis.json')
        rows=[r for r in data['metrics'] if r['window']=='body']
        self.assertEqual([r['changed_fraction'] for r in rows],[.5,.5,0,0])
        self.assertEqual([r['expression_group'] for r in rows],['up','down','not_significant','not_tested'])
        self.assertEqual({r['gene_id'] for r in data['selected']},{'g1.1','g2.1'})
        allrows=[r for r in data['aggregates'] if r['window']=='body' and r['group_kind']=='all']
        self.assertEqual([r['bp'] for r in allrows],[100,200,0,500])
        self.assertAlmostEqual(sum(r['entity_mean_fraction'] for r in allrows),1)
        go=load(root/'enrichment_KO_vs_WT/enrichment.json')['results']
        self.assertAlmostEqual(go[0]['pvalue'],1/3);self.assertAlmostEqual(go[0]['padj'],2/3)
        self.assertEqual(go[0]['background_genes'],3)
        self.assertTrue((root/'COMPLETE.txt').is_file());self.assertFalse((root/'INCOMPLETE.txt').exists())
        self.assertIn('E1\tE1',(root/'comparison_KO_vs_WT/all_overlap_tiles.bed').read_text())

    def test_preflight_no_outputs(self):
        self.run_cli('check','--config',self.config)
        self.assertFalse((self.base/'results').exists())

    def test_changed_segmentation_rejected(self):
        with (self.base/'ko.bed').open('a') as f:f.write('chrToy\t1000\t1100\tE1\n')
        self.run_cli('check','--config',self.config,ok=False)

    def test_tampered_comparison_matrix_rejected(self):
        p=self.base/'comparison/overlap.bp.tsv';p.write_text(p.read_text().replace('200.0','201.0',1))
        # Exact integer formatting differs across NumPy versions; deliberately alter one numeric cell.
        lines=p.read_text().splitlines();row=lines[1].split('\t');row[1]=str(float(row[1])+1);lines[1]='\t'.join(row);p.write_text('\n'.join(lines)+'\n')
        self.run_cli('check','--config',self.config,ok=False)

    def test_state_mapping_required(self):
        self.edit(lambda c:c['models']['shared'].pop('state_id_map'))
        self.run_cli('check','--config',self.config,ok=False)

    def test_unknown_configuration_key(self):
        self.edit(lambda c:c.update(min_covergae=.5))
        self.run_cli('check','--config',self.config,ok=False)

    def test_reference_mismatch(self):
        (self.base/'chrom.sizes').write_text('chrToy\t500\n')
        self.run_cli('check','--config',self.config,ok=False)

    def test_strand_windows_and_exon_bodies(self):
        self.edit(lambda c:c.update(selections=[],enrichment=None,windows=[dict(name='p',kind='window',anchor='tss',start=-20,end=10),dict(name='e',kind='body',mode='exons')]))
        cfg=configuration(self.config);regions,_=build(cfg)
        self.assertEqual(regions[0]['blocks'],[(80,110)])
        self.assertEqual(regions[2]['blocks'],[(490,520)])
        self.assertEqual(regions[1]['blocks'],[(100,150),(250,300)])

    def test_short_and_boundary_clipped_windows(self):
        self.edit(lambda c:c.update(selections=[],enrichment=None,windows=[dict(name='p',kind='window',anchor='tss',start=-200,end=0)]))
        regions,_=build(configuration(self.config))
        self.assertEqual(regions[0]['requested_bp'],200);self.assertEqual(regions[0]['clipped_bp'],100)

    def test_transcript_level_exact_ids(self):
        self.edit(lambda c:c.update(entity_level='transcript'))
        p=self.base/'expression.tsv';p.write_text(p.read_text().replace('\ng','\nt'))
        self.run_cli('run','--config',self.config)
        r=load(self.base/'results/comparison_KO_vs_WT/analysis.json')['metrics'][0]
        self.assertEqual(r['entity_id'],'t1.1');self.assertEqual(r['gene_id'],'g1.1')

    def test_expression_unmatched_rejected(self):
        p=self.base/'expression.tsv';p.write_text(p.read_text().replace('\ng','\nunknown'))
        self.run_cli('check','--config',self.config,ok=False)

    def test_expression_duplicate_rejected(self):
        p=self.base/'expression.tsv';p.write_text(p.read_text()+'g1.1\t1\t0.01\t4\n')
        self.run_cli('check','--config',self.config,ok=False)

    def test_explicit_expression_id_map(self):
        p=self.base/'expression.tsv';p.write_text(p.read_text().replace('\ng','\nx'))
        (self.base/'ids.tsv').write_text('expression_id\tentity_id\n'+''.join(f'x{i}.1\tg{i}.1\n' for i in range(1,5)))
        self.edit(lambda c:c.update(expression=dict(id_map='ids.tsv')))
        self.run_cli('check','--config',self.config)

    def test_ambiguous_id_map_rejected(self):
        (self.base/'ids.tsv').write_text('expression_id\tentity_id\na\tg1.1\nb\tg1.1\n')
        self.edit(lambda c:c.update(expression=dict(id_map='ids.tsv')))
        self.run_cli('check','--config',self.config,ok=False)

    def test_resume_reuses_verified_outputs(self):
        self.run_cli('run','--config',self.config)
        file=self.base/'results/comparison_KO_vs_WT/analysis.json';before=file.stat().st_mtime_ns
        result=self.run_cli('run','--config',self.config,'--resume')
        self.assertIn('Reusing verified stage',result.stdout)
        self.assertEqual(before,file.stat().st_mtime_ns)

    def test_resume_detects_modified_outputs(self):
        self.run_cli('run','--config',self.config)
        (self.base/'results/comparison_KO_vs_WT/entity_metrics.tsv').write_text('modified')
        self.run_cli('run','--config',self.config,'--resume',ok=False)

    def test_resume_detects_config_change(self):
        self.run_cli('run','--config',self.config)
        self.edit(lambda c:c.update(min_coverage=.5))
        self.run_cli('run','--config',self.config,'--resume',ok=False)

    def test_refuses_existing_output(self):
        (self.base/'results').mkdir();(self.base/'results/keep').write_text('important')
        self.run_cli('run','--config',self.config,ok=False)
        self.assertEqual((self.base/'results/keep').read_text(),'important')

    def test_report_recovery_without_original_inputs(self):
        self.run_cli('run','--config',self.config)
        (self.base/'wt.bed').unlink();(self.base/'expression.tsv').unlink()
        self.run_cli('report','--results',self.base/'results','--outdir',self.base/'recovered')
        self.assertTrue((self.base/'recovered/report.html').exists())

    def test_failed_report_resume_without_reanalysis(self):
        with patch('chromatin_state.pipeline.report',side_effect=RuntimeError('simulated plotting failure')):
            with self.assertRaises(RuntimeError):pipeline.run(self.config)
        file=self.base/'results/comparison_KO_vs_WT/analysis.json';before=file.stat().st_mtime_ns
        self.assertFalse((self.base/'results/COMPLETE.txt').exists())
        self.run_cli('run','--config',self.config,'--resume')
        self.assertEqual(before,file.stat().st_mtime_ns)

    def test_multiple_comparisons(self):
        def change(c):
            other=copy.deepcopy(c['comparisons'][0]);other['name']='Second_vs_WT';c['comparisons'].append(other)
        self.edit(change);self.run_cli('run','--config',self.config)
        rows=load(self.base/'results/concordance/concordance.json')
        self.assertTrue(rows[0]['consistent_expression_direction'])
        self.assertIsNone(rows[0]['consistent_dominant_pair'])
        self.assertTrue(rows[1]['consistent_dominant_pair'])

    def test_genomic_context_base_conservation(self):
        self.run_cli('run','--config',self.config)
        rows=load(self.base/'results/comparison_KO_vs_WT/analysis.json')['genomic_context']
        self.assertEqual(sum(r['bp'] for r in rows if r['selection']=='all_shared'),1000)
        self.assertEqual(sum(r['bp'] for r in rows if r['selection']=='E1_to_E2'),200)

    def test_gene_sets_no_matching_ids_rejected(self):
        p=self.base/'terms.tsv';p.write_text(p.read_text().replace('\tg','\tunmapped'))
        self.run_cli('check','--config',self.config,ok=False)

    def test_distinct_models_have_no_retention_claim(self):
        cmd=[sys.executable,str(ROOT/'ChromHMMTools.py'),'compare','--wt',str(self.base/'wt.bed'),'--mt',str(self.base/'ko.bed'),
             '--state-mode','distinct','--outdir',str(self.base/'distinct')]
        subprocess.run(cmd,check=True,capture_output=True)
        def change(c):
            c['models']['other']=copy.deepcopy(c['models']['shared'])
            c['comparisons'][0].update(directory='distinct',target_model='other')
        self.edit(change);self.run_cli('run','--config',self.config)
        self.assertTrue(all(r['changed_fraction'] is None for r in load(self.base/'results/comparison_KO_vs_WT/analysis.json')['metrics']))

    def test_model_unobserved_state_retained(self):
        with (self.base/'emissions.tsv').open('a') as f:f.write('3\t0.5\t0.5\n')
        with (self.base/'state_ids.tsv').open('a') as f:f.write('3\tE3\n')
        self.edit(lambda c:(c['models']['shared'].pop('state_annotations'),c['models']['shared'].pop('context_tables')))
        self.run_cli('run','--config',self.config)
        data=load(self.base/'results/comparison_KO_vs_WT/analysis.json')
        self.assertEqual(data['wt_states'],['E1','E2','E3'])
        self.assertTrue(all(r['p_mt_given_wt'] is None for r in data['aggregates'] if r['wt_state']=='E3'))

    def test_invalid_emission_probabilities(self):
        p=self.base/'emissions.tsv';p.write_text(p.read_text().replace('0.9','1.9'))
        self.run_cli('check','--config',self.config,ok=False)

    def test_custom_nearest_tss_is_explicit(self):
        (self.base/'custom.bed').write_text('chrToy\t0\t10\tdistal\n')
        self.edit(lambda c:c.update(custom_regions=[dict(name='distal',path='custom.bed',nearest_tss=True)]))
        regions,_=build(configuration(self.config));r=regions[-1]
        self.assertEqual(r['gene_id'],'g1.1');self.assertEqual(r['association'],'nearest_tss_proximity')

    def test_explicit_nearest_expression_join(self):
        (self.base/'custom.bed').write_text('chrToy\t0\t10\tdistal\n')
        self.edit(lambda c:c.update(custom_regions=[dict(name='distal',path='custom.bed',nearest_tss=True,join_nearest_expression=True)]))
        self.run_cli('run','--config',self.config)
        row=load(self.base/'results/comparison_KO_vs_WT/analysis.json')['metrics'][-1]
        self.assertEqual(row['expression_group'],'up')
        self.assertEqual(row['expression_entity_id'],'g1.1')
        self.assertEqual(row['association'],'nearest_tss_proximity')

    def test_relocated_upstream_inputs(self):
        (self.base/'wt.bed').rename(self.base/'relocated.bed')
        self.edit(lambda c:c['comparisons'][0].update(input_overrides=dict(wt='relocated.bed')))
        self.run_cli('check','--config',self.config)

    def test_explicit_mapped_models(self):
        (self.base/'ko_other.bed').write_text((self.base/'ko.bed').read_text().replace('E1','A').replace('E2','B'))
        (self.base/'mapping.tsv').write_text('A\tE1\nB\tE2\n')
        subprocess.run([sys.executable,str(ROOT/'ChromHMMTools.py'),'compare','--wt',str(self.base/'wt.bed'),
                        '--mt',str(self.base/'ko_other.bed'),'--state-mode','mapped','--state-map',str(self.base/'mapping.tsv'),
                        '--outdir',str(self.base/'mapped')],check=True,capture_output=True)
        (self.base/'other_ids.tsv').write_text('emission_state\tsegment_state\n1\tA\n2\tB\n')
        def change(c):
            c['models']['other']=dict(emissions='emissions.tsv',state_id_map='other_ids.tsv')
            c['comparisons'][0].update(directory='mapped',target_model='other')
        self.edit(change);self.run_cli('run','--config',self.config)
        rows=load(self.base/'results/comparison_KO_vs_WT/analysis.json')['metrics']
        self.assertEqual(rows[0]['changed_fraction'],.5)

    def test_masked_shared_coverage(self):
        (self.base/'mask.bed').write_text('chrToy\t150\t300\n')
        subprocess.run([sys.executable,str(ROOT/'ChromHMMTools.py'),'compare','--wt',str(self.base/'wt.bed'),
                        '--mt',str(self.base/'ko.bed'),'--state-mode','shared','--exclude-bed',str(self.base/'mask.bed'),
                        '--outdir',str(self.base/'masked')],check=True,capture_output=True)
        self.edit(lambda c:c['comparisons'][0].update(directory='masked'))
        self.run_cli('run','--config',self.config)
        row=load(self.base/'results/comparison_KO_vs_WT/analysis.json')['metrics'][0]
        self.assertEqual(row['shared_bp'],50);self.assertFalse(row['eligible'])

    def test_coverage_and_entity_weighting(self):
        regions=[dict(region_id='r1',entity_id='a',gene_id='a',window='body',association='gene',chrom='c',blocks=[(0,100)],effective_bp=100,requested_bp=100,clipped_bp=0),
                 dict(region_id='r2',entity_id='b',gene_id='b',window='body',association='gene',chrom='c',blocks=[(100,400)],effective_bp=300,requested_bp=300,clipped_bp=0)]
        comp=dict(tiles=[('c',0,100,'E1','E2'),('c',100,400,'E1','E1')],info=dict(state_mode='shared'),wt_states=['E1'],mt_states=['E1','E2'])
        cfg=dict(min_coverage=.8,selections=[])
        rows,counts,agg,states,sel=summarize(regions,comp,{},cfg,False)
        changed=next(r for r in agg if r['mt_state']=='E2' and r['group_kind']=='all')
        self.assertEqual(changed['bp_fraction'],.25);self.assertEqual(changed['entity_mean_fraction'],.5)
        comp['tiles']=[('c',0,20,'E1','E2')]
        rows,*_=summarize(regions,comp,{},cfg,False)
        self.assertEqual(rows[0]['status'],'insufficient_coverage');self.assertFalse(rows[1]['eligible'])

    def test_hypergeometric_against_integer_oracle(self):
        for N in range(2,15):
            for K in range(1,N):
                for n in range(1,N):
                    for k in range(0,min(K,n)+2):
                        exact=sum(math.comb(K,i)*math.comb(N-K,n-i) for i in range(k,min(K,n)+1) if 0<=n-i<=N-K)/math.comb(N,n)
                        self.assertAlmostEqual(hypergeom_tail(k,N,K,n),exact,places=11)
        self.assertEqual(bh([.01,.04,.03]),[.03,.04,.04])

    def test_randomized_interval_oracle(self):
        import random
        random.seed(81)
        tiles=[('c',i,i+10,'E1','E2' if i%30 else 'E1') for i in range(0,200,10)]
        index=TileIndex(tiles)
        for _ in range(80):
            a=random.randrange(200);b=random.randrange(a+1,230)
            hits=list(index.intersect(dict(chrom='c',blocks=[(a,b)])))
            expected=sum(1 for pos in range(a,b) if 0<=pos<200)
            self.assertEqual(sum(y-x for x,y,w,m in hits),expected)


if __name__=='__main__':unittest.main()
