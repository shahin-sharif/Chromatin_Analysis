# Complete post-segmentation example

From the repository root:

~~~bash
python -m pip install -r requirements.txt
python examples/post_segmentation/run_example.py --outdir post-example
~~~

The destination must not exist. The runner generates tiny artificial input files,
executes ChromHMMTools compare, runs the new tool's deep preflight, runs the full
analysis, and checks known answers. It writes example-validation.json only after
those checks pass. No ChromHMM Java training is performed in this software test.
Use --no-plots for a faster numerical run.

The readable generate.py is the source for every fixture. It creates:

- A 1,000-base chrToy reference, with two states in WT and KO.
- Two-mark emissions, an explicit 1→E1/2→E2 label map and reviewed toy labels/colors.
- Small neighborhood and annotation-overlap tables.
- Four genes/transcripts with two exons each and both strands represented.
- A DE table with up, down, not-significant and not-tested examples.
- Two artificial gene sets with explicit memberships. These are not real GO terms.
- A complete configuration referencing newly generated upstream results.

Known answers:

| Measurement | Expected |
|---|---|
| Global WT E1→KO E1 overlap | 200 bp |
| Global WT E1→KO E2 overlap | 200 bp |
| Global WT E2→KO E2 overlap | 600 bp |
| Changed body fractions for g1, g2, g3, g4 | 0.5, 0.5, 0, 0 |
| Genes selected by any E1→E2 body support | g1.1 and g2.1 |
| Eligible annotated enrichment background | g1.1, g2.1, g3.1; g4 lacks a tested DE result |
| Toy process A over-representation p-value | 1/3 |
| BH adjusted p-value for that term | 2/3 (not significant) |

Open post-example/results/report/report.html after completion. The report embeds
its images; the same figures are available as PNG files. All numerical outputs
remain available in the stage folders. Expected p-values are software checks,
not evidence of biological enrichment.

To generate inputs without running the new pipeline:

~~~bash
python examples/post_segmentation/generate.py --outdir another-example
python ChromatinStateAnalysis.py check --config another-example/config.json
python ChromatinStateAnalysis.py run --config another-example/config.json
~~~

To exercise verified recovery:

~~~bash
python ChromatinStateAnalysis.py run --config another-example/config.json --resume
python ChromatinStateAnalysis.py report \
  --results another-example/results --outdir recovered-report
~~~

Resume rejects changed source/input/configuration or modified saved stage files.
The report-only command reads saved numerical data; original GTF/BED/DE files are
not required for that operation.
