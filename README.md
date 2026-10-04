# Chromatin_Analysis

Tools for inspecting and comparing chromatin-state segmentations, with explicit
state correspondence, reproducible examples and interpretable reports.

The comparison tool is **[ChromHMMTools.py](ChromHMMTools.py)**. It consolidates the
comparison and reporting functionality previously developed from
CompareChromHMM_new.py and ReportChromHMM.py into three subcommands:

| Command | Purpose |
|---|---|
| compare | Calculate exact overlaps, annotation-stratified counts and empirical spatial transitions |
| report | Build a portable HTML report from saved comparison outputs |
| run | Perform comparison and reporting together |

The new **[ChromatinStateAnalysis.py](ChromatinStateAnalysis.py)** consumes the
modeling inputs and verified ChromHMMTools outputs, then connects explicit
gene/transcript windows to expression, selected state pairs and gene-set enrichment.
Its [complete guide](docs/ChromatinStateAnalysis.md) explains configuration,
coordinate rules, backgrounds, recovery and every output.

The repository does not train ChromHMM models or perform statistical differential
chromatin testing. Gene-set enrichment is a separate over-representation analysis
with an explicit universe and multiple-testing correction.

## Download and install

Run these commands in Bash. GitHub authentication is required while the
repository is private.

~~~bash
git clone https://github.com/shahin-sharif/Chromatin_Analysis.git
cd Chromatin_Analysis
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
~~~

Python 3.8+ and NumPy are required. Matplotlib generates report images.
The tool also supports reports without plots through --no-plots.

## Start with a tiny example

~~~bash
python examples/run_example.py --outdir sample-results
~~~

The example verifies:

- Overlap counts of 200, 200, 0 and 400 bp.
- 800 shared bp and 75% matching-state retention.
- Annotation allocation of 300 genic, 200 promoter and 300 intergenic bp.
- The changed interval at chr1:200–400, plus embedded report images.

All example inputs are synthetic. They are not hg38 or experimental data.
Use a new output directory for each run. Open sample-results/report.html in
your browser; it is self-contained.

A [pre-generated synthetic report](examples/example-report.html) is included.
Download it or open it locally; GitHub's file viewer shows HTML source.

## Continue from segmentation comparison to biological interpretation

~~~bash
# Full known-answer demonstration, including upstream comparison and plots.
python examples/post_segmentation/run_example.py --outdir post-example

# For a configured real analysis:
python ChromatinStateAnalysis.py check --config analysis.json
python ChromatinStateAnalysis.py run --config analysis.json
~~~

The new tool verifies original upstream input hashes and reproduces saved overlap
counts before using them. It then measures exact promoter/body/TES or custom
windows, joins a supplied DE table with an ID audit, produces base-pair and
equal-entity summaries, selects named transitions and optionally tests local,
versioned GO/gene-set memberships. A self-contained HTML report and all numerical
outputs are saved. Resume verifies source/input/configuration and stage checksums;
report-only recovery does not require the original input files.

The [example guide](examples/post_segmentation/README.md) includes expected
numerical results. The [implementation report](docs/ChromatinStateAnalysis_ImplementationReport.md)
maps the historical R analysis to the new workflow and states validation limits.

## Compare your segmentations

~~~bash
python ChromHMMTools.py run \
  --wt WT_segments.bed \
  --mt MT_segments.bed \
  --state-mode shared \
  --bin-size 200 \
  --write-change-bed \
  --outdir comparison01
~~~

Use shared only if both files use the same state definitions, such as
segmentations from one fitted model. For independently trained models, use
distinct, or mapped with an explicit correspondence justified from the models.
Identical state numbers alone do not establish equivalence.

Inputs must use the same genome assembly and compatible chromosome names.
BED coordinates are zero-based, half-open. The bin size controls empirical
adjacent-bin summaries and should match the segmentation setup.

Read the [complete tool guide](docs/ChromHMMTools.md) for masks, annotations,
state mapping, output definitions and separate compare/report commands.

## Outputs and interpretation

Outputs include exact base-pair overlap matrices, normalized overlaps,
annotation-stratified overlaps, observed adjacent-bin transitions, optional
changed intervals, provenance metadata and an HTML report.

The observed spatial transitions are **not** ChromHMM's fitted HMM transition
parameters. Rankings and enrichments are descriptive, not significance tests.
Reports contain run paths and parameters; inspect them before sharing.

## Validation

~~~bash
python -m unittest discover -s tests -v
python examples/run_example.py --outdir sample-results-validation
~~~

The original twenty regression tests cover interval/bin oracles, gaps, masks, annotations,
state mapping, unsupported rows, rectangular matrices, report consistency,
escaping and file safety. The example runner also checks report images.

The new post-segmentation tests additionally cover known-answer windows,
expression joins, state correspondence, upstream integrity, coverage weighting,
enrichment against an independent integer oracle, resume and report recovery.
See the implementation report for the final test count and validation evidence.

No representative real ChromHMM segmentation pair has been validated yet.
Synthetic checks do not establish biological suitability for every analysis.

## Migration from ToolKit

ChromHMMTools.py, its tests, guide and synthetic fixtures moved here from
[ToolKit](https://github.com/shahin-sharif/ToolKit), source commit
9dde4ddbfa7bb75926a765fc1bd6aa5d5625b1f9, on 3 October 2026.
The implementation and command-line interface are unchanged in this migration.
Historical source remains available in ToolKit's Git history.

ToolKit remains the home for general-purpose interval, alignment, matrix and
bedGraph utilities. Use this repository for future chromatin-specific modules.
Original scripts on the source drive have not been changed.

Maintained by Shahin Behrouz Sharif.
