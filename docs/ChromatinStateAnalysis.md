# ChromatinStateAnalysis.py

A reproducible post-segmentation workflow connecting chromatin states to genomic
features, transcript/gene expression, selected state changes and functional gene
sets. The workflow accepts ChromHMM modeling files and verified results from
ChromHMMTools.py. It does not train models, call differential expression or treat
genomic bins as biological replicates.

## Where it belongs in the analysis

~~~text
Histone-mark data
       |
       v
ChromHMM modeling
  segmentations + emissions + optional context enrichment tables
       |
       v
ChromHMMTools.py compare/run
  exact WT/KO overlaps + masks/state correspondence + RUNINFO.json
       |
       v
ChromatinStateAnalysis.py
  explicit gene/transcript windows
  expression-stratified state changes and composition
  selected transitions + genomic annotation + enrichment
  reproducible tables, figures and report
~~~

ChromHMMTools supplies the trusted global comparison. Its overlap matrices use
exact shared base pairs, with row/column/global normalization. It also describes
adjacent-bin spatial transitions within each condition, exports optional changed
intervals, and builds a report. Those spatial transitions are not ChromHMM's fitted
HMM transition parameters.

The new tool verifies upstream input hashes and reconstructs coordinate-level
**all-pair overlaps, including unchanged regions**, using the existing comparison
engine. It requires those overlaps to reproduce the saved global count matrix
exactly. A changed-only BED or global matrix alone is insufficient for per-gene
window denominators. Original segmentation files must therefore remain available,
or be relocated using hash-verified input overrides.

The existing ChromHMMTools implementation and interface are unchanged.

## Install and run the complete example

Clone/download the complete Chromatin_Analysis repository. Keep ChromHMMTools.py,
ChromatinStateAnalysis.py and the chromatin_state/ directory together.

~~~bash
cd Chromatin_Analysis
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

# Creates toy modeling inputs, runs ChromHMMTools, checks and runs the new tool.
# Choose a directory that does not already exist.
python examples/post_segmentation/run_example.py --outdir post-example

# Open post-example/results/report/report.html in a web browser.
# Inspect the numerical output in Bash:
cat post-example/results/comparison_KO_vs_WT/entity_metrics.tsv
cat post-example/results/enrichment_KO_vs_WT/go_enrichment.tsv
~~~

All fixture coordinates, genes, histone marks and gene sets are artificial. The
example does not run ChromHMM model training; it creates small model-like inputs.
The synthetic enrichment terms are not real GO terms. Use --no-plots on the
example runner for a faster numerical check.

Python 3.8+ and NumPy are required. Matplotlib is needed when plots are enabled.
Segmentations, annotations and summary tables are held in memory; many
transcripts/windows and full plotting can produce large outputs. Genome-scale
resource usage has not yet been benchmarked.

This version implements the enrichment calculation in Python and does not require
R, ChIPseeker or a live web enrichment service. Local gene-set memberships are an
explicit, versioned input rather than an automatically changing online resource.

## Inputs

| Input | Requirement and purpose |
|---|---|
| ChromHMMTools output directory | Required per comparison: RUNINFO.json and overlap.bp.tsv from version 1.0.0, plus original hashed inputs |
| ChromHMM emission table | Required per model, numeric probabilities in [0,1], original row IDs and mark columns |
| Matching GTF | Required for exact gene/transcript identity, boundaries and exon unions |
| Chromosome sizes | Required to validate shared territory and clip constructed windows |
| Differential-expression table | Optional per comparison; required for expression-linked analyses |
| State annotation table | Optional reviewed labels, colors and plotting order |
| Emission-to-segmentation ID map | Required when labels differ, e.g. emission row 1 versus BED state E1 |
| Context enrichment tables | Optional ChromHMM overlap/neighborhood tables, with explicit units |
| Custom regions | Optional named BED4–BED6 files |
| Gene-set membership table | Required when enrichment is enabled; exact annotation gene IDs, versioned source |

A single comparison can omit expression and/or enrichment, retaining state and
window analysis. Multiple comparisons can represent separate KO clones, provided
separate segmentations actually exist. A pooled KO segmentation cannot recover
clone-specific chromatin results.

## Configure a real analysis

All relative paths in the configuration resolve relative to the JSON file, not
the shell's current directory. This is an example configuration; replace the
paths, model IDs and selected state pairs with those from the intended analysis.

~~~json
{
  "schema_version": 1,
  "outdir": "results/post01",
  "gtf": "references/annotation.gtf.gz",
  "chrom_sizes": "references/chrom.sizes",
  "entity_level": "gene",
  "min_coverage": 0.8,
  "plots": true,
  "annotation_promoter": [-1000, 200],
  "models": {
    "joint13": {
      "emissions": "model/emissions_13.txt",
      "state_id_map": "model/emission_ids.tsv",
      "state_annotations": "model/state_annotations.tsv"
    }
  },
  "comparisons": [
    {
      "name": "B113_vs_WT",
      "directory": "chromhmm_comparisons/B113_vs_WT",
      "reference_model": "joint13",
      "target_model": "joint13",
      "expression": "expression/B113_vs_WT.tsv"
    },
    {
      "name": "B142_vs_WT",
      "directory": "chromhmm_comparisons/B142_vs_WT",
      "reference_model": "joint13",
      "target_model": "joint13",
      "expression": "expression/B142_vs_WT.tsv"
    }
  ],
  "expression": {
    "id_column": "gene_id",
    "log2fc_column": "log2FoldChange",
    "padj_column": "padj",
    "abundance_column": "baseMean",
    "fdr": 0.05,
    "min_abs_log2fc": 0.5
  },
  "abundance_breaks": [100, 500],
  "windows": [
    {"name": "tss_core", "kind": "window", "anchor": "tss", "start": -25, "end": 25},
    {"name": "upstream", "kind": "window", "anchor": "tss", "start": -1000, "end": -25},
    {"name": "downstream", "kind": "window", "anchor": "tss", "start": 25, "end": 1000},
    {"name": "broad_promoter", "kind": "window", "anchor": "tss", "start": -3000, "end": 3000},
    {"name": "body", "kind": "body", "mode": "span"},
    {"name": "tes", "kind": "window", "anchor": "tes", "start": -500, "end": 500}
  ],
  "selections": [
    {"name": "E10_to_E8_E9", "pairs": [["E10", "E8"], ["E10", "E9"]],
     "windows": ["upstream", "tss_core", "downstream"], "min_bp": 1, "min_fraction": 0.1}
  ],
  "enrichment": {
    "gene_sets": "references/go_memberships.tsv.gz",
    "source": "Describe the annotation database and identifier mapping used",
    "release": "Record the actual database release",
    "membership": "explicit_expanded",
    "min_size": 5,
    "max_size": 500
  }
}
~~~

The thresholds shown are configurable analysis choices. Use one comparison if
only one WT/KO pair is available. A model name denotes a declared fitted model,
not an arbitrary plotting group. In upstream shared mode both sides must reference
one model entry. Distinct models require upstream distinct or mapped mode.

Unknown configuration keys are errors, so misspelled options cannot silently use
defaults. Optional sections can be omitted. Omit enrichment until a compatible
membership table is available; the tool does not silently skip a requested module.

### State IDs, labels and context tables

An emission ID mapping file is a tab-separated table:

~~~text
emission_state segment_state
1              E1
2              E2
~~~

Use actual tabs, not the alignment spaces shown above. Mapping must cover every
emission row exactly once. It does not establish correspondence between separately
trained models; that is the upstream state-mode/mapping decision.

State annotations require columns state, label, color, order. Each state appears
once, colors use #RRGGBB and order values are unique integers. Unobserved model
states remain on summary axes with zero support; unsupported conditional
probabilities are NA. IDs are not inferred from row positions.

Each model can have context_tables entries, for example:

~~~json
{"name":"WT_TSS", "path":"WT_13_RefSeqTSS_neighborhood.txt",
 "kind":"neighborhood", "value_unit":"fold_enrichment"}
~~~

Accepted units: fold_enrichment, fraction, percent, count. Plots show raw values
with their stated units; log10(0) is never relabeled as zero. Optional columns lists
select exact column names. For mixed-unit overlap files, exclude Genome % from a
fold-enrichment panel and import it separately as percent. Recognized background
rows are retained and identified; unexpected state rows are rejected. Context
values are imported, not recalculated or used to automatically assign state names.

### Coordinate and annotation conventions

BED coordinates are zero-based, half-open. GTF coordinates are converted from
one-based inclusive. For a positive-strand entity [s,e), the TSS boundary is s;
for a negative-strand entity it is e. A window [a,b) becomes [TSS+a,TSS+b) on +
and [TSS-b,TSS-a) on -. TES uses the opposite boundary. Thus [-25,25) is 50 bp
on either strand before clipping.

This explicit boundary convention differs from the existing ChromHMMTools GTF
promoter option, which includes the TSS base and has up+down+1 bases. The new
window tables define their own documented measurements; upstream annotation
matrices are not silently substituted for them.

entity_level is gene or transcript. Declared boundaries are preferred; without
a declaration, the span of available features is used and recorded. Duplicate
or inconsistent locus identities are errors. The optional entity_ids file contains
exact IDs, one per line, for a user-selected subset. No automatic canonical or
longest-transcript choice is made. Version suffixes and PAR_Y identifiers remain
unchanged.

Body mode span includes introns and the whole annotated span; exons measures the
union of annotated exons. A missing exon model yields an unavailable window.
Windows may overlap: broad promoter and body are not disjoint analysis groups.
The region export is BED6 with one row per block, and the JSON retains the
entity/window/block relationship.

Separate genome-context annotation divides shared bases into disjoint categories,
with promoter > exon > intron > genic_unresolved > intergenic priority.
annotation_promoter sets its TSS-boundary offsets. Genes without exon annotation
are genic_unresolved, not presumed intronic. This context annotation uses all GTF
genes even when the entity analysis selects a subset. It is a reproducible feature
partition, not an exact emulation of ChIPseeker's nearest-transcript algorithm.

### Custom regions and proximity assignments

~~~json
"custom_regions": [
  {"name":"enhancers", "path":"enhancers.bed", "nearest_tss":true,
   "join_nearest_expression":true}
]
~~~

Each BED4–BED6 row needs a unique name. BED12 is deliberately rejected for this
custom-region interface; exon blocks come from GTF mode. nearest_tss defaults to
false. When enabled, it assigns a gene by distance from the interval to a TSS
boundary, with no maximum distance cutoff; distance ties are ambiguous and are
not assigned. This is recorded as proximity, not a demonstrated target-gene link.

Expression is not joined to nearest genes unless join_nearest_expression is
explicitly true; that option requires gene-level analysis. The expression entity
ID and association method remain in every output. Custom regions retain their
own region identities, so summaries over these rows are equal-region summaries,
not necessarily equal-gene summaries.

### Expression input and identifier mapping

Supply the full tested differential-expression table when possible, rather than
only significant results. The tool reads saved results and does not fit DESeq2. Positive log2FoldChange
must mean higher expression in the target/KO relative to the reference/WT for
that comparison; column names alone cannot verify contrast direction.
Default columns are entity_id, log2FoldChange, padj and baseMean; names can be
configured. abundance_column may be null to omit abundance measurements.

Groups are up, down, not_significant, not_tested, unmatched, or not_supplied.
Up/down require finite fold change and adjusted p-value, padj <= fdr, and absolute
fold change >= min_abs_log2fc. Zero fold change is not called up/down. Missing
statistics remain not_tested; missing IDs remain unmatched. Not significant is
not evidence of biological equivalence.

Use expression.id_map for an explicit two-column TSV: expression_id, entity_id.
Mappings must be one-to-one; ambiguous aggregation of multiple transcripts into
a gene is rejected. An optional gene_column checks gene identity against the GTF.
The expression_id_audit.tsv includes unmatched expression IDs and annotation IDs
without results. Zero matching IDs is an error. No symbol-based merging or silent
version stripping occurs.

Abundance bins use the declared measurement and configurable breaks. Default
[100,500] yields [0,100), [100,500), [500,infinity). baseMean is not TPM, and these
bins are descriptive and contrast/input dependent.

## Measurements and interpretation

The eligible territory is the intersection of the two retained segmentations,
after the upstream include/exclude masks and state mapping. Missing/unshared
coverage is never interpreted as quiescence or unchanged state.

For every entity-window:

- effective_bp is its clipped genomic/block length;
- shared_bp is the number of bases covered by both retained segmentations;
- coverage_fraction is shared_bp/effective_bp;
- eligibility requires nonzero shared coverage and coverage_fraction >= min_coverage;
- pair fractions use shared_bp as their denominator;
- changed fraction is available only in shared/mapped state mode.

Excluded mask territory reduces shared coverage; it is not removed from the
original window's effective length. Clipping and unavailable status are audited.
Low-coverage rows remain available in the per-entity outputs but do not enter
group summaries, selections or enrichment backgrounds.

Each group matrix reports attributed bp, attributed bp fraction, equal-entity
mean fraction and P(KO state | WT state). Attribution may reuse genomic bases
across overlapping genes/transcripts, so totals are explicitly not unique-genome
coverage. Equal-entity fractions first divide by each entity's shared coverage
and then average. This prevents longer entities dominating that summary.
Within any one entity-window, blocks and overlap bases are counted once.

All, expression and abundance groups are separate views of the same data and
must not be added together. No statistical differential chromatin test is run.
For multiple comparisons, concordance.tsv keeps per-comparison values side by
side and flags consistent up/down expression direction. Dominant-pair agreement is
reported only for eligible, untied results with shared/mapped states tied to the
same reference model; other cases are NA. This does not test complete composition
equivalence. Shared WT inputs do not
make those comparisons independent biological replicates.

## Transition selections and GO/gene-set enrichment

A selection specifies state pairs and optional windows, expression_groups,
min_bp and min_fraction. Pair support is summed over the selected pairs within
each entity-window. At least min_bp shared bases and the minimum shared-coverage
fraction must satisfy the selection. State labels refer to the upstream comparison
axes, including mapped labels when state-mode is mapped.

selected_entities.tsv applies all these filters. selected_transitions.bed and
its genome-context chart instead show **all genomic tiles matching the chosen
state pairs**, before expression/window filters. This distinction is intentional
and stated in the report.

Membership TSV columns are term_id, term_name, gene_id. Genes must match GTF gene
IDs exactly, even in transcript-level analysis. Duplicate memberships are collapsed.
For GO, supply an explicitly expanded annotation including the desired ancestor
memberships, consistent ontology/evidence filters, and a recorded database release.
This release does not download GO or propagate an ontology. The same engine also
accepts other explicitly defined gene sets; the historical output filename
here is go_enrichment.tsv but term IDs need not be GO IDs.

For each comparison × selection × window:

1. Form the eligible background from covered entities matching the selection's
   window/expression criteria. If expression was supplied, retain tested entities
   (up/down/not_significant), excluding unmatched and not_tested rows.
2. Collapse to unique GTF gene IDs, then intersect with genes represented in the
   supplied membership resource. Export both losses and the exact gene universe.
3. Collect selected genes within that same universe. A gene is selected if at
   least one eligible entity meets the selection; this matters in transcript mode.
4. Test over-representation using the one-sided hypergeometric tail.
5. Apply Benjamini–Hochberg correction across every eligible term in that family,
   including terms with zero selected hits. Term-size limits use background counts.

Adjusted p-values are per stated family, not a single correction across all
windows/comparisons. Empty selections/backgrounds produce audited empty result
tables. A resource with no annotation gene-ID matches fails preflight. Longer
or more transcript-rich genes can still differ in opportunity for selection;
GO results are descriptive association follow-ups with their stated background,
not proof of a causal mechanism or a bias-adjusted genomic enrichment model.

## Run, resume and rebuild a report

~~~bash
# Deep preflight reads inputs, verifies hashes and checks exact upstream counts.
python ChromatinStateAnalysis.py check --config analysis.json

# Start in a new output directory specified by the JSON.
python ChromatinStateAnalysis.py run --config analysis.json

# Recover after a failed later stage, preserving validated earlier stages.
python ChromatinStateAnalysis.py run --config analysis.json --resume

# Rebuild a report using saved numerical data; original inputs are not required.
python ChromatinStateAnalysis.py report \
  --results results/post01 --outdir results/post01-report-rebuilt
~~~

Resume compares configuration, input SHA-256 hashes and source-code hashes, then
checks every reused stage's file checksums. Changed inputs/code/configuration
require a new output directory. Preflight still rereads/reconstructs upstream
inputs on resume, but completed per-entity/enrichment stages are reused. Reports
can be rebuilt separately from verified saved results without original inputs.

Each stage is generated in a temporary directory and renamed only after success.
INCOMPLETE.txt remains until every requested stage succeeds; COMPLETE.txt is then
written. Existing outputs are never replaced through a force flag. A run lock
prevents concurrent writers. After an abrupt process kill, inspect the process
before manually removing a stale .run.lock. Interrupted temporary folders can be
preserved for inspection; they are not treated as completed stages.

If upstream inputs moved, provide comparison.input_overrides, for example:

~~~json
"input_overrides": {"wt":"relocated/WT.bed", "mt":"relocated/KO.bed"}
~~~

Hashes must still match. Overrides can cover all names recorded in upstream
RUNINFO.json, such as include, exclude, gtf, annotation or mapping. Relative
paths originally recorded by ChromHMMTools require explicit overrides rather
than guessing its original working directory.

## Output layout

~~~text
results/post01/
  manifest.json                    configuration, versions, input/code hashes
  COMPLETE.txt                     only on complete success
  regions/
    regions.json, regions.bed, region_qc.tsv
  models/
    models.json                    emissions, labels, optional context tables
  comparison_B113_vs_WT/
    all_overlap_tiles.bed           chromosome/start/end/WT-state/KO-state
    global_state_coverage.tsv       per-condition retained coverage
    entity_metrics.tsv             coverage, expression, changed fraction
    entity_state_pairs.tsv         observed per-entity pair support
    group_state_pairs.tsv          complete state-pair grids per group
    group_states.tsv               WT/KO state proportions
    selected_entities.tsv          selected entities with expression and IDs
    selected_transitions.bed        genomic pair selection, before entity filters
    genomic_context.tsv            disjoint feature bp counts/fractions
    expression_id_audit.tsv         identifier matching outcomes
    analysis.json                  saved numerical report inputs
  enrichment_B113_vs_WT/
    go_enrichment.tsv
    background_audit.tsv
    gene_universes.tsv              exact background/selected genes
    enrichment.json
  concordance/
    concordance.tsv, concordance.json
  report/
    report.html                    self-contained with embedded figures
    *.png                          separately reusable figures
~~~

all_overlap_tiles.bed is a five-column interval table with two state labels;
its fifth column is **not a conventional BED score**. Use the documented schema.
Every stage also has stage.json checksums. Reports include paths/configuration;
review those before sharing outside the project.

## Relationship to the historical R analysis

| Historical analysis | New reproducible implementation |
|---|---|
| Manual emission reordering/naming | Explicit model state-ID map and label/color/order table |
| WT/KO transition heatmaps | Verified upstream comparison plus complete coordinate-level overlaps |
| ChIPseeker annotation and distance filters | Explicit GTF windows and separate disjoint genomic-context annotation |
| Up/down transcript subsets from workspace objects | Configured saved DE table, thresholds and ID audit |
| TSS/upstream/downstream and broad promoter selections | Named strand-aware windows with fixed coordinates |
| Manually padded state matrices | State-ID keyed complete model axes; unsupported probabilities remain NA |
| Frequency/percentage plots | Base-pair weighted and equal-entity composition with recorded denominators |
| baseMean high/moderate/low grouping | Configurable abundance bins and declared measurement column |
| E10→E8/E9 annotation and GO | Explicit pair selection, selected entities, context charts and saved enrichment universe |
| Manual interactive plotting | Saved PNGs and one portable HTML report |
| Dependence on an old R workspace | Explicit inputs, provenance, staged completion and verified resume |

The intended biological questions are preserved. Old numerical results may change
because identifiers, coordinate definitions, selection units and enrichment
backgrounds are now explicit. These changes must be reviewed on real data before
comparing conclusions with an earlier figure.
