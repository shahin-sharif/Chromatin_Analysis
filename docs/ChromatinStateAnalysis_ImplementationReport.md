# ChromatinStateAnalysis: implementation and review report

Prepared 4 October 2026. Version 0.1.0.

## Delivered tool

ChromatinStateAnalysis.py is a Python command-line entry point backed by the
chromatin_state/ modules. It is located alongside the existing ChromHMMTools.py
in Chromatin_Analysis. The complete folder should be retained together.

The workflow covers the post-ChromHMM questions recorded in AllSegmentation.R
and ChrHMM.R: interpreting model emissions, measuring WT/KO state relationships,
locating changes around genes/transcripts, connecting them to expression groups,
examining expression abundance, annotating selected transitions and testing
functional gene sets. It replaces dependence on an interactive R workspace with
explicit files, configuration, saved stages and an auditable report.

The approved build is included in the Chromatin_Analysis repository. No experimental
server analysis was launched during implementation.

## Responsibilities of the three analysis levels

| Level | Responsibility | Important distinction |
|---|---|---|
| Original ChromHMM | Fit chromatin-state models and produce segmentation, emissions and optional context enrichment | This new tool does not train or select an HMM |
| ChromHMMTools.py | Exact global WT/KO overlap, normalized overlaps, annotation-stratified counts, empirical spatial adjacency, changed intervals and report | It establishes the reproducible segmentation comparison, not expression-linked inference |
| ChromatinStateAnalysis.py | Verified downstream windows, expression integration, grouped compositions, selected transition annotations, enrichment and report | It consumes the first two levels and adds biological interpretation with explicit counting rules |

The existing ChromHMMTools code was not modified. Its tested interval engine is
reused to reconstruct complete coordinate-level WT/KO tiles from the original,
hash-verified inputs. The resulting count matrix must exactly match the saved
upstream overlap.bp.tsv. Unchanged tiles are retained because they are part of
the denominator for statements such as “half of this promoter changed state.”

Upstream adjacent-bin transitions describe spatial adjacency within a segmentation.
WT→KO overlap describes differences between conditions at matching coordinates.
Neither is silently substituted for ChromHMM's fitted transition parameters.

## Implemented stages

### 1. Deep preflight and provenance

The check command validates configuration keys, required files, original upstream
SHA-256 hashes, saved comparison counts, state labels, emission probabilities,
chromosome compatibility, annotation features, expression columns/IDs, selections
and optional gene-set resources. Unknown options and zero-match expression or
gene-set inputs fail rather than silently producing plausible empty results.

Models are explicitly associated with comparisons. Shared state mode requires
one declared model. Distinct and explicitly mapped upstream comparisons are
supported. Emission labels such as 1 and segmentation labels such as E1 require
an explicit one-to-one map; numerical position is never used as identity.

The run records input and source-code hashes, resolved configuration, Python and
package versions. A relocated input may be supplied through input_overrides only
when its content hash matches the upstream run.

### 2. Model interpretation

The report presents emission probabilities on their original 0–1 scale, with
reviewed state labels, colors and order. Unobserved model states remain represented.
Optional ChromHMM neighborhood and annotation-overlap tables are imported with
explicit units. Context tables with mixed genome percentages and fold enrichments
must be separated using named column selection.

Biological state names are provided annotations. The program does not infer
“super-enhancer,” “poised” or other biological classifications automatically from
a state number.

### 3. Explicit genomic windows

A consistent GTF defines gene or transcript identity, strand, boundaries and
exons. Default windows cover TSS core, upstream/downstream promoter, broad
promoter, body and TES. Windows can be replaced by named user definitions;
body mode can be genomic span or exon union. An exact-ID selection file can
restrict the entities being studied.

Coordinates are half-open, with strand-reflected TSS/TES boundary definitions.
For example, [-25,25) is a 50-base TSS window. Clipping, empty regions and missing
chromosomes remain visible in region_qc.tsv. Custom named BED regions can be
included, with optional explicitly labeled nearest-TSS association. Ties are
ambiguous rather than arbitrarily assigned.

### 4. Per-entity state measurements

Each window is intersected with the complete WT/KO tile set. Outputs retain
effective length, shared coverage, pair-support bp, pair fractions, changed
fraction where state correspondence permits it, and dominant-pair ties.

The default minimum shared coverage is 80%, configurable. Low-coverage windows
remain in the entity table but do not enter summaries or gene selection. Missing
territory is not counted as quiescent state or as unchanged sequence.

### 5. Expression integration

Saved gene- or transcript-level differential-expression results are joined by
explicit IDs. This does not rerun DESeq2. Fold-change, adjusted-p-value and
abundance columns are configurable, with an optional explicit one-to-one ID map.

Expression classes are up, down, not_significant, not_tested, unmatched and
not_supplied. Significance and effect-size thresholds are configurable; absent
statistics never become an “unchanged” claim. Identifier audits preserve both
unmatched expression rows and annotation entities lacking expression results.

The original baseMean grouping is supported through configurable abundance
breaks. Default bins are [0,100), [100,500) and [500,infinity). The configured
measurement is recorded; baseMean is not renamed TPM.

### 6. Grouped state composition and transitions

For each region type and expression/abundance group, the workflow produces:

- Complete WT × KO state-pair matrices.
- Attributed base-pair fractions.
- Equal-entity mean fractions.
- P(KO state | WT state), with NA for unsupported rows.
- WT and KO state-composition tables and figures.

Long regions dominate base-pair weighting but not equal-entity means. Overlapping
genes/transcripts can attribute the same genomic bases to multiple entities;
those totals are clearly distinguished from unique-genome coverage.

Multiple comparisons have a concordance table with separate values and a flag
for consistent up/down expression direction. A separate dominant-pair agreement
flag is available only for eligible, untied results on corresponding state axes.
This is descriptive agreement,
not a new statistical test or evidence that shared-reference comparisons are
independent.

### 7. Selected transitions and genomic annotation

A selection can specify one or more pairs, such as E10→E8 and E10→E9, plus minimum
bp/fraction support, region types and expression groups. Selected entity tables
retain gene IDs, coordinates through region IDs, expression statistics, association
method and the actual supporting fraction.

Genomic pair tiles are also exported before entity/expression filters. Their
context is summarized by exact bp under a disjoint promoter > exon > intron >
genic_unresolved > intergenic policy, with pie charts. Genes without exon data
are not automatically labeled intronic. This fulfills the historical context
annotation question with an explicit partition rather than reproducing accidental
R data-frame assignments.

Nearest-TSS associations for custom distal regions remain proximity assignments.
Expression joins through them require a separate explicit option and are labeled.
They are not presented as demonstrated enhancer–gene interactions.

### 8. Functional enrichment

The Python implementation accepts local, versioned term-to-gene memberships.
For GO, the input must already contain the intended ancestor-expanded memberships
and filtering choices. The current tool does not download an ontology, propagate
GO relationships or infer annotation IDs.

For each comparison × selection × window, the program saves the eligible gene
universe, annotation losses and selected genes, then calculates a one-sided
hypergeometric over-representation test. Benjamini–Hochberg correction includes
all size-eligible terms, including zero-hit terms, within that declared family.
Term overlaps, effect sizes, p-values, adjusted p-values and contributing genes
are exported. Ranked dot plots are not restricted to significant terms; the report
explicitly identifies them as ranked tested terms.

When expression is supplied, GO backgrounds exclude untested and unmatched
entities. Transcript selections are collapsed to unique gene IDs. The resource
and background definitions remain in the output rather than relying on an
implicit package-wide universe.

### 9. Figures and report

One self-contained HTML report embeds figures and previews the output tables.
Separate PNG files are available for reuse. It contains model emissions, optional
context panels, global state overlap, expression/abundance-stratified transition
heatmaps, state-composition bars, genomic-context pies, enrichment dot plots,
coverage/ID audits and configuration/provenance.

No zero-to-log substitution is used. Probabilities, fractions, bp counts and
fold enrichments retain their stated units. Full tables are saved even when
only a preview is included in HTML.

### 10. Recovery

Each numerical/report stage is published only after completion, with checksums.
A failed plotting stage does not invalidate completed window/enrichment stages.
Resume verifies input/configuration/source fingerprints and saved stage files.
Report-only regeneration works from saved data without the original inputs.
An existing unrelated directory is not overwritten; concurrent runs are locked.

## Validation evidence

The local suite passed **54 tests: 20 existing ChromHMMTools tests and 34 new
post-segmentation tests**. Python 3.12 was used with the available NumPy and
matplotlib runtime.

Tests cover:

- Whole-pipeline known-answer overlaps, windows, DE classes, selected genes and GO.
- Both strands, exon unions, clipped windows and transcript/gene identities.
- Base-pair versus equal-entity weighting and low/missing coverage.
- Shared, distinct and explicit mapped models; unsupported model states.
- Upstream mask reuse, relocated hashed inputs, altered input/count rejection.
- Duplicate and unmatched expression IDs, explicit maps and ambiguity rejection.
- Optional nearest-TSS assignment and explicit proximity expression joining.
- Genomic-context base conservation.
- Randomized interval intersections checked against an independent per-base oracle.
- Hypergeometric tails checked against exact integer-combination calculations over
  many small populations, plus a known BH correction example.
- Multiple comparisons, overwrite protection, modified-output detection, verified
  resume, simulated report failure and report recovery after source input deletion.

The plotted synthetic end-to-end runner completed and verified expected outputs.
The example's two selected genes yield a known enrichment p-value of 1/3 and BH
adjusted p-value of 2/3. This intentionally nonsignificant example tests arithmetic,
not a biological discovery. Representative PNGs were visually inspected.

GitHub Actions runs the regression suite and both complete examples on pushes and
pull requests. Local validation evidence is recorded above; per-commit remote
results are available in the repository’s Actions tab.

## What remains to validate with experimental data

No full BAZ1B experimental analysis or genome-scale performance benchmark has
been run with the new tool. Before interpreting biological results, use the
actual current model, a ChromHMMTools comparison generated from its segmentations,
a matching GTF and the intended full DE tables. Review state correspondence,
identifier-match rates, coverage exclusions, window definitions and gene-set
mapping before comparing with historical figures.

The method is intentionally descriptive for chromatin differences. It does not
implement ChromDiff, fit differential chromatin models, adjust gene-set tests
for gene-length/transcript-count selection bias, prove regulatory causation or
recover biological replicates from pooled segmentations.

## Next discussion points

- Gene versus transcript as the primary analysis unit, with the other optionally
  run as a separate configuration.
- Which current 13-state model and reviewed state labels should be authoritative.
- Whether separate B113/B142 chromatin segmentations are available or only pooled KO.
- Final promoter/TES windows, minimum shared coverage and DE effect-size thresholds.
- The GO membership source/release and exact mapping to the selected annotation.

These are scientific configuration choices. The software foundation, reproducible
example, tests and full command documentation are now available for that review.
