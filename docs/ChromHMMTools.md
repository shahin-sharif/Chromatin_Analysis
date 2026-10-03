# ChromHMMTools.py

Compare two ChromHMM segmentation BED files, quantify state overlaps and spatial adjacency, and create a portable HTML report. This consolidates the roles of `CompareChromHMM_new.py` and `ReportChromHMM.py` into one tool with `compare`, `report`, and `run` subcommands.

## Installation and first run

Python 3.8+ and NumPy are required. Matplotlib is required for report images; `--no-plots` produces tables and HTML without it.

```bash
python -m pip install numpy matplotlib
python ChromHMMTools.py run \
  --wt WT_segments.bed --mt MT_segments.bed \
  --state-mode shared --bin-size 200 \
  --write-change-bed --outdir chromhmm_comparison
```

Use `shared` **only when the two segmentations use the same state definitions**, for example from the same fitted model. The tool cannot infer this from state names. WT/MT are the labels for the two input conditions, not requirements on experimental design.

### Choosing state correspondence

- `--state-mode shared`: matching labels represent corresponding states. The union of labels forms both matrix axes; unsupported rows are NA.
- `--state-mode distinct`: independently trained or unknown state definitions. Separate WT and MT axes may have different sizes. Identical state numbers are not treated as conservation. All positive overlap pairs can appear in the ranking.
- `--state-mode mapped --state-map mapping.tsv`: explicitly map MT labels to WT labels, using a two-column, headerless `MT_state<TAB>WT_state` file. Every observed MT label must occur exactly once; targets must be distinct observed WT labels. Many-to-one state collapsing is deliberately unsupported.

For independently trained models, inspect their emission profiles and biological annotations before asserting correspondence. Genomic overlap alone does not establish equivalent state meaning.

## Separate computation and reporting

```bash
python ChromHMMTools.py compare \
  --wt WT_segments.bed.gz --mt MT_segments.bed.gz \
  --state-mode distinct --outdir comparison

python ChromHMMTools.py report --indir comparison --outdir report \
  --min-bp 10000 --min-fraction 0.05 --top-k 20
```

`report` consumes this tool's count matrices and `RUNINFO.json`; it is not a reader for arbitrary legacy output folders. It recomputes normalized values from counts. To regenerate a report in the same directory, use `report --indir comparison --outdir comparison --force`. `run` does both stages in a single command.

`--min-bp` is in **base pairs**, not megabases. `--min-fraction` is the minimum P(MT state | WT state). Rankings are descriptive: decreasing overlap bp, then conditional probability. They are not statistical significance tests.

## Input and coordinate rules

- Segmentation inputs are tab-delimited BED with at least four fields: chromosome, start, end, state. Extra fields are ignored. BED coordinates are zero-based, half-open.
- Plain text and `.gz` inputs are accepted. Blank lines, `#` comments, and standard `track`/`browser` lines are skipped.
- Inputs are sorted internally. Overlapping segmentation intervals, empty labels, negative coordinates, and nonpositive lengths are errors. Adjacent intervals with the same state are merged.
- Chromosome names are preserved and must match between files. Default behavior retains all chromosomes. `--autosomes-only` selects human 1–22 / chr1–chr22; it is not a general species-aware filter.
- Both files must describe the same reference assembly. Assembly identity cannot be established automatically from BED coordinates.
- Segments and annotations are held in memory; this avoids expanding a genome into per-base arrays but is not a streaming implementation.

## Masks, annotations, and change intervals

Use `--include-bed regions.bed` and/or `--exclude-bed blacklist.bed`. Only the first three BED fields matter. Overlapping mask intervals are unioned. Include is applied first, exclusion second. An empty include file selects nothing and produces an error.

Optional annotation is either:

```bash
--gtf gencode.annotation.gtf.gz --promoter-up 1000 --promoter-down 200
```

or `--anno-bed classes.bed` (fourth field is the class label).

GTF annotation uses **gene** rows. Gene bodies are `genic`; promoter windows are strand-aware and include the TSS base. Thus up=1000/down=200 gives 1201 bases before clipping at coordinate zero. No upper chromosome boundary is inferred; overlap is naturally restricted to segmentation coverage. Use BED annotation for other feature definitions.

Overlap tiles are split at annotation boundaries. Overlapping classes are resolved using `--anno-priority promoter,genic,intergenic`; listed labels take precedence, and unlisted labels follow alphabetically. Unannotated territory is `intergenic`. Each shared base is assigned once, so annotation matrices sum to the global overlap matrix.

`--write-change-bed` writes merged adjacent intervals with the same WT→MT pair. `--min-change-len N` filters after merging. In shared/mapped mode only different states are included. In distinct mode **all state mappings** are written; the filename `state_changes.bed` does not imply biological change in that mode.

## What the numbers mean

**Overlap matrices:** count exact shared covered bases after masks. WT states are rows; MT states are columns. Row normalization gives P(MT | WT), column normalization P(WT | MT), and global normalization the fraction of all shared bases. Gaps/unshared bases are excluded from these denominators and recorded separately in metadata.

**Spatial transitions:** for each condition independently, use complete bins of `--bin-size` (default 200), anchored to coordinate zero. Count successive bins only when both are wholly contained in a retained state interval and are contiguous. Partial bins at boundaries are excluded; no transition crosses a gap, excluded region, or chromosome boundary. Direction is increasing genomic coordinate. These counts describe the observed segmentation, **not ChromHMM's learned HMM transition parameters**. They use each condition's retained coverage, which may differ between conditions; restrict with an appropriate common include mask if equal territory is required. Choose the bin size used to generate your segmentation.

**Zero support:** normalized rows/columns with no observations are NA, not zero or an inferred probability. Self-transition probabilities with no outgoing bin pairs are also NA.

**Report summaries:** retention is the matching-label fraction where correspondence is declared. Shannon entropy and `2**entropy` describe the spread across MT states, conditional on each WT state. Log2 enrichment is `log2((observed_bp + 1)/(expected_bp + 1))`, where expectation comes from overlap-matrix marginals. It is a descriptive, pseudocount-dependent measure.

## Outputs

| File | Contents |
| --- | --- |
| `overlap.bp.tsv` | Integer WT × MT overlap bases |
| `overlap.row.tsv`, `.col.tsv`, `.global.tsv` | Normalized overlap matrices |
| `WT.trans.counts.tsv`, `MT.trans.counts.tsv` | Empirical adjacent-bin counts |
| `WT.trans.probs.tsv`, `MT.trans.probs.tsv` | Row-normalized adjacency probabilities |
| `annotation_N.overlap.bp.tsv` | Exact overlap counts by annotation; filenames mapped to classes in metadata |
| `state_changes.bed` | Optional BED6 of changed/mapped intervals |
| `RUNINFO.json` | Parameters, input paths/SHA-256 hashes, shared/unshared coverage and bin counts |
| `state_summary.tsv` | Report: per-WT-state support, retention, entropy and effective targets |
| `transition_summary.tsv` | Report: outgoing bin-pair counts and self-transition probabilities |
| `top_pairs.tsv` | Report: ranked positive overlap pairs and conditional probabilities |
| `report.html` | Self-contained HTML with embedded PNGs, tables and provenance |

HTML embeds its images and can be moved or shared as one file. It includes input paths and parameters from the run. No JavaScript or network resources are needed.

## File safety and migration

Inputs are unchanged. Outputs are prepared in a temporary directory and validated before installation. Existing filenames require `--force`; output/input aliases are rejected. Individual files are installed by atomic replacement, but the entire collection is not a filesystem transaction. Use a fresh directory for each comparison: `--force` replaces generated matching filenames and does not clean up unrelated or obsolete files from earlier runs.

This is a revised interface, not a drop-in replacement for legacy command lines or matrix filenames. Use the documented subcommands and regenerate comparisons before using `report`. The original scripts on your drive are preserved. This version focuses on exact overlap matrices, empirical spatial transitions, annotation-stratified overlap, changed intervals and reports; it does not fit HMMs, align emission models automatically, or conduct differential significance testing.

## Validation

Twenty regression tests cover independently known overlap counts, a randomized complete-bin oracle, gaps and partial bins, mask subtraction, empty includes, annotation base conservation, GTF strand conversion, zero support, rectangular matrices, explicit mapping, label-aligned reports, escaping, and overwrite protection. A plotting smoke check also exercises embedded report images. No representative real ChromHMM segmentation pair has been tested yet.
