# Run REDItools2 RNA editing detection

Detects RNA editing events in a BAM file by running REDItools2, a
pysam-based pileup engine that identifies positions with nucleotide
substitutions relative to the reference genome. The canonical use case
is A-to-I (A\>G) editing detection in RNA-seq data.

## Usage

``` r
CYAN_run_reditools(
  bam,
  reference,
  output_dir,
  work_dir = tempdir(),
  sample_name = NULL,
  strand_mode = 0L,
  min_coverage = 1L,
  min_base_quality = 30L,
  min_read_quality = 20L,
  min_read_length = 30L,
  min_edits = 0L,
  strict = FALSE,
  region = NULL,
  omopolymeric_file = NULL,
  splicing_file = NULL,
  bed_file = NULL,
  require_dup_marked = TRUE,
  keep_tmp = FALSE,
  save_rds = TRUE,
  verbose = TRUE,
  debug = FALSE
)
```

## Arguments

- bam:

  Character. Path to a sorted, indexed BAM file.

- reference:

  Character. Path to the reference FASTA (must be indexed with
  `samtools faidx`).

- output_dir:

  Character. Directory for final result files (created if absent).

- work_dir:

  Character. Directory for intermediate/temporary files (created if
  absent). Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html). Can point to an
  external drive for large datasets.

- sample_name:

  Character or `NULL`. Label used in output filenames. Derived from the
  BAM filename if `NULL`.

- strand_mode:

  Integer. Library strandedness: 0 = unstranded (default), 1 =
  second-strand, 2 = first-strand.

- min_coverage:

  Integer. Minimum read coverage at a position (default 1, matching
  REDItools2 native default). Increase (e.g. 10) to suppress
  low-coverage noise post-hoc or here to reduce output size.

- min_base_quality:

  Integer. Minimum per-base quality score (default 30, matching
  REDItools2 native default).

- min_read_quality:

  Integer. Minimum read mapping quality (default 20, matching REDItools2
  native default).

- min_read_length:

  Integer. Minimum read length; shorter reads are discarded (default 30,
  matching REDItools2 native default).

- min_edits:

  Integer. Minimum number of edited reads at a position (default 0,
  matching REDItools2 native default). Filtering on this value post-hoc
  is straightforward:
  `result[result$frequency * result$coverage >= 3, ]`.

- strict:

  Logical. If `TRUE`, only positions with detected substitutions are
  written to the output. Default `FALSE` (matching REDItools2 native
  default) — all covered positions are reported. Set to `TRUE` to
  substantially reduce output size when only editing sites are of
  interest.

- region:

  Character or `NULL`. Restrict analysis to a genomic region in samtools
  format (e.g. `"chr21:1-10000000"`).

- omopolymeric_file:

  Character or `NULL`. Path to homopolymeric regions file. Recommended
  to reduce false positives.

- splicing_file:

  Character or `NULL`. Path to splice-site positions file. Recommended
  to reduce mismatch artefacts near junctions.

- bed_file:

  Character or `NULL`. Path to a BED file of target regions. Analysis is
  restricted to these regions when supplied.

- require_dup_marked:

  Logical. If `TRUE` (default), hard-stops before running REDItools2
  unless the BAM has at least one read flagged as a duplicate (SAM FLAG
  0x400). REDItools2's own duplicate filter relies entirely on this flag
  already being set upstream (e.g. via `samtools markdup` or Picard
  `MarkDuplicates`) – on a BAM that was never duplicate-marked, that
  filter is a silent no-op. Set to `FALSE` only if you are confident
  this BAM genuinely has no duplicates to flag.

- keep_tmp:

  Logical. If `FALSE` (default), intermediate files in `work_dir` are
  removed after completion.

- save_rds:

  Logical. If `TRUE` (default), the parsed result is saved as an RDS
  file in `output_dir`.

- verbose:

  Logical. Print `[CYAN]` R-level progress messages (default `TRUE`).
  Does not affect REDItools2's own output.

- debug:

  Logical. If `TRUE`, passes `-V` to REDItools2 and routes its
  stdout/stderr to the console (default `FALSE`). Use this when a run
  fails and you need to see the raw Python output to diagnose the
  problem.

## Value

An object of class `cyan_editing_result`: a `data.frame` with one row
per reported genomic position and columns:

- region:

  Chromosome / contig name

- position:

  1-based genomic position

- reference:

  Reference base

- strand:

  Strand (0/1/2)

- coverage:

  Read coverage after quality filtering

- mean_quality:

  Mean base quality at this position

- count_A, count_C, count_G, count_T:

  Per-base read counts

- all_subs:

  Space-separated substitution codes (e.g. `"AG"`, `"AG TC"`); `"-"` if
  no substitution detected

- frequency:

  Fraction of reads supporting the most frequent substitution (0–1)

The object additionally carries a `params` attribute listing the call
parameters for reproducibility.

## Details

REDItools2 is run as a subprocess inside a managed basilisk Python
environment. The raw output (tab-delimited text) is written to
`output_dir` as a checkpoint before parsing, so the expensive
computation is never lost. The parsed `data.frame` is saved alongside as
an RDS file (opt-out with `save_rds = FALSE`).

**File outputs** (in `output_dir`):

- `{sample}_{date}_reditools.txt` — full raw REDItools2 table (all
  covered positions)

- `{sample}_{date}_reditools_subs.txt` — substitution sites only (rows
  where `all_subs != "-"`)

- `{sample}_{date}_reditools.rds` — parsed `cyan_editing_result` (if
  `save_rds = TRUE`)

## Parameter notes

- `strict`:

  When `TRUE` (default), only positions with at least one detected
  substitution are emitted. Set to `FALSE` to receive all covered
  positions (large output).

- `min_coverage`:

  Minimum number of reads covering a position (after base-quality
  filtering) for it to be reported.

- `min_edits`:

  Minimum number of edited reads at a position.

- `strand_mode`:

  0 = unstranded; 1 = second-strand (dUTP/fr-firststrand); 2 =
  first-strand (fr-secondstrand). Match to your library protocol.

- `omopolymeric_file`:

  Path to a file of homopolymeric regions to exclude — strongly
  recommended to reduce false positives. Can be generated by REDItools2
  with the `-c` flag (not yet exposed here).

- `splicing_file`:

  Path to a file of splice-site positions to exclude — recommended to
  avoid mismatch artefacts near splice junctions.

## See also

[`CYAN_check_bam`](https://ylefol.github.io/L-Atelier/CYAN/reference/CYAN_check_bam.md),
[`CYAN_check_references`](https://ylefol.github.io/L-Atelier/CYAN/reference/CYAN_check_references.md)
