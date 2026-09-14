# Prepare BigWig Signal Data for 2D Density Scatter Comparison

Extracts mean BigWig signal over genomic regions for two groups, handles
replicate averaging, optional log2 ratio computation, and low-signal
filtering. Designed to feed directly into
[`AETHER_plot_density_scatter`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_density_scatter.md).

## Usage

``` r
ARTEMIS_prepare_bw_comparison(
  sample_sheet,
  group_x,
  group_y,
  regions,
  bin_size = 1000L,
  genome = NULL,
  ratio = FALSE,
  pseudocount = 1,
  min_signal = NULL,
  sample_col = "sample_id",
  bw_col = "bw_path",
  coverage = NULL
)
```

## Arguments

- sample_sheet:

  data.frame. Must include columns identified by `sample_col` and
  `bw_col`.

- group_x:

  Named list defining the x-axis. See Details.

- group_y:

  Named list defining the y-axis. See Details.

- regions:

  One of: `"bins"` for genome-wide fixed-width bins, a path to a BED
  file (3-column, 0-based half-open), or a `GRanges` object.

- bin_size:

  Integer. Bin width in bp. Only used when `regions = "bins"`. Default
  1000.

- genome:

  Character or `NULL`. Genome identifier (e.g. `"hg38"`). Used for
  logging when `regions = "bins"`. Chromosome sizes are derived from the
  BigWig seqinfo; this parameter does not load a BSgenome package.
  Default `NULL`.

- ratio:

  Logical. If `TRUE`, compute log2 ratio for any axis whose list has
  exactly 2 names. Default `FALSE`.

- pseudocount:

  Numeric. Added to group means before log2 ratio:
  `log2((numerator + pseudocount) / (denominator + pseudocount))`.
  Prevents `log(0)` and dampens instability at near-zero signal.
  Symmetric — does not introduce directional bias. Default 1.

- min_signal:

  Numeric or `NULL`. Regions where the representative signal on *both*
  axes falls below this threshold are removed before ratio computation.
  Representative = single group mean (1-name list) or
  `max(group1, group2)` (2-name list). Regions enriched in only one
  group are retained. Default `NULL` (no filtering).

- sample_col:

  Character. Column in `sample_sheet` holding sample IDs. Default
  `"sample_id"`.

- bw_col:

  Character. Column in `sample_sheet` holding BigWig file paths. Default
  `"bw_path"`.

- coverage:

  `NULL` or an `artemis_bw_coverage` object returned by
  [`ARTEMIS_load_bw_coverage`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_load_bw_coverage.md).
  When provided, BigWig files are not read from disk — signal is
  extracted from the pre-loaded RLE coverage objects. Speeds up
  multi-region workflows by loading each file only once. Default `NULL`.

## Value

data.frame with columns `x`, `y`, `chr`, `start`, `end` (one row per
region). Attributes `x_label` and `y_label` store descriptive axis
labels that `AETHER_plot_density_scatter` uses as defaults.

## Details

**Input list logic:**  
Each of `group_x` and `group_y` is a named list mapping group labels to
character vectors of sample IDs:

- **1 name** — direct mean signal across all samples. E.g.
  `list(CTRL = c("s1", "s2"))`.

- **2 names + `ratio = TRUE`** — log2 ratio. The *first name is the
  denominator*, the *second is the numerator*:
  `log2((2nd + pc) / (1st + pc))`. E.g.
  `list(CTRL = c("s1","s2"), TREAT = c("s3","s4"))` yields
  `log2((TREAT + pc) / (CTRL + pc))`.

Mixed axes are valid (one direct, one ratio).

**Order of operations:** validate inputs → resolve paths → define
regions → extract mean signal per BigWig → average replicates per group
→ filter (min_signal) → compute log2 ratio → return data frame.
