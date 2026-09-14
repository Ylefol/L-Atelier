# Compare multiple decoupleR methods

Runs multiple statistical methods on the same data and compares results.
Useful for assessing robustness and identifying consensus activities.

## Usage

``` r
ARTEMIS_decoupler_compare_methods(
  mat,
  network,
  methods,
  minsize = 5,
  norm_method = "log2cpm",
  consensus = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- mat:

  Numeric matrix of gene expression data (genes x samples), or an
  `artemis_norm` object from
  [`ARTEMIS_normalize_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_counts.md).
  Rownames must be gene identifiers matching the network.

- network:

  Prior knowledge network.

- methods:

  Character vector. Methods to compare. REQUIRED - no default.

- minsize:

  Integer. Minimum targets per source. Default: 5.

- norm_method:

  Character. Normalization to apply before running methods. `"log2cpm"`
  (default): log2(CPM + 1). `"none"`: pass through as-is. Normalization
  is applied once to the full matrix before any method runs to ensure
  consistency. If `mat` is an `artemis_norm` object, DESeq2 normalized
  counts are extracted and log2-transformed automatically.

- consensus:

  Logical. Compute consensus scores across methods. Default: TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

- ...:

  Additional arguments passed to each method.

## Value

A list with class "decoupler_comparison" containing:

- results:

  Named list of decoupler_result objects (one per method)

- summary:

  Data.frame summarizing results across methods, built from each
  method's raw activity scores. Comparable "mean activity" across
  sources, but NOT comparable across methods when methods differ in
  output scale (e.g. wsum's unbounded weighted-sum scores vs. ulm/mlm's
  t-like statistics vs. viper's NES-like scores) – a method with a much
  larger raw scale will dominate `mean_activity`/`min_activity`/
  `max_activity` here regardless of actual cross-method agreement. Use
  `summary_zscored` instead when comparing/ranking across methods with
  different native scales.

- summary_zscored:

  Same columns as `summary`, but each method's activity matrix is
  z-scored (matrix-wide mean/sd, same normalization `consensus` uses)
  before the per-source mean/sd/min/max/ agreement_score are computed –
  puts all methods on a comparable scale first, so `mean_activity`
  ranking and the min-max range reflect actual cross-method agreement
  rather than one method's raw output magnitude. Used by
  `AETHER_plot_method_agreement(scale = "zscore")`.

- correlations:

  Method-method correlation matrix

- method_stats:

  Per-method statistics

- consensus:

  Consensus activity scores (if consensus = TRUE)

- methods:

  Methods used

- n_methods:

  Number of methods

## Details

The summary table includes for each source (TF/pathway):

- Mean activity across methods

- Standard deviation across methods

- Number of methods where significant (p \< 0.05)

- Agreement score (1 - normalized SD)

This helps identify which activities are robust across different
statistical approaches vs. method-dependent.

## Examples

``` r
if (FALSE) { # \dontrun{
network <- APOLLO_get_collectri()
comparison <- ARTEMIS_decoupler_compare_methods(
  expr_matrix, network,
  methods = c("ulm", "mlm", "wsum")
)

# View summary
head(comparison$summary)

# View method correlations
comparison$correlations

} # }
```
