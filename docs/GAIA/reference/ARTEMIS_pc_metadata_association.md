# PC-Metadata Association Test

Runs PCA on a log2-scale matrix, then tests the association of each
principal component with each metadata variable. Returns a p-value
matrix that summarises which variables explain sample-level variance.

Intended use: run after loading and filtering, before differential
analysis, to distinguish expected biological drivers (group, timepoint)
from unexpected technical ones (plate, batch). If a technical variable
is significant, consider correcting with `POSEIDON` before proceeding.

## Usage

``` r
ARTEMIS_pc_metadata_association(
  matrix,
  sample_meta,
  n_pcs = 10,
  ntop = NULL,
  sample_col = NULL,
  metadata_cols = NULL,
  scale = TRUE,
  categorical_threshold = 10L,
  verbose = TRUE
)
```

## Arguments

- matrix:

  Numeric matrix, features x samples (log2-scale). Must have no NA
  values — filter with
  [`HADES_filter_olink_proteins()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_olink_proteins.md)
  first if needed. Rownames = feature IDs, colnames = sample IDs.

- sample_meta:

  Data.frame with one row per sample. Matched to `colnames(matrix)` by
  `sample_col` (or rownames if NULL).

- n_pcs:

  Integer. Number of top PCs to test. Default: 10.

- ntop:

  Integer or NULL. Number of most variable features (by row variance) to
  use for PCA. NULL = use all features. Default: NULL.

- sample_col:

  Character or NULL. Column in `sample_meta` holding sample IDs matching
  `colnames(matrix)`. NULL = use rownames. Default: NULL.

- metadata_cols:

  Character vector or NULL. Columns of `sample_meta` to test. NULL = all
  columns except `sample_col`. Default: NULL.

- scale:

  Logical. Scale features to unit variance before PCA. Recommended for
  NPX data where proteins have different dynamic ranges. Default: TRUE.

- categorical_threshold:

  Integer. Variables with at most this many unique non-NA values are
  treated as categorical (Kruskal-Wallis). Variables with more unique
  values are treated as continuous (Spearman). Default: 10.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_pc_assoc"` containing:

- pvalues:

  Numeric matrix (n_pcs x n_vars). Raw p-values.

- effect_sizes:

  Numeric matrix (n_pcs x n_vars). Spearman rho for continuous
  variables; eta-squared for categorical variables.

- test_used:

  Named character vector. Test applied per variable: "spearman",
  "kruskal", or "skipped" (degenerate variable).

- var_explained:

  Numeric vector. Proportion of variance explained per PC (e.g. 0.25 =
  25%).

- pca:

  The `prcomp` object for further use or custom plotting.

- params:

  List of parameters used.

## Details

Test selection per variable:

- Continuous (\> `categorical_threshold` unique values): Spearman rank
  correlation between PC scores and variable values.

- Categorical (\<= `categorical_threshold` unique values, or
  character/factor): Kruskal-Wallis rank sum test across groups.

- Skipped: variables with only one unique value, or where fewer than 2
  complete observations per group exist.

Effect sizes: Spearman rho (range -1 to 1); eta-squared for
Kruskal-Wallis (H - k + 1) / (n - k), where H is the test statistic, k
is the number of groups, and n is the number of observations. Negative
eta-squared is set to 0.

## Examples

``` r
if (FALSE) { # \dontrun{
ol <- HADES_filter_olink_proteins(HADES_filter_olink(ol_raw))

assoc <- ARTEMIS_pc_metadata_association(
  matrix      = ol$wide,
  sample_meta = ol$sample_meta,
  sample_col  = "SampleID",
  n_pcs       = 10
)

print(assoc)
AETHER_plot_pc_association(assoc)
} # }
```
