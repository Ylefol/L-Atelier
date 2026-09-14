# Filter Top Variable Features

Selects the top N most variable features from each dataset. Essential
for reducing dimensionality before methods like sCCA that don't scale
well with very high feature counts (\>10k features).

Works with any numeric data: counts, normalized values, log-transformed,
etc.

## Usage

``` r
POSEIDON_filter_top_variable(X, n_top = 5000, verbose = TRUE)
```

## Arguments

- X:

  Data matrix (samples x features) or a named list of data matrices. For
  sCCA input, this should be the transposed format (samples in rows).

- n_top:

  Integer. Number of top variable features to keep per dataset. If a
  single value, applied to all datasets. Can also be a named list
  matching dataset names for dataset-specific values (default = 5000).

- verbose:

  Logical. Print filtering summary (default = TRUE).

## Value

Filtered data in same format as input (matrix or list of matrices).
Column names are preserved to track which features were retained.

## Details

Variance is calculated per feature (column) using
[`var()`](https://rdrr.io/r/stats/cor.html). Features are ranked by
variance and the top N are retained.

For multi-omics sCCA, typical values:

- ATAC-seq: 5,000 - 10,000 peaks

- ChIP-seq: 1,000 - 5,000 regions (often fewer to start)

- RNA-seq: 5,000 - 10,000 genes

If n_top exceeds the number of features in a dataset, all features are
kept with a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single matrix
X_reduced <- POSEIDON_filter_top_variable(X, n_top = 5000)

# List of matrices (for sCCA)
scca_data <- list(ATAC = atac_matrix, ChIP = chip_matrix, RNA = rna_matrix)
scca_reduced <- POSEIDON_filter_top_variable(scca_data, n_top = 5000)

# Different n_top per dataset
scca_reduced <- POSEIDON_filter_top_variable(
  scca_data,
  n_top = list(ATAC = 10000, ChIP = 2000, RNA = 5000)
)

} # }
```
