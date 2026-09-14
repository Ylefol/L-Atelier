# Remove Low Variance Features

Removes the bottom percentage of features by variance. A conservative
noise-removal approach that doesn't bias toward high-variance features -
it simply removes features that are clearly uninformative.

Preferred over
[`POSEIDON_filter_top_variable()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_filter_top_variable.md)
for methods like sCCA where you want to preserve cross-dataset
correlations, not just variance.

## Usage

``` r
POSEIDON_filter_low_variance(X, bottom_pct = 0.5, verbose = TRUE)
```

## Arguments

- X:

  Data matrix (samples x features) or a named list of data matrices.

- bottom_pct:

  Numeric. Percentage of lowest-variance features to remove (default =
  0.5, removes bottom 50%). Value between 0 and 1.

- verbose:

  Logical. Print filtering summary (default = TRUE).

## Value

Filtered data in same format as input (matrix or list of matrices).
Column names are preserved to track which features were retained.

## Details

Unlike
[`POSEIDON_filter_top_variable()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_filter_top_variable.md)
which selects features based on high variance (potentially biasing
results), this function only removes features that are clearly
uninformative (very low variance).

This is more appropriate for integration methods like sCCA where the
goal is to find cross-dataset correlations, not necessarily
high-variance features. A feature with moderate variance but strong
cross-dataset correlation is valuable and should be retained.

Suggested values:

- 0.25 - Conservative, removes only bottom 25%

- 0.50 - Moderate, removes bottom half (default)

- 0.75 - Aggressive, keeps only top 25%

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove bottom 50% of features by variance
X_filtered <- POSEIDON_filter_low_variance(X, bottom_pct = 0.5)

# For sCCA with multiple datasets
scca_data <- list(ATAC = atac_matrix, ChIP = chip_matrix, RNA = rna_matrix)
scca_filtered <- POSEIDON_filter_low_variance(scca_data, bottom_pct = 0.5)

} # }
```
