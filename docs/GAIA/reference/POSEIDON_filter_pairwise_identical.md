# Filter Features with Pairwise Identical Values (Debug Filter)

Aggressively filters features where ANY pair of samples has identical
values. Such features will produce NaN when scaled in CV folds
containing that pair (due to zero variance → division by zero).

## Usage

``` r
POSEIDON_filter_pairwise_identical(X, tolerance = 1e-10, verbose = TRUE)
```

## Arguments

- X:

  A matrix (samples x features) or a list of matrices.

- tolerance:

  Numeric. Values within this tolerance are considered identical
  (default = 1e-10).

- verbose:

  Logical. Print filtering summary (default = TRUE).

## Value

Filtered matrix or list of matrices with problematic features removed.

## Details

This is an AGGRESSIVE filter intended for debugging small-sample
scenarios. It removes features where any two samples have identical (or
near-identical) values, which would cause scaling to fail in CV when
those samples end up in the same training fold.

For a feature to pass this filter, ALL pairwise comparisons between
samples must show some variance. This can remove a substantial number of
features, especially with:

- Count data with many zeros

- Averaged technical replicates that converge to similar values

- Small sample sizes

WARNING: This filter is NOT recommended for production analyses as it
may remove biologically meaningful features. Use only for debugging or
when small sample sizes make CV scaling problematic.

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter a single matrix
X_filtered <- POSEIDON_filter_pairwise_identical(X)

# Filter a list of matrices (e.g., for multi-omics)
data_filtered <- POSEIDON_filter_pairwise_identical(scca_data)

} # }
```
