# Filter Low-Count Features (Peaks/Genes)

Filters features (peaks, genes, etc.) with insufficient counts across
samples. Commonly used to remove unreliable features before differential
analysis or integration.

## Usage

``` r
POSEIDON_filter_low_counts(
  quant_result,
  min_count = 10,
  min_samples = 2,
  verbose = TRUE
)
```

## Arguments

- quant_result:

  A list containing at minimum:

  - counts: Matrix of counts (features x samples)

  - annotation: Data.frame with feature annotations (optional)

- min_count:

  Integer. Minimum count threshold (default = 10).

- min_samples:

  Integer. Minimum number of samples that must meet min_count threshold
  (default = 2).

- verbose:

  Logical. Print filtering summary (default = TRUE).

## Value

The input list with filtered counts and annotation.

## Details

A feature is retained if at least `min_samples` samples have counts \>=
`min_count`. This filtering removes features that are:

- Not detected in most samples

- Too low-count for reliable statistical analysis

- Likely to introduce noise in downstream analyses

Common thresholds:

- ATAC-seq peaks: min_count = 10, min_samples = 3

- RNA-seq genes: min_count = 10, min_samples = 2-3

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter peaks with at least 10 counts in at least 3 samples
counts_filtered <- POSEIDON_filter_low_counts(counts, min_count = 10, min_samples = 3)

} # }
```
