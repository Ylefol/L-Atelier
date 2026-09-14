# Batch Correction for Count Data

Removes batch effects from count data using ComBat or limma. Preserves
biological variation while removing technical batch effects.

For mass spec data (`massspec_data` objects, already log2-scale with
NAs), use
[`POSEIDON_correct_batch_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_correct_batch_massspec.md)
instead – it does not share this function's count-data assumptions (raw
counts, log2(x+1) transform, no missingness).

## Usage

``` r
POSEIDON_correct_batch_counts(
  quant_result,
  batch_col = "batch",
  group_col = "group",
  method = "combat",
  log_transform = TRUE,
  verbose = TRUE
)
```

## Arguments

- quant_result:

  A list containing:

  - counts: Matrix of counts (features x samples)

  - targets: Data.frame with sample metadata including batch information

- batch_col:

  Character string. Column name in targets containing batch information
  (default = "batch").

- group_col:

  Character string. Optional column name for biological groups to
  preserve during correction (default = "group").

- method:

  Character string. Batch correction method: "combat" (default) or
  "limma".

- log_transform:

  Logical. Should data be log-transformed before ComBat? ComBat assumes
  approximately normal data. If TRUE, applies log2(counts + 1) before
  correction and back-transforms after (default = TRUE).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

The input list with batch-corrected counts matrix.

## Details

Two methods are available:

**ComBat** (sva package):

- Empirical Bayes framework for batch correction

- Recommended for most applications

- Can preserve biological groups during correction

- Assumes approximately normal data (use log_transform = TRUE for
  counts)

**limma::removeBatchEffect**:

- Linear model-based batch removal

- Faster than ComBat

- Also preserves biological groups

- Works well with log-transformed data

For count data (ATAC-seq, RNA-seq), log transformation is recommended
before batch correction since these methods assume approximately normal
distributions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic batch correction
counts_corrected <- POSEIDON_correct_batch_counts(atac_counts)

# Preserve group differences while correcting batch
counts_corrected <- POSEIDON_correct_batch_counts(atac_counts,
                                            batch_col = "batch",
                                            group_col = "group")

# Use limma instead of ComBat
counts_corrected <- POSEIDON_correct_batch_counts(atac_counts, method = "limma")

} # }
```
