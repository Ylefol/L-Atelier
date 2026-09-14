# RUVg Correction Using Control Genes

Removes unwanted (technical) variation from a counts matrix via RUVSeq's
RUVg method, which estimates factors of unwanted variation from a
user-supplied set of negative control genes assumed to be constant
across the samples' biological conditions (e.g. housekeeping genes).

## Usage

``` r
POSEIDON_correct_ruvg(quant_result, control_genes, k = 1, verbose = TRUE)
```

## Arguments

- quant_result:

  A list containing:

  - counts: Matrix of raw counts (features x samples)

  - targets: Data.frame with sample metadata

- control_genes:

  Character vector of feature IDs (matching
  `rownames(quant_result$counts)`) to use as negative controls. RUVg
  does not ship or assume any particular list – housekeeping-gene panels
  or empirically-derived (least-DE) controls must be supplied by the
  caller.

- k:

  Integer. Number of factors of unwanted variation to estimate (default
  = 1).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

The input list with RUVg-corrected counts, plus:

- ruvg_W: matrix of estimated unwanted-variation factors (samples x k)

- ruvg_control_genes: the control genes actually used (i.e. found in
  `quant_result$counts`)

## Details

Requires the `RUVSeq` Bioconductor package:
`BiocManager::install("RUVSeq")`.

## Examples

``` r
if (FALSE) { # \dontrun{
housekeeping <- c("ENSG00000075624", "ENSG00000111640")  # ACTB, GAPDH, ...
combined_data <- POSEIDON_correct_ruvg(combined_data, control_genes = housekeeping)
} # }
```
