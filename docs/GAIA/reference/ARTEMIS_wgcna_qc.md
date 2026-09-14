# Quality control for WGCNA data

Checks for genes and samples with too many missing values or zero
variance. Wraps WGCNA's goodSamplesGenes() with additional reporting.

## Usage

``` r
ARTEMIS_wgcna_qc(wgcna_data, verbose = TRUE)
```

## Arguments

- wgcna_data:

  A wgcna_data object from ARTEMIS_wgcna_prepare(), or a matrix with
  samples as rows and genes as columns.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

If wgcna_data object provided, returns updated object with QC applied.
If matrix provided, returns filtered matrix. In both cases, attributes
record what was removed.

## Details

Uses WGCNA::goodSamplesGenes() which checks for:

- Genes with too many missing values

- Samples with too many missing values

- Genes with zero variance

## Examples

``` r
if (FALSE) { # \dontrun{
wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
wgcna_data <- ARTEMIS_wgcna_qc(wgcna_data)

} # }
```
