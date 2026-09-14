# Get Subset of Sample Sheet by Omics Type

Extracts samples of a specific omics type from a validated sample sheet.

## Usage

``` r
ELEUTHIA_get_omics_subset(sample_sheet, omics)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame (with sample_id column).

- omics:

  Character string. The omics type to extract (e.g., "ATACseq",
  "CHIPseq", "RNAseq", "5hmu", or any custom omics label).

## Value

A data.frame containing only rows matching the specified omics type.
Returns empty data.frame if no matches found.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
atac_samples <- ELEUTHIA_get_omics_subset(result$sample_sheet, "ATACseq")
custom_samples <- ELEUTHIA_get_omics_subset(result$sample_sheet, "5hmu")

} # }
```
