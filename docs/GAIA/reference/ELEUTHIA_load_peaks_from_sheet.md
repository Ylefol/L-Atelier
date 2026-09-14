# Load Multiple Peak Files from Sample Sheet

Loads all peak files for a given omics type from a validated sample
sheet into a named list.

## Usage

``` r
ELEUTHIA_load_peaks_from_sheet(
  sample_sheet,
  omics,
  min_score = 0,
  min_qvalue = 0,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame with sample_id column.

- omics:

  Character string. Omics type to load ("ATACseq" or "CHIPseq").

- min_score:

  Numeric. Minimum peak score to retain (default = 0).

- min_qvalue:

  Numeric. Minimum -log10(qvalue) to retain (default = 0).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A named list of peak data.frames, with names corresponding to sample_id.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
atac_peaks <- ELEUTHIA_load_peaks_from_sheet(result$sample_sheet, "ATACseq")

} # }
```
