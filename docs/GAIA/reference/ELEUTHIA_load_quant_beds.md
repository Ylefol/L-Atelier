# Load Quantification BED Files from Sample Sheet

Loads fragment BED files for quantification from the bed_loc column of a
validated sample sheet.

## Usage

``` r
ELEUTHIA_load_quant_beds(sample_sheet, omics, verbose = TRUE)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame with bed_loc column.

- omics:

  Character string. Omics type to load ("ATACseq" or "CHIPseq").

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A named list of BED data.frames, with names corresponding to sample_id.

## Details

This function reads fragment BED files (typically created with bedtools
bamtobed) for use with ELEUTHIA_quantify_bed(). It uses the bed_loc
column from the sample sheet rather than file_loc/file_name.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
atac_frags <- ELEUTHIA_load_quant_beds(result$sample_sheet, "ATACseq")
chip_frags <- ELEUTHIA_load_quant_beds(result$sample_sheet, "CHIPseq")

} # }
```
