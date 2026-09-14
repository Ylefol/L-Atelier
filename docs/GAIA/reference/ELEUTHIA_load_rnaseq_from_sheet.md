# Load RNA-seq Counts from Sample Sheet

Loads all RNA-seq count files from a validated sample sheet and combines
them into a count matrix.

## Usage

``` r
ELEUTHIA_load_rnaseq_from_sheet(
  sample_sheet,
  omics = "RNAseq",
  header = NA,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame.

- omics:

  Character string. Omics type to load (default = "RNAseq").

- header:

  Logical or NA. Passed through to
  [`ELEUTHIA_load_count_file`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_count_file.md)
  for each sample. NA (default) auto-detects per file.

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A list containing:

- counts:

  Matrix of counts (genes x samples)

- targets:

  Data.frame with sample metadata

## Details

All count files must have the same genes (rows). The function will use
the intersection of genes if there are differences, with a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
rna_data <- ELEUTHIA_load_rnaseq_from_sheet(result$sample_sheet)
count_matrix <- rna_data$counts

} # }
```
