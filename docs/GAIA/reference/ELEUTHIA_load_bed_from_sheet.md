# Load ChIP/CUT&TAG BED Files from Sample Sheet

Loads fragment BED files for ChIP-seq or CUT&TAG data from a validated
sample sheet.

## Usage

``` r
ELEUTHIA_load_bed_from_sheet(sample_sheet, omics = "CHIPseq", verbose = TRUE)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame.

- omics:

  Character string. Omics type to load (default = "CHIPseq").

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A named list of BED data.frames, with names corresponding to sample_id.
