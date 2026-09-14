# Load multiple HOMER result directories into a batch object

Scans a base directory for HOMER output subdirectories and loads them
into a combined homer_motif_batch object.

## Usage

``` r
APOLLO_load_homer_batch(output_dir, set_names = NULL, verbose = TRUE)
```

## Arguments

- output_dir:

  Base directory containing subdirectories per peak set.

- set_names:

  Optional character vector of subdirectory names to load. If NULL,
  auto-detects from subdirectories containing knownResults.txt.

- verbose:

  Logical. Print loading messages. Default TRUE.

## Value

S3 object of class "homer_motif_batch" (same as batch enrichment).
