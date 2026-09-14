# Run HOMER motif enrichment on multiple BED files

Iterates over a named vector of BED file paths, running
APOLLO_homer_motif_enrichment() for each and combining results.

## Usage

``` r
APOLLO_homer_motif_enrichment_batch(
  bed_files,
  genome,
  output_dir,
  ...,
  verbose = TRUE
)
```

## Arguments

- bed_files:

  Named character vector of BED file paths. Names are used as peak set
  identifiers.

- genome:

  HOMER genome string (e.g., "hg38").

- output_dir:

  Base output directory. Subdirectories are created per peak set.

- ...:

  Additional arguments passed to APOLLO_homer_motif_enrichment() (size,
  mask, denovo, denovo_n, bg, nproc, extra_args).

- verbose:

  Logical. Print progress messages. Default TRUE.

## Value

S3 object of class "homer_motif_batch" with components:

- results:

  Named list of homer_motif objects

- combined:

  Data.frame of all known motif results with peak_set column

- summary:

  Data.frame with per-set summary statistics

- metadata:

  List with genome, parameters, bed_files

## Examples

``` r
if (FALSE) { # \dontrun{
bed_files <- c(
  up   = "results/peaks_up.bed",
  down = "results/peaks_down.bed"
)
batch <- APOLLO_homer_motif_enrichment_batch(
  bed_files, "hg38", "results/motifs/"
)

} # }
```
