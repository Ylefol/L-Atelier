# Convert Chromosome Sizes to Seqinfo Object

Converts a chromosome sizes data.frame to a Seqinfo object for use with
GenomicRanges and related packages.

## Usage

``` r
APOLLO_chr_sizes_to_seqinfo(chr_sizes, genome = NA_character_)
```

## Arguments

- chr_sizes:

  A data.frame with columns 'chr' and 'size', typically from
  APOLLO_get_chromosome_sizes().

- genome:

  Optional character string for genome name (e.g., "hg38", "T2T").

## Value

A Seqinfo object.
