# Expand Genomic Regions by a Window Size

Symmetrically expands genomic regions by a specified window size.
Handles chromosome boundaries by clamping start to 0 and optionally
clamping end to chromosome length if sizes are provided.

## Usage

``` r
ELEUTHIA_expand_regions(
  regions,
  window_size = 0,
  chrom_sizes = NULL,
  verbose = TRUE
)
```

## Arguments

- regions:

  A data.frame with at minimum columns: chr, start, end. Typically from
  ELEUTHIA_create_consensus_peaks() or
  ELEUTHIA_call_regions_from_fragments().

- window_size:

  Integer. Number of base pairs to extend on each side (default = 0, no
  expansion). Must be \>= 0.

- chrom_sizes:

  Optional. Either:

  - A named numeric vector where names are chromosome names and values
    are chromosome lengths

  - A data.frame with columns 'chr' and 'size'

  - NULL (default) - end coordinates are not clamped to chromosome
    length

- verbose:

  Logical. Print summary messages (default = TRUE).

## Value

A data.frame with the same structure as input, but with start and end
coordinates expanded by window_size (respecting chromosome boundaries).
An additional column 'original_width' is added showing the pre-expansion
width, and 'window_size' records the expansion used.

## Details

For each region:

- new_start = max(0, start - window_size)

- new_end = end + window_size (or min(end + window_size, chrom_length)
  if sizes provided)

If window_size = 0, the function returns the regions unchanged (with
added metadata columns). This allows consistent usage in pipelines that
test multiple window sizes including zero.

## Examples

``` r
if (FALSE) { # \dontrun{
# No expansion (window_size = 0)
regions_0 <- ELEUTHIA_expand_regions(chip_regions, window_size = 0)

# Expand by 500bp on each side
regions_500 <- ELEUTHIA_expand_regions(chip_regions, window_size = 500)

# With chromosome size clamping
chrom_sizes <- c(chr1 = 248956422, chr2 = 242193529, chr3 = 198295559)
regions_500 <- ELEUTHIA_expand_regions(chip_regions, window_size = 500,
                                        chrom_sizes = chrom_sizes)

# Use in anchored analysis workflow
expanded <- ELEUTHIA_expand_regions(chip_regions, window_size = 250)
atac_at_anchors <- ELEUTHIA_quantify_bed(sample_sheet, expanded, "ATACseq")

} # }
```
