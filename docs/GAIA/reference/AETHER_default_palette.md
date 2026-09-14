# Default Color Palettes

Returns default color palettes for common omics types and experimental
groups. These can be used as base colors for gradient generation.

## Usage

``` r
AETHER_default_palette(type = "omics")
```

## Arguments

- type:

  Character string. Type of palette: "omics" = colors by omics type
  (ATACseq, CHIPseq, RNAseq), "groups" = colors by common group names
  (WT, KO, Control, Treatment).

## Value

Named character vector of hex color codes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get omics palette
omics_colors <- AETHER_default_palette("omics")
# Returns: ATACseq="#2ca02c", CHIPseq="#1f77b4", RNAseq="#ff7f0e"

# Get group palette
group_colors <- AETHER_default_palette("groups")

} # }
```
