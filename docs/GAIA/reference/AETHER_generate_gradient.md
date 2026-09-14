# Generate Color Gradient

Generates a gradient of colors from a dark base color to lighter shades.
Useful for plotting multiple samples within the same group where visual
distinction is needed but with a consistent color theme.

## Usage

``` r
AETHER_generate_gradient(base_color, n = 5, direction = "to_light")
```

## Arguments

- base_color:

  Character string. The base (darkest) color as a hex code (e.g.,
  "#08519c") or R color name (e.g., "darkblue").

- n:

  Integer. Number of colors to generate (default = 5).

- direction:

  Character string. Direction of gradient: "to_light" (default) = dark
  to light, "to_dark" = light to dark.

## Value

Character vector of n hex color codes.

## Details

The function converts the base color to HSL (Hue, Saturation, Lightness)
color space, then generates n colors by interpolating the lightness
value while preserving hue and saturation.

For circos plots with multiple samples per group:

- Use "to_light" direction

- Plot lightest color first, darkest last (dark on top)

- This creates clean visuals where similar data is "hidden" under darker
  lines

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate 6 shades of blue
blues <- AETHER_generate_gradient("#08519c", n = 6)

# Generate 4 shades of red, light to dark
reds <- AETHER_generate_gradient("#d62728", n = 4, direction = "to_dark")

# Preview colors
barplot(rep(1, 6), col = blues, border = NA)

} # }
```
