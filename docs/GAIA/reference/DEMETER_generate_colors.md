# Generate distinct visible colors

Generates a set of distinct colors suitable for plotting, filtering out
near-white colors that would be invisible on light backgrounds. Useful
for remapping WGCNA module colors or any categorical color assignment.

## Usage

``` r
DEMETER_generate_colors(n, palette = "Dark 3", max_brightness = 200)
```

## Arguments

- n:

  Integer. Number of colors needed.

- palette:

  Character. HCL color palette name. Default: "Dark 3". See
  `grDevices::hcl.pals("qualitative")` for options.

- max_brightness:

  Numeric (0-255). Maximum mean RGB brightness allowed. Colors brighter
  than this are excluded. Default: 200.

## Value

Character vector of `n` hex color codes.

## Examples

``` r
if (FALSE) { # \dontrun{
cols <- DEMETER_generate_colors(10)
cols <- DEMETER_generate_colors(5, palette = "Set 2")

# Remap WGCNA module colors for circos
mod_names <- rownames(circos_data$module_df)
mod_names <- mod_names[mod_names != "dummy"]
new_colors <- setNames(DEMETER_generate_colors(length(mod_names)), mod_names)
identity_rename <- setNames(mod_names, mod_names)
circos_data <- DEMETER_rename_circos_modules(circos_data, identity_rename, new_colors)

} # }
```
