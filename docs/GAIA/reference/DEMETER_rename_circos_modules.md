# Rename and recolor WGCNA modules in circos data

Rename and recolor WGCNA modules in circos data

## Usage

``` r
DEMETER_rename_circos_modules(
  circos_data,
  rename_vector,
  recolor_vector = NULL
)
```

## Arguments

- circos_data:

  A `part_wgcna_circos` object.

- rename_vector:

  Named character vector: old names as names, new names as values.

- recolor_vector:

  Named character vector: new names as names, hex colors as values. If
  NULL, colors are auto-generated for renamed modules.

## Value

Updated `part_wgcna_circos` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Rename WGCNA color-based names to generic labels
rename_vec <- c(blue = "module_1", brown = "module_2", turquoise = "module_3",
                dummy = "dummy")
recolor_vec <- c(module_1 = "dodgerblue2", module_2 = "#E31A1C",
                 module_3 = "green4", dummy = "#CCCCCC")
circos_data <- DEMETER_rename_circos_modules(circos_data, rename_vec, recolor_vec)

} # }
```
