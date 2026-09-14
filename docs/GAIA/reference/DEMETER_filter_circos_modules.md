# Filter circos data to specific WGCNA modules

Keeps selected modules and assigns all other genes to a "dummy" module.
The dummy module appears in the circos but is visually muted, allowing
focus on modules of interest while preserving all PART cluster data.

## Usage

``` r
DEMETER_filter_circos_modules(
  circos_data,
  selected_modules,
  dummy_name = "dummy",
  dummy_color = "#CCCCCC"
)
```

## Arguments

- circos_data:

  A `part_wgcna_circos` object.

- selected_modules:

  Character vector of module names to keep.

- dummy_name:

  Character. Name for the catch-all module. Default: "dummy".

- dummy_color:

  Character. Color for the dummy module. Default: "#CCCCCC".

## Value

Updated `part_wgcna_circos` object.
