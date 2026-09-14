# Generate WGCNA analysis report

Creates a formatted text report summarizing the WGCNA analysis results.

## Usage

``` r
ELEUTHIA_wgcna_report(
  output_dir,
  wgcna_data = NULL,
  power_result = NULL,
  modules = NULL,
  trait_cor = NULL,
  gene_sig = NULL,
  hubs = NULL,
  prefix = "wgcna",
  verbose = TRUE
)
```

## Arguments

- output_dir:

  Output directory for report.

- wgcna_data:

  Optional. wgcna_data object (data summary).

- power_result:

  Optional. wgcna_power object (power selection).

- modules:

  Optional. wgcna_modules object (module detection).

- trait_cor:

  Optional. wgcna_trait_cor object (trait correlations).

- gene_sig:

  Optional. wgcna_gene_sig object (gene significance).

- hubs:

  Optional. wgcna_hubs object (hub genes).

- prefix:

  Filename prefix. Default: "wgcna".

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisible path to report file.
