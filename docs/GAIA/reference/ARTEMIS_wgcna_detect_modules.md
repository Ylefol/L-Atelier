# Detect co-expression modules

Identifies modules of co-expressed genes using WGCNA's blockwiseModules.
This is the core network construction and module detection step.

## Usage

``` r
ARTEMIS_wgcna_detect_modules(
  wgcna_data,
  power,
  network_type = "unsigned",
  tom_type = "unsigned",
  min_module_size = 30,
  merge_cut_height = 0.25,
  reassign_threshold = 0,
  max_block_size = 20000,
  n_threads = 1,
  save_tom = FALSE,
  tom_file_base = "TOM",
  verbose = TRUE
)
```

## Arguments

- wgcna_data:

  A wgcna_data object, wgcna_cluster object, or expression matrix
  (samples as rows, genes as columns).

- power:

  Soft threshold power. Can be:

  - Numeric value (e.g., 6)

  - A wgcna_power object from ARTEMIS_wgcna_pick_power()

  - NULL to auto-select (calls ARTEMIS_wgcna_pick_power internally)

- network_type:

  Network type: "unsigned", "signed", or "signed hybrid". Default:
  "unsigned".

- tom_type:

  TOM type: "unsigned", "signed", "signed Nowick", "none". Default:
  "unsigned".

- min_module_size:

  Minimum number of genes in a module. Default: 30.

- merge_cut_height:

  Dendrogram cut height for merging similar modules. Lower values = more
  merging. Default: 0.25.

- reassign_threshold:

  Threshold for reassigning genes to modules. Default: 0 (no
  reassignment).

- max_block_size:

  Maximum block size for network calculation. Increase for more genes
  but requires more memory. Default: 20000.

- n_threads:

  Number of threads for parallel computation. Default: 1.

- save_tom:

  Logical. Save TOM matrix (large, but needed for some analyses).
  Default: FALSE.

- tom_file_base:

  Base name for TOM file if save_tom = TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list with class "wgcna_modules" containing:

- module_names:

  Named vector: gene -\> module name (module_0, module_1, ...)

- module_labels:

  Named vector: gene -\> module number (0, 1, 2, ...)

- module_colors:

  Named vector: module_name -\> hex display color

- color_names:

  Named vector: module_name -\> WGCNA color name (for reference)

- module_eigengenes:

  Data.frame of module eigengenes (MEmodule_1, ...)

- module_summary:

  Data.frame: module, n_genes, color_name

- gene_module_df:

  Data.frame: gene, module_name, module_label, color_name

- dendrograms:

  List of gene dendrograms

- net:

  Raw blockwiseModules output (for advanced use)

- power:

  Power used

- datExpr:

  Expression matrix used

- n_modules:

  Number of modules detected (excluding module_0)

## Details

The function wraps WGCNA::blockwiseModules() with sensible defaults and
returns results in a more accessible format.

Module colors are assigned by WGCNA. The "grey" module contains genes
that couldn't be assigned to any module.

## Examples

``` r
if (FALSE) { # \dontrun{
wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
power_result <- ARTEMIS_wgcna_pick_power(wgcna_data)
modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = power_result)

# Or with auto power selection
modules <- ARTEMIS_wgcna_detect_modules(wgcna_data, power = NULL)

} # }
```
