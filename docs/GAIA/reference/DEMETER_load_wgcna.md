# Load WGCNA results from export directory

Loads previously exported WGCNA results. Prioritizes RDS files (complete
objects) and falls back to CSV reconstruction if RDS not available.

## Usage

``` r
DEMETER_load_wgcna(
  results_dir,
  prefix = "wgcna",
  what = "all",
  strip_version = FALSE,
  verbose = TRUE
)
```

## Arguments

- results_dir:

  Path to the WGCNA export directory (output of
  ELEUTHIA_export_wgcna_results).

- prefix:

  Filename prefix used during export. Default: "wgcna".

- what:

  Character vector specifying which results to load. Options: "modules",
  "trait_cor", "gene_sig", "hubs", "enrichment", "all". Default: "all".

- strip_version:

  Logical. If TRUE, strips Ensembl version suffixes (e.g., ".12") from
  gene IDs. Default: FALSE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A named list containing the requested result objects. Each element will
be the appropriate class (e.g., wgcna_modules, wgcna_trait_cor) if
loaded from RDS, or a reconstructed object if loaded from CSV.

## Details

The function first looks for RDS files (e.g., wgcna_modules.rds) which
contain the complete R objects with all attributes and class
information. If RDS files are not found, it attempts to reconstruct the
objects from CSV files.

For full compatibility with ARTEMIS functions (e.g.,
ARTEMIS_wgcna_module_traits), ensure that `save_rds = TRUE` was used
during export.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load all WGCNA results
results <- DEMETER_load_wgcna("results/wgcna_output")

# Use modules with decoupleR integration
tf_trait_cor <- ARTEMIS_wgcna_module_traits(
  modules = results$modules,
  traits = t(tf_activities)
)

# Load only modules
results <- DEMETER_load_wgcna("results/wgcna_output", what = "modules")

} # }
```
