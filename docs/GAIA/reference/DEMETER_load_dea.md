# Load DEA results from export directories

Loads previously exported DEA results from one or more experiment
subdirectories. Prioritizes RDS files (complete objects) and falls back
to CSV reconstruction. Optionally loads associated PART results.

## Usage

``` r
DEMETER_load_dea(
  dea_dir,
  experiments = NULL,
  prefix = "dea",
  what = "dea",
  strip_version = FALSE,
  verbose = TRUE
)
```

## Arguments

- dea_dir:

  Path to the parent directory containing DEA result subdirectories (one
  per experiment/comparison).

- experiments:

  Character vector. Specific subdirectory names to load. If NULL
  (default), all subdirectories are scanned.

- prefix:

  Character. Filename prefix used during export. Default: "dea".

- what:

  Character vector specifying what to load per experiment:

  - "dea": the DEA result (default)

  - "part": the PART clustering result

  - "all": both dea and part

- strip_version:

  Logical. If TRUE, strips Ensembl version suffixes from gene IDs.
  Default: FALSE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

A named list (one element per experiment). Each element is a list
containing the requested objects:

- dea:

  DEA result list with: results (data.frame), summary, comparison,
  norm_counts. If loaded from RDS, includes dds object.

- part:

  artemis_part object (if requested and available)

If only "dea" is requested, each element is the DEA result directly (not
wrapped in a sub-list) for convenience.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load all DEA results
dea_list <- DEMETER_load_dea("results/DEA/")
# Returns: list(C5RO = <dea_result>, CS1AN = <dea_result>, ...)

# Pass directly to circos prep
circos_data <- DEMETER_prepare_part_wgcna_circos(
  part_result, wgcna$modules, wgcna$trait_cor,
  dea_results = dea_list
)

# Load specific experiments with PART results
results <- DEMETER_load_dea("results/DEA/",
                             experiments = c("C5RO", "CS1AN"),
                             what = "all")
# results$C5RO$dea, results$C5RO$part

} # }
```
