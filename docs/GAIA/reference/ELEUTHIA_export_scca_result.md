# Export sCCA Direct Run Results

Exports results from a direct (non-CV) sCCA run using
[`HEPHAESTUS_multi_convCCA()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_multi_convCCA.md)
or
[`HEPHAESTUS_multi_relPMDCCA()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_multi_relPMDCCA.md).
Saves per-dataset feature weight tables, a sparsity summary, the feature
weight plot, a metadata file, and optionally the full R object as RDS.

## Usage

``` r
ELEUTHIA_export_scca_result(
  result,
  output_dir,
  dataset_names = NULL,
  method_name = "sCCA",
  prefix = "scca",
  top_n = NULL,
  save_plots = TRUE,
  plot_format = "png",
  save_rds = TRUE,
  verbose = TRUE
)
```

## Arguments

- result:

  List with element `$W`: a list of canonical weight vectors, one per
  dataset. Direct output of
  [`HEPHAESTUS_multi_convCCA()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_multi_convCCA.md)
  or
  [`HEPHAESTUS_multi_relPMDCCA()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_multi_relPMDCCA.md).

- output_dir:

  Character. Directory to save results. Created if needed.

- dataset_names:

  Character vector. Names for each dataset (e.g.,
  `c("ATAC", "ChIP", "RNA")`). Default: "Dataset_1", "Dataset_2", etc.

- method_name:

  Character. Method label used in plot titles and metadata. Default:
  "sCCA".

- prefix:

  Character. Prefix for all output filenames. Default: "scca".

- top_n:

  Integer or NULL. If NULL, exports all non-zero weighted features per
  dataset. If an integer, exports the top N features by absolute weight
  (regardless of zero status). Default: NULL.

- save_plots:

  Logical. Save the feature weight plot. Default: TRUE.

- plot_format:

  Character. Plot format: "png", "pdf", or "both". Default: "png".

- save_rds:

  Logical. Save full result object as RDS. Default: TRUE.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Invisibly returns a character vector of file paths created.

## Details

Creates the following files in `output_dir`:

    output_dir/
      <prefix>_<dataset_name>_weights.csv   (one per dataset)
      <prefix>_sparsity_summary.csv
      <prefix>_feature_weights.png / .pdf
      <prefix>_metadata.txt
      <prefix>_result.rds                   (optional)

Weight tables contain: `feature`, `weight`, `abs_weight`, `rank` (by
absolute weight). Feature names are taken from `names(result$W[[i]])` if
available, otherwise integer indices are used.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- HEPHAESTUS_multi_convCCA(X = scca_data, tau = as.list(rep(0.1, 4)))
ELEUTHIA_export_scca_result(
  result,
  output_dir  = "results/scca/",
  dataset_names = c("ATAC_1", "ATAC_2", "ChIP", "RNA"),
  method_name = "multi.convCCA"
)
} # }
```
