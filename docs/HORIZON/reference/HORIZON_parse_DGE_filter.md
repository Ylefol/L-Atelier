# Filter velocity AnnData using PARSE DGE output and add metadata

Takes the raw `adata_vel.h5ad` from
[`HORIZON_run_parse_velocity`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_parse_velocity.md)
and filters it to the cells and genes present in the split-pipe DGE
output (typically `DGE_filtered`), then adds cell and gene metadata from
`cell_metadata.csv` and `all_genes.csv`.

## Usage

``` r
HORIZON_parse_DGE_filter(
  velocity_h5ad,
  dge_dir,
  conda_env,
  output_file = NULL,
  force = FALSE
)
```

## Arguments

- velocity_h5ad:

  Character. Path to the raw `adata_vel.h5ad` produced by
  [`HORIZON_run_parse_velocity`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_parse_velocity.md).

- dge_dir:

  Character. Path to the split-pipe DGE directory (e.g.
  `"<output_dir>/combined/all-sample/DGE_filtered"`). Must contain
  `cell_metadata.csv` and `all_genes.csv`.

- conda_env:

  Character. Path to the conda environment used for
  [`HORIZON_run_parse_velocity`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_parse_velocity.md).

- output_file:

  Character or NULL. Output path for the filtered h5ad. Defaults to
  `adata_vel_filtered.h5ad` in the same directory as `velocity_h5ad`.

- force:

  Logical. Reprocess even if `output_file` already exists.

## Value

Character. Path to the filtered h5ad, invisibly.
