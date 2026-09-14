# Generate spliced/unspliced velocity matrices from PARSE split-pipe output

Calls `inst/python/parse_velocity.py` to build per-run spliced and
unspliced sparse matrices from the `tscp_assignment.csv` files produced
by
[`HORIZON_run_splitpipe`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_splitpipe.md),
concatenates them across all runs, and writes a single raw
`adata_vel.h5ad`. No cell/gene filtering or metadata addition is
performed — use
[`HORIZON_parse_DGE_filter`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_parse_DGE_filter.md)
for that step.

## Usage

``` r
HORIZON_run_parse_velocity(
  sample_sheet,
  conda_env,
  output_file = NULL,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_parse_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_parse_sheet.md).
  All rows must share the same `output_dir`.

- conda_env:

  Character. Path to the conda environment containing the required
  Python packages.

- output_file:

  Character or NULL. Path for the raw `adata_vel.h5ad`. Defaults to
  `<output_dir>/adata_vel.h5ad`.

- force:

  Logical. If `FALSE` (default), skip processing when `output_file`
  already exists. Set `TRUE` to reprocess.

## Value

Character. Path to the written `adata_vel.h5ad`, invisibly.

## Details

**Prerequisites:**

- All runs in `sample_sheet` must have been processed by
  [`HORIZON_run_splitpipe`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_splitpipe.md)
  and combined with
  [`HORIZON_combine_splitpipe`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_combine_splitpipe.md).

- `tscp_assignment.csv` files may be gzip-compressed (`.csv.gz`); this
  function will decompress them automatically in-place before
  processing.

- Required Python packages in `conda_env`: `scanpy`, `scvelo`,
  `anndata`, `dask`, `pandas`, `scipy`, `numpy`.
