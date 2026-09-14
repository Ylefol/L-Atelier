# Combine multiple PARSE split-pipe runs

Wraps the `split-pipe --mode comb` step to merge all per-run outputs
produced by
[`HORIZON_run_splitpipe`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_splitpipe.md)
into a single combined dataset. This step is required before running
[`HORIZON_run_parse_velocity`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_parse_velocity.md)
and before loading data into CAULDRON via `TALARIA_load_parse()`.

## Usage

``` r
HORIZON_combine_splitpipe(
  sample_sheet,
  conda_env,
  output_subdir = "combined",
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_parse_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_parse_sheet.md).
  All rows must share the same `output_dir`. Each run's split-pipe
  output directory (`output_dir/run_id/`) is passed as a sublibrary to
  split-pipe.

- conda_env:

  Character. Path to the conda environment containing `split-pipe` — the
  same environment used for
  [`HORIZON_run_splitpipe`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_splitpipe.md).

- output_subdir:

  Character. Name of the subdirectory within `output_dir` where the
  combined output will be written. Default `"combined"`. The resulting
  path (`output_dir/combined/all-sample/`) is what
  [`HORIZON_run_parse_velocity`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_parse_velocity.md)
  expects for `cell_metadata_path`.

- force:

  Logical. If `FALSE` (default), skip if the combined output directory
  already exists. Set `TRUE` to reprocess.

- verbose:

  Logical. If `TRUE` (default), prints the full split-pipe command
  before executing.

## Value

Character. Path to the combined output directory, invisibly.
