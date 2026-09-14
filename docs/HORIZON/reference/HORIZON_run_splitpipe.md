# Run PARSE Biosciences split-pipe

Wraps the `split-pipe` CLI to perform demultiplexing, alignment, and
cell/nucleus barcode counting for a single PARSE Biosciences run.
Supports both single-cell and single-nucleus RNA-seq data.

## Usage

``` r
HORIZON_run_splitpipe(
  sample_sheet,
  run_id,
  sample_layout,
  conda_env,
  mode = "all",
  threads = 8L,
  force = FALSE,
  kit_score_skip = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_parse_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_parse_sheet.md).
  The output location for each run is controlled by the `output_dir`
  column in this sheet — a sub-folder named `run_id` is created inside
  it. There is no separate `output_dir` parameter on this function.

- run_id:

  Character. Run ID to process (must match a row in `sample_sheet`).

- sample_layout:

  Named character vector. Maps sample names to split-pipe well range
  strings (e.g. `c(control = "A1-B6", treatment = "C1-D6")`).

- conda_env:

  Character. Path to the conda environment containing `split-pipe` (e.g.
  `"/path/to/conda/envs/parse_env"`).

- mode:

  Character. split-pipe run mode. Default `"all"` (demultiplex + align +
  count in one pass). Use `"split"` for demultiplexing only, or
  `"count"` to recount from a prior split output.

- threads:

  Integer. Number of threads passed to split-pipe (`--nthreads`).
  Default `8L`.

- force:

  Logical. If `FALSE` (default), skip processing when the output
  directory already exists. Set `TRUE` to reprocess and overwrite.

- kit_score_skip:

  Logical. If `TRUE`, passes `--kit_score_skip` to split-pipe to bypass
  the automatic kit barcode check. Useful when the detected kit does not
  exactly match the specified kit but processing should proceed anyway.
  Default `FALSE`.

- verbose:

  Logical. If `TRUE`, prints the full split-pipe command before
  executing it. Useful for verifying parameter passing. Default `TRUE`.

- ...:

  Additional split-pipe flags passed verbatim to the command.

## Value

Character. Path to the split-pipe output directory, invisibly.

## Details

`split-pipe` must be installed in a dedicated conda environment passed
via `conda_env`. It cannot share the standard `horizon_cli` environment
due to conflicting Python dependencies.

Example installation:

      conda create -n parse_env python=3.10
      conda activate parse_env
      pip install parsebiosciences

**Sample layout format:** A named character vector mapping sample names
to split-pipe well range strings. Names become the split-pipe sample
identifiers.

      # All wells as a single sample
      sample_layout = c(all_cells = "A1-D12")

      # Two conditions across separate well ranges
      sample_layout = c(control = "A1-B6", treatment = "C1-D6")

The split-pipe output directory (`output_dir/run_id/`) contains a
per-sample DGE matrix compatible with `TALARIA_load_parse()` in
CAULDRON.
