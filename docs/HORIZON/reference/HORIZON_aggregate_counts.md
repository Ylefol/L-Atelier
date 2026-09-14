# Aggregate per-sample count files into a single count matrix

Reads the per-sample count files produced by
[`HORIZON_run_count`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_count.md)
and combines them into a single gene x sample count matrix. Samples
missing a count file are skipped with a warning.

## Usage

``` r
HORIZON_aggregate_counts(
  sample_sheet,
  output_root = NULL,
  save_rds = TRUE,
  save_metadata = TRUE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).

- output_root:

  Character or NULL. Root output directory — the parent of the
  per-sample folders. If NULL, inferred from the first row's
  `output_dir`. Provide explicitly if samples have different
  `output_dir` values and you want the aggregated folder somewhere
  specific.

- save_rds:

  Logical. Save the count matrix as an RDS file in addition to CSV.
  Default TRUE.

- save_metadata:

  Logical. Save the sample sheet as `sample_metadata.csv` alongside the
  count matrix. Default TRUE.

## Value

data.frame. The aggregated count matrix (genes x samples), invisibly.

## Details

Output is written to `<output_root>/aggregated/counts/`:

- `count_matrix.csv` — genes x samples count matrix

- `count_matrix.rds` — same matrix as an R object (if `save_rds = TRUE`)

- `sample_metadata.csv` — sample sheet columns carried through as
  metadata (if `save_metadata = TRUE`)
