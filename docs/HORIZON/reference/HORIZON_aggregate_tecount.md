# Aggregate per-sample TEcount output into a TE count matrix

Reads the per-sample `.cntTable` files produced by
[`HORIZON_run_tecount`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_tecount.md)
and combines the **TE rows only** (gene rows are dropped – see
[`HORIZON_run_tecount`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_tecount.md)
for why) into a single TE-subfamily x sample count matrix.

## Usage

``` r
HORIZON_aggregate_tecount(
  sample_sheet,
  te_gtf,
  output_root = NULL,
  save_rds = TRUE,
  save_metadata = TRUE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).

- te_gtf:

  Character. Path to the same TE GTF passed to
  [`HORIZON_run_tecount`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_tecount.md)
  for every sample being aggregated.

- output_root:

  Character or NULL. Root output directory. If NULL, inferred from the
  first row's `output_dir`.

- save_rds:

  Logical. Also save the matrix as an RDS file. Default TRUE.

- save_metadata:

  Logical. Save the sample sheet as `sample_metadata.csv` alongside the
  matrix. Default TRUE.

## Value

The aggregated TE count matrix (TE subfamily x sample), invisibly.

## Details

TE rows are identified by cross-referencing each `.cntTable`'s row IDs
against the ground-truth set of composite TE IDs parsed directly from
`te_gtf`, rather than assuming a fixed ID-format/delimiter convention.

Output is written to `<output_root>/aggregated/tecount/`:

- `te_count_matrix.csv` / `.rds` – TE subfamily x sample count matrix

- `sample_metadata.csv` – sample sheet columns carried through (if
  `save_metadata = TRUE`)
