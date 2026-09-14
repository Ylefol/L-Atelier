# Aggregate per-sample Salmon quant.sf files into TPM / count matrices

Reads the per-sample `quant.sf` files produced by
[`HORIZON_run_salmon`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_salmon.md)
and combines them into transcript x sample matrices of TPM and estimated
read counts (`NumReads`). Samples missing a `quant.sf` file are skipped
with a warning.

Output is written to `<output_root>/aggregated/salmon/`:

- `tpm_matrix.csv` / `.rds` — transcript x sample TPM

- `counts_matrix.csv` / `.rds` — transcript x sample estimated counts
  (`NumReads`)

- `sample_metadata.csv` — sample sheet columns carried through (if
  `save_metadata = TRUE`)

## Usage

``` r
HORIZON_aggregate_salmon(
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

  Character or NULL. Root output directory. If NULL, inferred from the
  first row's `output_dir`.

- save_rds:

  Logical. Also save matrices as RDS. Default TRUE.

- save_metadata:

  Logical. Save sample sheet as `sample_metadata.csv` alongside the
  matrices. Default TRUE.

## Value

A list with elements `tpm` and `counts` (transcript x sample matrices),
invisibly.
