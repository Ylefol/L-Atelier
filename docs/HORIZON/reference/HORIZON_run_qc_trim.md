# Run QC and adapter trimming with fastp

Wraps [`rfastp`](https://rdrr.io/pkg/Rfastp/man/rfastp.html) to perform
quality control and adapter trimming in a single pass. fastp
auto-detects adapters by default: for paired-end data it uses read-pair
overlap analysis; for single-end data it falls back to a built-in list
of common Illumina adapters. Explicit adapter sequences can be provided
to override auto-detection.

## Usage

``` r
HORIZON_run_qc_trim(
  sample_sheet,
  sample_id,
  threads = 4,
  adapter_r1 = NULL,
  adapter_r2 = NULL,
  force = FALSE,
  ...
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).

- sample_id:

  Character. Sample ID to process. Must match a value in the `sample_id`
  column of `sample_sheet`.

- threads:

  Integer. Number of threads for fastp. Default 4.

- adapter_r1:

  Character or NULL. Explicit R1 adapter sequence. Default NULL
  (auto-detect via `adapterSequenceRead1 = "auto"`).

- adapter_r2:

  Character or NULL. Explicit R2 adapter sequence. Default NULL
  (auto-detect). Ignored for single-end samples.

- force:

  Logical. If `FALSE` (default), skip trimming when the output FASTQs
  already exist and return their paths directly. Set to `TRUE` to re-run
  and overwrite existing files.

- ...:

  Additional arguments passed directly to
  [`rfastp`](https://rdrr.io/pkg/Rfastp/man/rfastp.html) (e.g.,
  `qualityFiltering`, `qualityFilterPhred`, `lengthFiltering`,
  `minReadLength`).

## Value

Named list with elements:

- trimmed_r1:

  Path to trimmed R1 FASTQ.

- trimmed_r2:

  Path to trimmed R2 FASTQ, or NULL for single-end.

- report_html:

  Path to HTML QC report.

- report_json:

  Path to JSON QC report.

- json_data:

  The JSON report object returned by rfastp, or NULL when skipped.

## Details

Outputs are written to `<output_dir>/<sample_id>/qc/`. Rfastp uses
`outputFastq` as a path prefix; trimmed FASTQs are named
`<prefix>_R1.fastq.gz` / `<prefix>_R2.fastq.gz` and QC reports are
`<prefix>.html` / `<prefix>.json`.
