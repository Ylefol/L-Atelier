# Create a PARSE Biosciences sample sheet template

Generates a CSV template for PARSE split-pipe runs. Each row represents
one sequencing run (library). Fill in paths and metadata, then pass to
[`HORIZON_validate_parse_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_parse_sheet.md)
before running
[`HORIZON_run_splitpipe`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_splitpipe.md).

## Usage

``` r
HORIZON_create_parse_sheet(output_path = NULL, n_runs = 1)
```

## Arguments

- output_path:

  Character or NULL. Path to write the template CSV. If NULL, the
  data.frame is returned without writing. Default NULL.

- n_runs:

  Integer. Number of placeholder rows to include. Default 1.

## Value

A data.frame with the template, invisibly.

## Details

Required columns:

- run_id:

  Unique run identifier. Used to name the output sub-folder.

- fastq_r1:

  Absolute path to R1 FASTQ file.

- fastq_r2:

  Absolute path to R2 FASTQ file.

- output_dir:

  Absolute path to root output directory. A sub-folder named `run_id`
  will be created inside it.

- genome_dir:

  Absolute path to the split-pipe genome directory (built with
  `split-pipe --mode mkref`).

- chemistry:

  Library chemistry version: `"v2"`, `"v3"`, or `"v4"`.

- kit:

  Library kit: `"WT"`, `"WT_mini"`, `"WT_mega"`, `"WT_mega_384"`,
  `"WT_penta"`, or `"WT_penta_384"`.

Additional columns can be added freely; they are carried through as
metadata.
