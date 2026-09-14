# Create a sample sheet template

Generates a CSV template with all required columns for a HORIZON
pipeline run. Fill in the paths and metadata, then pass to
[`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md)
before running any pipeline steps.

## Usage

``` r
HORIZON_create_sample_sheet(output_path = NULL, n_samples = 1)
```

## Arguments

- output_path:

  Character or NULL. Path to save the template CSV. If NULL, the
  data.frame is returned without writing to disk. Default NULL.

- n_samples:

  Integer. Number of placeholder rows to include. Default 1.

## Value

A data.frame with the sample sheet template, invisibly.

## Details

Required columns:

- sample_id:

  Unique sample identifier. Used to name output files and folders.

- fastq_r1:

  Absolute path to R1 (or single-end) FASTQ file.

- fastq_r2:

  Absolute path to R2 FASTQ file. Leave empty for single-end.

- paired_end:

  Logical. TRUE for paired-end libraries.

- output_dir:

  Absolute path to the root output directory for this sample. A
  sub-folder named `sample_id` will be created inside it.

- strandedness:

  Library strandedness: `"unstranded"`, `"forward"`, or `"reverse"`.

- condition:

  Experimental condition or group label. Passed through to the
  aggregated count matrix as metadata.

Additional columns can be added freely; they are carried through as
metadata.
