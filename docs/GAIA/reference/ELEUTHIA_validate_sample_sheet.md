# Eleuthia - Sample Sheet Functions

Functions for validating and processing sample sheets used in
multi-omics analysis pipelines. Validate Sample Sheet

Validates a sample sheet for use in multi-omics analysis. Checks for
required columns, valid values, file existence, and generates sample IDs
programmatically.

## Usage

``` r
ELEUTHIA_validate_sample_sheet(
  sample_sheet,
  check_files = TRUE,
  check_beds = TRUE,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  Either a data.frame containing sample information, or a character
  string path to a CSV file.

- check_files:

  Logical. If TRUE, verifies that data files exist at specified paths
  (default = TRUE).

- check_beds:

  Logical. If TRUE, verifies that BED files for quantification exist
  where required by format (peaks, bed). Only checked if check_files =
  TRUE (default = TRUE).

- verbose:

  Logical. If TRUE, prints validation progress and warnings (default =
  TRUE).

## Value

A list containing:

- sample_sheet:

  The validated and cleaned sample sheet with added sample_id column

- valid:

  Logical. TRUE if validation passed with no errors

- errors:

  Character vector of error messages (empty if valid)

- warnings:

  Character vector of warning messages

## Details

The function operates in two modes depending on whether a `sample_id`
column is present:

**Experimental mode** (no `sample_id` column): all of `bio_rep`,
`tech_rep`, and `batch` are required, and `sample_id` is auto-generated
as `<omics>_<group>_<bio_rep><tech_rep>_b<batch>`.

**Cohort mode** (`sample_id` column present): `bio_rep`, `tech_rep`, and
`batch` are optional. The provided `sample_id` values are used directly
and validated for uniqueness. Suitable for public datasets (e.g., TCGA)
where each sample is a unique individual rather than a replicate within
an experimental design.

Required columns in both modes:

- file_loc: Directory containing the data file

- file_name: Name of the data file

- group: Experimental group (e.g., "WT", "SMUG1_KO", "male", "female")

- omics: Type of omics data (free-form label, e.g., "ATACseq",
  "CHIPseq", "RNAseq", "5hmu", "CUT&RUN", etc.)

- format: File format ("peaks", "bed", "counts")

- bed_loc: Path to fragment BED file for quantification (required for
  peaks and bed formats). Typically created with bedtools bamtobed.

Additional required columns in experimental mode only:

- bio_rep: Biological replicate identifier

- tech_rep: Technical replicate identifier

- batch: Batch number

The function performs the following validations:

1.  Checks all required columns are present

2.  Removes empty rows

3.  Validates format values against allowed types

4.  Checks data files exist (if check_files = TRUE)

5.  Checks quantification BED files exist for peak/bed formats (if
    check_beds = TRUE)

6.  Generates unique sample_id for each row

Note: The omics column accepts any string value. This allows flexibility
for custom omics types (e.g., "5hmu", "CUT&RUN", "MeDIP") that may be
processed similarly to standard types. The format column determines how
data is loaded.

In experimental mode, sample IDs are generated as:
`<omics>_<group>_<bio_rep><tech_rep>_b<batch>`

If a `timepoint` column is present, it is appended for uniqueness:
`<omics>_<group>_<bio_rep><tech_rep>_b<batch>_t<timepoint>`

## Examples

``` r
if (FALSE) { # \dontrun{
# From CSV file
result <- ELEUTHIA_validate_sample_sheet("data/sample_sheet.csv")

# From data.frame
result <- ELEUTHIA_validate_sample_sheet(my_sample_df)

# Check validation status
if (result$valid) {
  sample_sheet <- result$sample_sheet
} else {
  stop(paste(result$errors, collapse = "\n"))
}

# Skip file checks for quick validation
result <- ELEUTHIA_validate_sample_sheet(df, check_files = FALSE)

} # }
```
