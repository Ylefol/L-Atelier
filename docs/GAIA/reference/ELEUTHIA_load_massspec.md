# Load DIA Mass Spectrometry Data (DIA-NN report.pg_matrix)

Reads a DIA-NN `report.pg_matrix` file (.xlsx or tab-delimited text),
parses sample identifiers from the column names (which are Windows or
POSIX file paths in DIA-NN output), log2-transforms intensities, and
returns a structured `massspec_data` object.

This function performs no QC filtering. To remove controls, pool
samples, and high-missingness samples/proteins, pass the result to
[`HADES_filter_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_massspec.md)
and
[`HADES_filter_massspec_proteins()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_massspec_proteins.md).

## Usage

``` r
ELEUTHIA_load_massspec(
  ms_file,
  sheet = 1,
  metadata_file = NULL,
  metadata_sheet = 1,
  sample_col = "SampleID",
  n_annotation_cols = NULL,
  verbose = TRUE
)
```

## Arguments

- ms_file:

  Character. Path to the DIA-NN `report.pg_matrix` file. Accepted
  formats: `.xlsx`, `.tsv`, `.txt`, `.csv`.

- sheet:

  Integer or character. Sheet to read if `ms_file` is an xlsx workbook.
  Default: 1.

- metadata_file:

  Character or NULL. Path to a sample metadata file (.csv, .tsv, or
  .xlsx). Joined onto `$sample_meta` by `sample_col`. Default: NULL.

- metadata_sheet:

  Integer or character. Sheet for xlsx metadata files. Default: 1.

- sample_col:

  Character. Column in `metadata_file` used as the join key. Parsed
  sample IDs (numeric strings, SEP###-# labels, etc.) are matched
  against this column. Default: `"SampleID"`.

- n_annotation_cols:

  Integer or NULL. Number of leading columns that are protein
  annotations rather than sample intensities. `NULL` (default) uses
  auto-detection: annotation columns are those whose names do not
  contain a file-path separator (`/` or `\`) or a `.d` extension.

- verbose:

  Logical. Print loading progress. Default: TRUE.

## Value

An S3 object of class `"massspec_data"` containing:

- wide:

  Numeric matrix: proteins (rows) x samples (cols). Rownames =
  Protein.Group. Colnames = parsed SampleID. Values are log2-transformed
  intensities; zeros and negatives become NA.

- sample_meta:

  Per-sample data.frame: SampleID, SampleType, PlateID, plus any columns
  joined from `metadata_file`.

- protein_meta:

  Per-protein data.frame with annotation columns from the input file
  (e.g. Protein.Group, Protein.Names, Genes, First.Protein.Description).
  Rownames match `rownames($wide)`.

- params:

  List of loader parameters used (for provenance).

## Details

**Expected input format:** Column names must follow the standardised
convention produced by the project's `prep_MS_data.R` script:
`{SampleID}_Plate {N}` (e.g. `"1027_Plate 1"`, `"SEP001-3_Plate 7"`,
`"pool1_Plate 12"`). Raw DIA-NN file paths are not supported; data must
be prepared through `prep_MS_data.R` first.

**Sample ID parsing:** `PlateID` is extracted from the `_Plate N`
suffix; `SampleID` is everything before it. Column names that do not
match the `_Plate N` suffix pattern yield `PlateID = NA` and
`SampleID = original column name`.

**Sample type classification:** identity is not guessed from ID shape –
a column's SampleID can be any format. Instead:

- `"POOL"` — pooled QC injections (`pool*`). Present on every plate for
  batch drift monitoring. No metadata entry expected.

- `"SAMPLE"` — any non-pool column whose parsed SampleID is found in
  `metadata_file` (by `sample_col`). Its biological distinction (group,
  timepoint, cohort) lives in the joined metadata columns, not in
  SampleType.

- `"OTHER"` — any non-pool column that does *not* match a row in
  `metadata_file`. Typically wells belonging to a different sub-study
  run on the same plates/DIA-NN batch (e.g. a sample-stability QC arm
  with its own SubjectID/Group scheme tracked only in the plate layout,
  not in the main clinical metadata) rather than malformed IDs. These
  columns cannot be joined to metadata and are dropped by
  [`HADES_filter_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_massspec.md)
  by default. If `metadata_file` is `NULL`, there is nothing to check a
  match against, so every non-pool column is classified `"SAMPLE"`.

**Log2 transformation:** Raw DIA intensities are on a linear scale.
Intensities \\\le 0\\ (zeros = not detected; negatives = artefact) are
set to NA before [`log2()`](https://rdrr.io/r/base/Log.html) is applied.
Do NOT transform again before passing to
[`ARTEMIS_limma_de()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_limma_de.md).

**Normalization:** Not performed at load time. Pass `$wide` to
[`POSEIDON_normalize_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_normalize_massspec.md)
after filtering.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load only
ms_raw <- ELEUTHIA_load_massspec("report.pg_matrix.xlsx")

# Load with Olink-matched metadata, then filter
ms_raw <- ELEUTHIA_load_massspec(
  ms_file       = "report.pg_matrix.xlsx",
  metadata_file = "olink_sample_meta.csv",
  sample_col    = "SampleID"
)
ms <- HADES_filter_massspec(ms_raw, keep_types = "SAMPLE")
ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
ms$wide <- POSEIDON_normalize_massspec(ms$wide)
} # }
```
