# Load Count Data from Directory

Loads count files from a directory and merges them into a single data
frame. Expects files with gene_id in first column and counts in second
column. Handles file naming conventions with numeric prefixes.

## Usage

``` r
ELEUTHIA_load_counts_from_dir(
  count_path,
  file_pattern = NULL,
  gene_col = 1,
  count_col = 2,
  clean_sample_names = TRUE
)
```

## Arguments

- count_path:

  Character string. Path to directory containing count files. Each file
  should be a tab-delimited text file with gene_id and count columns.

- file_pattern:

  Optional regex pattern to filter files (default = NULL, uses all
  files)

- gene_col:

  Name or index of gene ID column (default = 1)

- count_col:

  Name or index of count column (default = 2)

- clean_sample_names:

  Logical, should sample names be cleaned? (default = TRUE) Removes
  numeric prefixes like "123_456-" from filenames.

## Value

Data frame with genes as rows and samples as columns. Row names are gene
IDs.

## Details

This function:

- Reads all files in the specified directory (or files matching
  file_pattern)

- Extracts sample names from filenames (removes file extensions)

- Optionally cleans sample names by removing numeric prefixes

- Merges all count files by gene_id

- Returns a matrix with genes (rows) × samples (columns)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load count data from directory
counts <- ELEUTHIA_load_counts_from_dir("path/to/counts/")

# With file pattern filtering
counts <- ELEUTHIA_load_counts_from_dir(
  "path/to/counts/",
  file_pattern = "\\.txt$"
)

} # }
```
