# Load Multiple Count Datasets from Directories

Convenience wrapper to load multiple count datasets from different
directories. Useful for loading multiple omics types.

## Usage

``` r
ELEUTHIA_load_multi_counts(count_paths, ...)
```

## Arguments

- count_paths:

  Named list of directory paths. Names will be used as dataset
  identifiers.

- ...:

  Additional arguments passed to ELEUTHIA_load_counts_from_dir

## Value

Named list of data frames, one per directory

## Examples

``` r
if (FALSE) { # \dontrun{
# Load RNA-seq and ATAC-seq data
datasets <- ELEUTHIA_load_multi_counts(
  count_paths = list(
    rnaseq = "path/to/rnaseq/",
    atacseq = "path/to/atacseq/"
  )
)

} # }
```
