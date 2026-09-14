# Eleuthia - RNA-seq Loading Functions

Functions for loading RNA-seq count data. Load a Single Count File

Reads a simple two-column count file (gene_id, count).

## Usage

``` r
ELEUTHIA_load_count_file(file_path, gene_col = 1, count_col = 2, header = NA)
```

## Arguments

- file_path:

  Character string. Path to the count file.

- gene_col:

  Integer. Column index for gene IDs (default = 1).

- count_col:

  Integer. Column index for counts (default = 2).

- header:

  Logical or NA. Does the file have a header row? NA (default)
  auto-detects by checking whether the first line's `count_col` field
  parses as a number – if it doesn't, the first line is treated as a
  header (e.g. HORIZON's `<sample_id>_counts.txt` files, which are
  always written with a `gene_id`/`count` header). Set explicitly to
  TRUE/FALSE to skip detection.

## Value

A named numeric vector of counts, with gene IDs as names.
