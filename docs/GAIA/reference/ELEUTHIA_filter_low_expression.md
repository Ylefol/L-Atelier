# Filter Low-Expressed Genes

Filters genes with low expression across samples. This is a convenience
wrapper that works with RNA-seq data structure.

## Usage

``` r
ELEUTHIA_filter_low_expression(
  rna_data,
  min_count = 10,
  min_samples = 2,
  verbose = TRUE
)
```

## Arguments

- rna_data:

  A list from ELEUTHIA_load_rnaseq_from_sheet().

- min_count:

  Integer. Minimum count threshold (default = 10).

- min_samples:

  Integer. Minimum samples meeting threshold (default = 2).

- verbose:

  Logical. Print progress (default = TRUE).

## Value

Filtered rna_data list.

## Details

Note: This is a simple count-based filter. For more sophisticated
filtering (e.g., CPM-based), use edgeR::filterByExpr() or similar.
