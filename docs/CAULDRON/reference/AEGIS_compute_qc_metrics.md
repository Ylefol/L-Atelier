# Compute per-cell QC metrics

Adds standard per-cell quality metrics to the `colData` of a
`SingleCellExperiment`:

## Usage

``` r
AEGIS_compute_qc_metrics(
  sce,
  assay_name = NULL,
  mito_prefix = "mt-",
  ribo_prefixes = c("Rps", "Rpl"),
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- assay_name:

  Character. Assay to compute metrics from. If `NULL` (default), tries
  `"counts"` then `"X"` then the first assay.

- mito_prefix:

  Character. Prefix identifying mitochondrial genes. Default `"mt-"`
  (mouse). Use `"MT-"` for human.

- ribo_prefixes:

  Character vector. Prefixes identifying ribosomal genes. Default
  `c("Rps", "Rpl")` (mouse). Use `c("RPS", "RPL")` for human. Pass
  `NULL` to skip ribosomal metrics.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

The input `SingleCellExperiment` with QC columns added to `colData`.

## Details

- `sum` — total UMI / read counts per cell

- `detected` — number of genes with non-zero counts

- `subsets_mt_percent` — \\

- `subsets_ribo_percent` — \\ (if ribosomal genes are found)

- `log10_genes_per_umi` — complexity score:
  `log10(detected) / log10(sum)`; values near 0.8 indicate
  high-complexity libraries
