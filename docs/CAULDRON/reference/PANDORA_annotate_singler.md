# Reference-based cell type annotation via SingleR

Annotates cells (or clusters) by correlating expression profiles against
a reference dataset using
[`SingleR::SingleR()`](https://rdrr.io/pkg/SingleR/man/SingleR.html). A
pre-loaded `SummarizedExperiment` reference or a shorthand string for a
`celldex` dataset can be supplied. Run
[`PANDORA_list_references`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_list_references.md)
to see available shorthand values.

## Usage

``` r
PANDORA_annotate_singler(
  sce,
  reference,
  labels_col = "label.main",
  cluster_col = NULL,
  assay_name = "logcounts",
  label_col = "singler_label",
  prune = TRUE,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"logcounts"` (or other) assay.

- reference:

  A `SummarizedExperiment` reference, or a character shorthand: one of
  `"HumanPrimaryCellAtlas"`, `"BlueprintEncode"`, `"DICE"`, `"Monaco"`,
  `"Novershtern"` (human), or `"ImmGen"`, `"MouseRNAseq"` (mouse).

- labels_col:

  Character. Column in the reference `colData` to use as labels. Default
  `"label.main"`. Use `"label.fine"` for finer resolution.

- cluster_col:

  Character or `NULL`. If `NULL` (default), annotation is at the cell
  level. If a `colData` column name is provided, expression is
  aggregated per cluster before annotation — one label per cluster is
  broadcast back to all cells in that cluster.

- assay_name:

  Character. Assay to use. Default `"logcounts"`.

- label_col:

  Character. `colData` column for the result. Default `"singler_label"`.

- prune:

  Logical. Replace low-confidence labels with `NA` using SingleR's
  delta-score pruning. Default `TRUE`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `colData(sce)[[label_col]]` populated.
