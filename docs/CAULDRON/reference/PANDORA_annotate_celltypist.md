# Cell type annotation via CellTypist

Annotates cells using a pre-trained CellTypist model. The model is
downloaded on first use to the celltypist cache (`~/.celltypist/`). All
Python execution is handled by basilisk in an isolated environment — no
manual Python setup is required.

## Usage

``` r
PANDORA_annotate_celltypist(
  sce,
  model,
  majority_voting = TRUE,
  cluster_col = NULL,
  force_update = FALSE,
  label_col = "celltypist_label",
  conf_col = "celltypist_conf",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"counts"` assay.

- model:

  Character. Model filename, e.g. `"Mouse_Whole_Brain.pkl"`. Use
  [`PANDORA_list_celltypist_models`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_list_celltypist_models.md)
  to browse available models.

- majority_voting:

  Logical. Use majority voting within over-clusters for more coherent
  labels. Default `TRUE`.

- cluster_col:

  Character or `NULL`. `colData` column with an existing per-cell
  clustering to use as CellTypist's `over_clustering` (passed through
  unchanged, aligned to cell order), instead of letting CellTypist build
  its own neighbour graph and over-cluster from scratch. Ignored when
  `majority_voting = FALSE`. Default `NULL`.

- force_update:

  Logical. Re-download the model even if already cached. Default
  `FALSE`.

- label_col:

  Character. `colData` column for the result. Default
  `"celltypist_label"`.

- conf_col:

  Character or `NULL`. If non-`NULL` and `majority_voting = FALSE`,
  stores the per-cell confidence score in this `colData` column. Default
  `"celltypist_conf"`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `colData(sce)[[label_col]]` populated.

## Details

Raw counts are passed to Python and normalised to 10,000 counts per cell
followed by log1p transformation (scanpy's `normalize_total` + `log1p`)
inside the basilisk environment. This is the exact normalisation
CellTypist requires, regardless of which R normalisation was applied to
`logcounts`.

With `majority_voting = TRUE` (default), CellTypist over-clusters the
data internally and assigns a consensus label per cluster by majority
vote, which typically improves coherence. With `FALSE`, each cell
receives an independent probabilistic prediction.

By default, CellTypist's internal over-clustering builds its own PCA +
neighbour graph from scratch on the (freshly normalised) counts passed
in — independent of any neighbour graph or clustering already computed
upstream in the SCE (e.g. after batch correction). Supplying
`cluster_col` passes an existing per-cell clustering straight to
CellTypist's `over_clustering` argument as an array aligned to cell
order, which skips its internal neighbour-graph construction entirely
(verified against the installed celltypist's
[`annotate()`](https://ggplot2.tidyverse.org/reference/annotate.html)/`Classifier.majority_vote()`
source: an array-like `over_clustering` bypasses `over_cluster()` and
goes directly to majority voting) — so the consensus vote reflects your
own established clusters/embedding rather than a second, independently
constructed graph.
