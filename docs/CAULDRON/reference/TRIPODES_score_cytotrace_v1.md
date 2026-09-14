# Score cells by differentiation state using CytoTRACE v1

Implements the CytoTRACE v1 algorithm (Gulati et al. 2020) to assign a
stemness score to each cell. A score of 1 indicates a highly
undifferentiated (stem-like) state; 0 indicates a fully differentiated
state.

## Usage

``` r
TRIPODES_score_cytotrace_v1(
  sce,
  assay_name = "logcounts",
  n_top_genes = 200L,
  n_neighbours = 10L,
  use_dimred = "PCA",
  n_pcs = 30L,
  recompute_graph = FALSE,
  score_col = "cytotrace_v1",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`. Should have been through the standard
  CAULDRON workflow (normalised, embedded, clustered).

- assay_name:

  Character. Assay to use for scoring. Default `"logcounts"`.

- n_top_genes:

  Integer. Number of top-correlated genes used to compute the raw score.
  Default `200`.

- n_neighbours:

  Integer. Number of nearest neighbours used for graph smoothing. Only
  used when `recompute_graph = TRUE` or when no SNN graph is found.
  Default `10`.

- use_dimred:

  Character. Reduced dimension to use for kNN construction when
  `recompute_graph = TRUE` or no graph is available. Default `"PCA"`.

- n_pcs:

  Integer. Number of PCs to use from `use_dimred`. Default `30`.

- recompute_graph:

  Logical. If `TRUE`, ignore the existing SNN graph and build a fresh
  kNN graph for smoothing. **This is not recommended** — it breaks
  consistency with clustering. A prominent warning is always issued when
  this is set to `TRUE`. Default `FALSE`.

- score_col:

  Character. Name of the `colData` column to store the score in. Default
  `"cytotrace_v1"`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

The input `sce` with:

- `colData(sce)[[score_col]]`: per-cell CytoTRACE score \\\[0, 1\]\\

- `metadata(sce)$cytotrace_v1_run = TRUE`

- `metadata(sce)$cytotrace_v1_col`: name of the score column

## Details

The algorithm exploits the empirical observation that less
differentiated cells express a greater diversity of genes. It correlates
each gene's expression with per-cell gene diversity, selects the most
correlated genes, scores cells by their mean expression of those genes,
then smooths across the neighbourhood graph.

**BPCells note:** If `assay_name` is backed by a BPCells
`IterableMatrix`, it is materialised to a dense matrix at the start of
the function. CytoTRACE v1 requires full matrix operations (row-wise
ranking for Spearman correlation) that cannot be performed on a BPCells
backend.

**Graph reuse:** By default the SNN graph built by
[`TALOS_build_graph`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md)
and stored in `metadata(sce)$snn_graph` is reused for neighbourhood
smoothing. This ensures the CytoTRACE score is smoothed in the same
neighbourhood structure used for clustering. Set
`recompute_graph = TRUE` only if you have a specific reason to use a
different graph — a prominent warning will be issued.

## References

Gulati GS et al. (2020) Single-cell transcriptional diversity is a
hallmark of developmental potential. *Science* 367(6476):405–411.
