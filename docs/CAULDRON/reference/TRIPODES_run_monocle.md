# Pseudotime and trajectory inference via Monocle3

Runs the Monocle3 trajectory inference pipeline on a
`SingleCellExperiment`, re-using an existing 2D embedding from TALOS
rather than re-computing it. Pseudotime and branch assignments are
written back into `colData`. The principal graph is stored in `metadata`
for downstream plotting.

## Usage

``` r
TRIPODES_run_monocle(
  sce,
  root_cells = NULL,
  root_cluster = NULL,
  root_by_cytotrace = FALSE,
  cytotrace_top_pct = 5,
  cytotrace_col = "cytotrace_v1",
  cluster_col = "cluster",
  use_partition = TRUE,
  close_loop = FALSE,
  pseudotime_col = "monocle_pseudotime",
  branch_col = "monocle_branch",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` that has been through the standard CAULDRON
  workflow (normalised, embedded, clustered via TALOS).

- root_cells:

  Character vector of cell barcodes to use as root. Exactly one of
  `root_cells`, `root_cluster`, or `root_by_cytotrace = TRUE` must be
  provided.

- root_cluster:

  Character. A single cluster label from `colData(sce)[[cluster_col]]`.
  All cells in that cluster are used as the root.

- root_by_cytotrace:

  Logical. Use the top `cytotrace_top_pct`\\ cells by CytoTRACE v1 score
  as root. Requires
  [`TRIPODES_score_cytotrace_v1`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_score_cytotrace_v1.md)
  to have been run. Default `FALSE`.

- cytotrace_top_pct:

  Numeric. Percentage of top-scoring cells to use as root when
  `root_by_cytotrace = TRUE`. Default `5`.

- cytotrace_col:

  Character. `colData` column containing CytoTRACE scores. Default
  `"cytotrace_v1"`.

- cluster_col:

  Character. `colData` column containing TALOS cluster labels. Used when
  `root_cluster` is specified. Default `"cluster"`.

- use_partition:

  Logical. If `TRUE`, Monocle3 learns separate graphs per partition
  (recommended for data with disconnected trajectories). If `FALSE`, a
  single graph is learned across all cells. Default `TRUE`.

- close_loop:

  Logical. Whether Monocle3 should attempt to close loops in the
  principal graph. Set `TRUE` only if your biology has cyclic dynamics.
  Default `FALSE`.

- pseudotime_col:

  Character. Name of the `colData` column for pseudotime values. Default
  `"monocle_pseudotime"`.

- branch_col:

  Character. Name of the `colData` column for Monocle3 partition
  assignments. Default `"monocle_branch"`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

The input `sce` with:

- `colData(sce)[[pseudotime_col]]`: per-cell pseudotime; `NA` for cells
  in partitions with no root.

- `colData(sce)[[branch_col]]`: Monocle3 partition assignment.

- `metadata(sce)$monocle_run = TRUE`

- `metadata(sce)$monocle_graph`: principal graph as an `igraph` object
  (topology only).

- `metadata(sce)$monocle_graph_nodes`: matrix of principal graph node
  coordinates in UMAP space (nodes \\\times\\ 2).

- `metadata(sce)$monocle_dimred`: always `"UMAP"`.

## Details

**Embedding reuse:** The pre-computed `"UMAP"` embedding from TALOS is
injected directly into the Monocle3 `cell_data_set`, bypassing
Monocle3's internal normalisation and PCA steps. This ensures the
trajectory is learned on exactly the same layout used for clustering.

**UMAP only:** Monocle3's `learn_graph()` and `cluster_cells()`
functions are hardcoded to require a `"UMAP"` reducedDim regardless of
any user-supplied `reduction_method` argument. This is a known, unfixed
limitation of the monocle3 package. `TRIPODES_run_monocle()` therefore
requires a UMAP to be present (run
[`TALOS_run_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_umap.md)
first) and will error immediately if it is absent rather than silently
failing inside Monocle3.

**Graph partitions:** Monocle3 clusters the cells internally on the
injected embedding (via `cluster_cells()`) to determine graph
partitions. These internal partitions are used only by `learn_graph()`
and do not replace or override the TALOS cluster labels. Root selection
by `root_cluster` refers to your TALOS cluster labels, not Monocle3's
internal partitions.

**Inf pseudotime:** Cells in partitions that do not contain a root cell
will receive `Inf` pseudotime from Monocle3. These are stored as `NA` in
`colData`.
