# Export a first-steps analysis

Saves all outputs of a CAULDRON first-steps pipeline run into a
structured directory tree. The function auto-detects what has been
computed from the SCE object and only exports what exists, so it works
whether the user ran the full pipeline or just loaded data and produced
a single UMAP. Sweep objects (if any were run) are passed explicitly.

## Usage

``` r
TALARIA_export_first_steps(
  sce,
  output_dir,
  hvg_sweep = NULL,
  pc_sweep = NULL,
  k_sweep = NULL,
  resolution_sweep = NULL,
  umap_sweep = NULL,
  tsne_sweep = NULL,
  colour_by = "cluster",
  sample_col = NULL,
  assay_name = "logcounts",
  format = c("png", "pdf"),
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object — the main analysis object. The
  function inspects its state (assays, reducedDims, colData, metadata)
  to determine what has been computed and should be exported.

- output_dir:

  Character. Path to the root output directory. Created (recursively) if
  it does not exist.

- hvg_sweep:

  A `pyri_hvg_sweep` object from
  [`PYRI_tune_hvg`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_tune_hvg.md),
  or `NULL` (default).

- pc_sweep:

  A `talos_pc_sweep` object from
  [`TALOS_tune_pc`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_pc.md),
  or `NULL` (default).

- k_sweep:

  A `talos_k_sweep` object from
  [`TALOS_tune_k`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_k.md),
  or `NULL` (default).

- resolution_sweep:

  A `talos_resolution_sweep` object from
  [`TALOS_tune_resolution`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_resolution.md),
  or `NULL` (default).

- umap_sweep:

  A `talos_embedding_sweep` object (type `"umap"`) from
  [`TALOS_tune_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_umap.md),
  or `NULL` (default).

- tsne_sweep:

  A `talos_embedding_sweep` object (type `"tsne"`) from
  [`TALOS_tune_tsne`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_tsne.md),
  or `NULL` (default).

- colour_by:

  Character. `colData` column used to colour embedding plots. Default
  `"cluster"`. If not found, falls back to the first factor or character
  column in `colData` with a console note.

- sample_col:

  Character or `NULL`. If provided, an additional QC plot grouped by
  this `colData` column is saved alongside the ungrouped version.
  Default `NULL`.

- assay_name:

  Character. Assay used when `colour_by` is a gene name in the embedding
  plots. Default `"logcounts"`.

- format:

  Character. Output format for plots: `"png"` (default, 600 DPI raster)
  or `"pdf"` (vector, no DPI).

- verbose:

  Logical. Print a step-by-step progress log and a final file manifest
  summary. Default `TRUE`.

## Value

Invisibly returns a list with `output_dir` (character) and `files`
(character vector of all paths successfully written).

## Details

**Output directory structure:**

    output_dir/
    ├── sce/
    │   └── sce.rds
    ├── plots/
    │   ├── qc_metrics.png
    │   ├── qc_metrics_by_sample.png  (if sample_col provided)
    │   ├── pca_elbow.png
    │   ├── sweep_hvg.png
    │   ├── sweep_pc.png
    │   ├── sweep_k.png
    │   ├── sweep_resolution.png
    │   ├── sweep_umap.png
    │   ├── sweep_tsne.png
    │   ├── grid_umap.png
    │   ├── grid_tsne.png
    │   ├── umap.png
    │   └── tsne.png
    └── tables/
        ├── analysis_parameters.csv
        ├── cell_metadata.csv
        ├── sweep_hvg_results.csv
        ├── sweep_pc_results.csv
        ├── sweep_k_results.csv
        ├── sweep_resolution_results.csv
        ├── sweep_umap_results.csv
        └── sweep_tsne_results.csv

**BPCells backend:** `saveRDS(sce, ...)` is used for both in-memory and
BPCells-backed objects — the BPCells matrix is already on disk and the
RDS stores only the SCE shell with a path reference. The BPCells
directory (stored in `metadata(sce)$bpcells_dir`) must remain accessible
alongside the RDS for the object to be reloadable.

**Analysis parameters table:** a two-column CSV (`parameter`, `value`)
summarising the key choices made during the analysis — `n_hvgs`,
`n_pcs_used`, `k`, `clustering_method`, `resolution`, `n_neighbors`,
`min_dist`, and `perplexity`. Values are populated from sweep objects
where provided; parameters whose sweep was not run are recorded as `NA`.
