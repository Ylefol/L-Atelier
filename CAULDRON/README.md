# CAULDRON

Single-cell and single-nuclei RNA-seq analysis toolkit. The single-cell counterpart to GAIA, designed to handle the full scRNA-seq / snRNA-seq workflow: from raw count matrices through QC, normalization, clustering, annotation, differential expression, trajectory inference, and cell-cell communication analysis.

Part of the [L-Atelier](../) repository. Built on `SingleCellExperiment` (Bioconductor). Follows the same modular architecture as GAIA, with subordinate modules named after Hephaestus's creations from Greek mythology.

> **Status:** Early development. Module files are scaffolded but most functions are not yet implemented.

---

## Modules

| Module   | Named after | Responsibility |
|----------|-------------|----------------|
| TALARIA  | Winged sandals of Hermes — pure transport | I/O: load 10x MEX, h5, h5ad, loom; export; format conversion |
| AEGIS    | Divine shield of Athena/Zeus — protects | QC: per-cell metrics, doublet detection, ambient RNA removal |
| PYRI     | Sacred forge fire — purifies raw ore | Preprocessing: normalization, HVG selection, scaling, batch correction |
| TALOS    | Bronze automaton that traversed Crete | Embedding: PCA, UMAP, tSNE, graph construction, clustering |
| PANDORA  | Fashioned by Hephaestus, gifted by all gods | Annotation: marker-based and reference-based cell type labelling |
| KERAUNOS | Thunderbolts forged by Hephaestus — decisive | Expression: cluster markers, condition DE, pseudobulk DE, GSEA |
| TRIPODES | Self-moving golden tripods — self-directed motion | Trajectory: pseudotime, RNA velocity, lineage inference |
| DIKTYON  | Invisible net/snare — captures connections | Communication: cell-cell signalling (CellChat, NicheNet, LIANA) |
| ASPIS    | Shield of Achilles — depicts everything | Visualization: UMAP, violin, dot, feature, heatmap, trajectory plots |
| KHALKOS  | Bronze, the foundational material of the forge | Utilities: shared helpers, object manipulation, configuration |

---

## Primary Object

CAULDRON is built on `SingleCellExperiment` (Bioconductor), chosen for its programmatic flexibility, native compatibility with the Bioconductor single-cell ecosystem (scran, scDblFinder, SingleR, Monocle3, edgeR/DESeq2), and structural near-equivalence with AnnData — making Scanpy/Python experience directly transferable. Seurat interoperability is handled via TALARIA conversion functions without a full Seurat dependency.

---

## Installation

```r
install.packages("devtools")
devtools::install("path/to/L-Atelier/CAULDRON")
library(CAULDRON)
```

---

## Development

```r
devtools::document("CAULDRON")
devtools::install("CAULDRON")
```

Implementation status will be tracked in `STATUS/` at the repository root once functions are implemented. All function names follow the same convention as GAIA: `MODULE_function_name` (e.g., `AEGIS_filter_cells()`, `TALOS_run_umap()`).
