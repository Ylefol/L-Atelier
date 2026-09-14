# CAULDRON

Single-cell and single-nuclei RNA-seq analysis toolkit. The single-cell
counterpart to GAIA, covering single cell omics. Currently this only
covers scRNAseq - it will be expanded as the need arises.

Part of the [L-Atelier](https://ylefol.github.io/L-Atelier/) repository.
Built on `SingleCellExperiment` (Bioconductor). Follows the same modular
architecture as GAIA, with subordinate modules named after Hephaestus’s
creations from Greek mythology.

> **Status:** This package, like the others, is under active
> development. Some elements are unlikely to change while others are
> very prone to change.

------------------------------------------------------------------------

## Modules

| Module   | Named after                                       | Responsibility                                                                                                                          |
|----------|---------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| TALARIA  | Winged sandals of Hermes — pure transport         | I/O: load 10x MEX/H5/H5AD/Loom, Seurat interop, BPCells backing, export (matrices, DE/GSEA/ORA results, interactive R Markdown widgets) |
| AEGIS    | Divine shield of Athena/Zeus — protects           | QC: per-cell metrics, doublet detection (scDblFinder), cell filtering                                                                   |
| PYRI     | Sacred forge fire — purifies raw ore              | Preprocessing: normalization, HVG selection and tuning                                                                                  |
| TALOS    | Bronze automaton that traversed Crete             | Embedding & clustering: PCA, Harmony, UMAP, tSNE, scVI, SNN graph construction/clustering, parameter sweeps                             |
| PANDORA  | Fashioned by Hephaestus, gifted by all gods       | Annotation: SingleR, scType, CellTypist, manual label assignment, ortholog conversion                                                   |
| KERAUNOS | Thunderbolts forged by Hephaestus — decisive      | Expression: cluster markers, cell-level and pseudobulk (DESeq2) DE, GSEA/ORA (single and pairwise), propeller proportion testing        |
| TRIPODES | Self-moving golden tripods — self-directed motion | Trajectory: Monocle3 pseudotime, RNA velocity (velociraptor/scVelo), CytoTRACE stemness scoring                                         |
| DIKTYON  | Invisible net/snare — captures connections        | Communication: cell-cell signalling (CellChat, NicheNet, LIANA) — **planned, not yet implemented**                                      |
| ASPIS    | Shield of Achilles — depicts everything           | Visualization: atlas plots, embeddings, QC, markers, composition, velocity, pseudotime, trajectory                                      |
| KHALKOS  | Bronze, the foundational material of the forge    | Utilities: SCE validation/merging, assay accessors, colour palettes, ggplot2 theme, logging                                             |

------------------------------------------------------------------------

## Primary Object

CAULDRON is built on `SingleCellExperiment` (Bioconductor), chosen for
its programmatic flexibility, native compatibility with the Bioconductor
single-cell ecosystem (scran, scDblFinder, SingleR, Monocle3,
edgeR/DESeq2), and structural near-equivalence with AnnData — making
Scanpy/Python experience directly transferable. Seurat interoperability
is handled via TALARIA conversion functions
([`TALARIA_to_seurat()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_to_seurat.md)
/
[`TALARIA_from_seurat()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_from_seurat.md))
without a full Seurat dependency.

------------------------------------------------------------------------

## Installation

``` r
install.packages("devtools")
devtools::install("path/to/L-Atelier/CAULDRON")
library(CAULDRON)
```

I haven’t set-up the required package list for this installation yet -
if you choose to use it you will have to install packages as they are
called by CAULDRON.

------------------------------------------------------------------------

## Typical Workflow (sc/sn RNAseq)

    TALARIA_load_10x() / TALARIA_load_h5ad() / TALARIA_load_h5() / TALARIA_load_loom()
      ↓
    AEGIS_compute_qc_metrics() → AEGIS_detect_doublets() → AEGIS_filter_cells()
      ↓
    PYRI_normalize() → PYRI_select_hvg()
      ↓
    TALOS_run_pca() → TALOS_run_harmony()  # if batch correction is needed
      ↓
    TALOS_run_umap() / TALOS_run_tsne()  +  TALOS_build_graph() → TALOS_cluster()
      ↓
    PANDORA_annotate_singler() / PANDORA_annotate_celltypist() / PANDORA_assign_labels()
      ↓
    KERAUNOS_find_markers() / KERAUNOS_de_pseudobulk() → KERAUNOS_gsea() / KERAUNOS_ora()
      ↓
    ASPIS_plot_umap() / ASPIS_plot_marker_dotplot() / ... → TALARIA_export_de() / TALARIA_export_first_steps()

Trajectory work
([`TRIPODES_run_monocle()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_monocle.md),
[`TRIPODES_run_velocity()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_velocity.md),
[`TRIPODES_score_cytotrace_v1()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_score_cytotrace_v1.md))
and propeller-based composition testing
([`KERAUNOS_propeller_proportions()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_propeller_proportions.md))
branch off after annotation, depending on the biological question.

------------------------------------------------------------------------

## References

Methods/tools wrapped or reimplemented within CAULDRON:

**QC & preprocessing** - **scDblFinder** (doublet detection) — Doublet
identification in single-cell sequencing data using scDblFinder. doi:
[10.12688/f1000research.73600.2](https://doi.org/10.12688/f1000research.73600.2) -
**scran / scater** (pooling-based normalization, HVG selection) — A
step-by-step workflow for low-level analysis of single-cell RNA-seq data
with Bioconductor. doi:
[10.12688/f1000research.9501.2](https://doi.org/10.12688/f1000research.9501.2) -
**scran pooling normalization** (`computeSumFactors`, specifically) —
Pooling across cells to normalize single-cell RNA sequencing data with
many zero counts. doi:
[10.1186/s13059-016-0947-7](https://doi.org/10.1186/s13059-016-0947-7)

**Embedding & clustering** - **UMAP** (uwot) — UMAP: Uniform Manifold
Approximation and Projection for Dimension Reduction. arXiv:
[1802.03426](https://arxiv.org/abs/1802.03426) - **t-SNE** (Rtsne) —
Visualizing Data using t-SNE. J Mach Learn Res. 9:2579-2605, 2008. -
**Harmony** (batch integration) — Fast, sensitive and accurate
integration of single-cell data with Harmony. doi:
[10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0) -
**scVI** (variational latent embedding, via a custom basilisk-managed
script wrapping `scvi-tools` directly) — Deep generative modeling for
single-cell transcriptomics. doi:
[10.1038/s41592-018-0229-2](https://doi.org/10.1038/s41592-018-0229-2) -
**Louvain clustering** (igraph) — Fast unfolding of communities in large
networks. doi:
[10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008) -
**Leiden clustering** (igraph) — From Louvain to Leiden: guaranteeing
well-connected communities. doi:
[10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)

**Annotation** - **SingleR** — Reference-based analysis of lung
single-cell sequencing reveals a transitional profibrotic macrophage.
doi:
[10.1038/s41590-018-0276-y](https://doi.org/10.1038/s41590-018-0276-y) -
**scType** *(marker-scoring algorithm reimplemented in R against a
bundled ScTypeDB, not the original authors’ package/script)* —
Fully-automated and ultra-fast cell-type identification using specific
marker combinations from single-cell transcriptomic data. doi:
[10.1038/s41467-022-28803-w](https://doi.org/10.1038/s41467-022-28803-w) -
**CellTypist** — Cross-tissue immune cell analysis reveals
tissue-specific features in humans. doi:
[10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197) -
**orthogene** (ortholog conversion) — orthogene: a package for
streamlining the interspecies gene mapping and conversion. doi:
[10.1093/bioinformatics/btac550](https://doi.org/10.1093/bioinformatics/btac550)

**Expression & differential analysis** - **DESeq2** (pseudobulk DE) —
Moderated estimation of fold change and dispersion for RNA-seq data with
DESeq2. doi:
[10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) -
**fgsea** (GSEA) — Fast gene set enrichment analysis. bioRxiv. doi:
[10.1101/060012](https://doi.org/10.1101/060012) - **gprofiler2** (ORA)
— g:Profiler: hierarchical, functional annotation and analysis of gene
lists. doi:
[10.12688/f1000research.24956.2](https://doi.org/10.12688/f1000research.24956.2) -
**MSigDB** (via msigdbr) — Gene set enrichment analysis: a
knowledge-based approach for interpreting genome-wide expression
profiles. doi:
[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102) -
**propeller** (speckle package, cell-type proportion testing) —
propeller: testing for differences in cell type proportions in single
cell data. doi:
[10.1093/bioinformatics/btac582](https://doi.org/10.1093/bioinformatics/btac582)

**Trajectory** - **Monocle3** — The single-cell transcriptional
landscape of mammalian organogenesis. doi:
[10.1038/s41586-019-0969-x](https://doi.org/10.1038/s41586-019-0969-x) -
**CytoTRACE v1** *(reimplemented natively in R by
[`TRIPODES_score_cytotrace_v1()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_score_cytotrace_v1.md),
not calling the original CytoTRACE package/CLI)* — Single-cell
transcriptional diversity is a hallmark of developmental potential. doi:
[10.1126/science.aax0249](https://doi.org/10.1126/science.aax0249) -
**scVelo** (RNA velocity, via a custom basilisk-managed script — see
note below) — Generalizing RNA velocity to transient cell states through
dynamical modeling. doi:
[10.1038/s41587-020-0591-5](https://doi.org/10.1038/s41587-020-0591-5)

**I/O interoperability** - **Seurat** (conversion only, via
[`TALARIA_to_seurat()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_to_seurat.md)/[`TALARIA_from_seurat()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_from_seurat.md))
— Integrated analysis of multimodal single-cell data. doi:
[10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048) -
**DropletUtils** (`read10xCounts()`/`write10xCounts()` only — CAULDRON
does not use this package’s `emptyDrops()` cell-calling method) —
EmptyDrops: distinguishing cells from empty droplets in droplet-based
single-cell RNA sequencing data. doi:
[10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y) -
**BPCells** (on-disk matrix backend) — GitHub-only tool
([bnprks/BPCells](https://github.com/bnprks/BPCells)). Treat as an
engineering dependency, not a citable method.

------------------------------------------------------------------------

## Development

All function names follow the same convention as GAIA:
`MODULE_function_name` (e.g.,
[`AEGIS_filter_cells()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_filter_cells.md),
[`TALOS_run_umap()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_umap.md)).
Per-module implementation status — the full, current function inventory,
in-progress work, and TODOs — is tracked in
`inst/status/<MODULE>_STATUS.txt` (e.g.,
`inst/status/KERAUNOS_STATUS.txt`), following the same format as GAIA’s
`STATUS/` files at the repository root.
