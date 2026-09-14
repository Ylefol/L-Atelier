# GAIA

Multi-omics integration toolkit for bulk level omics. Covers quality
control, preprocessing, annotation, enrichment, integration, machine
learning, statistical analysis, visualization, reporting, and more.

Part of the [L-Atelier](https://ylefol.github.io/L-Atelier/) repository.
Its nine subordinate modules are each named after a figure from Greek
mythology.

> **Status:** This package, like the others, is under active
> development. Some elements are unlikely to change while others are
> very prone to change.

------------------------------------------------------------------------

## Installation

See the
[Installation](https://ylefol.github.io/L-Atelier/GAIA/articles/installation.md)
guide for setup instructions and optional dependencies.

------------------------------------------------------------------------

## Modules

GAIA orchestrates nine subordinate modules. Individual modules can also
be used independently.

| Module     | Named after                                                          | Responsibility                                                                                                                                                                                                |
|------------|----------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Hades      | Deity of the underworld — removes the “dead”                         | Quality control: ATAC-seq QC, mass spec QC, Olink QC, mixed-data QC/filtering                                                                                                                                 |
| Poseidon   | God of the seas — controls and purifies waters                       | Preprocessing & normalization: batch correction (ComBat/limma), feature filtering, imputation, normalization, regression prep, sample management                                                              |
| Apollo     | God of knowledge and prophecy — brings wisdom through annotation     | Annotation & enrichment: peak annotation, gene enrichment (GO/KEGG/GSEA/ORA), motif enrichment, isoform annotation, peak-to-gene linking, activity networks (decoupleR/OmnipathR)                             |
| Hephaestus | Master craftsman — forges disparate materials into unified artifacts | Integration: sparse CCA (ConvCCA, RelPMDCCA), MOFA2 multi-omics factor analysis                                                                                                                               |
| Minerva    | Roman goddess of wisdom and strategic warfare                        | Machine learning & validation: sCCA cross-validation, penalty functions (lasso/group lasso)                                                                                                                   |
| Artemis    | Goddess of the hunt — precisely targets statistical insights         | Analysis & statistics: differential expression (DESeq2/limma), clustering (PAM/consensus/PART/WGCNA), discriminative analysis, FAMD, time series, CIBERSORT deconvolution, isoform switching, ATAC annotation |
| Aether     | Pure upper air breathed by gods — clarifies and illuminates results  | Visualization: coverage tracks, circos, PCA/clustering, enrichment, differential, time series, WGCNA, MOFA, motif and network plots                                                                           |
| Eleuthia   | Cradle facility that brought new life into the world                 | Import, export & reporting: sample sheets, data import/quantification, DEA/ATAC/time series/composition report export, BigWig, HOMER export, logging                                                          |
| Demeter    | Goddess of agriculture and harvest — cultivates the foundation       | Utilities: annotation/correlation/matrix/grid helpers, result loaders, simulated data generation                                                                                                              |

------------------------------------------------------------------------

## Development

All source files are flat under `GAIA/R/` with module-prefixed names
(e.g., `artemis_wgcna.R`, `aether_circos.R`). Function names follow the
same convention: `MODULE_function_name` (e.g.,
`ARTEMIS_wgcna_detect_modules`, `AETHER_plot_circos`).

Module implementation status is tracked in `STATUS/` at the repository
root (e.g., `STATUS/Artemis_STATUS.txt`).

------------------------------------------------------------------------

## References

Methods implemented or wrapped within GAIA:

**Integration** - **sCCA** *(methodology GAIA’s own sparse CCA
implementation follows)* — Integrating multi-OMICS data through sparse
canonical correlation analysis for the prediction of complex traits: a
comparison study. doi:
[10.1093/bioinformatics/btaa530](https://doi.org/10.1093/bioinformatics/btaa530) -
**PMA** (sparse CCA, via `minerva_sCCA_CV.R`) — A penalized matrix
decomposition, with applications to sparse principal components and
canonical correlation analysis. doi:
[10.1093/biostatistics/kxp008](https://doi.org/10.1093/biostatistics/kxp008) -
**MOFA2/MOFA+** — MOFA+: a statistical framework for comprehensive
integration of multi-modal single-cell data. doi:
[10.1186/s13059-020-02015-1](https://doi.org/10.1186/s13059-020-02015-1)

**Batch correction** - **ComBat** (sva) — Adjusting batch effects in
microarray expression data using empirical Bayes methods. doi:
[10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037) -
**limma** (`removeBatchEffect()`, and the DE workflow below) — limma
powers differential expression analyses for RNA-sequencing and
microarray studies. doi:
[10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007) - **RUVSeq** —
Normalization of RNA-seq data using factor analysis of control genes or
samples. doi: [10.1038/nbt.2931](https://doi.org/10.1038/nbt.2931)

**Differential expression & time series** - **DESeq2** — Moderated
estimation of fold change and dispersion for RNA-seq data with DESeq2.
doi:
[10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) -
**TiSA** — TiSA: TimeSeriesAnalysis — a pipeline for the analysis of
longitudinal transcriptomics data. doi:
[10.1093/nargab/lqad020](https://doi.org/10.1093/nargab/lqad020)

**Clustering & dimensionality reduction** - **WGCNA** — An R package for
weighted correlation network analysis. doi:
[10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559) -
**FactoMineR** — An R Package for Multivariate Analysis. doi:
[10.18637/jss.v025.i01](https://doi.org/10.18637/jss.v025.i01) -
**missMDA** (FAMD imputation) — missMDA: A Package for Handling Missing
Values in Multivariate Data Analysis. doi:
[10.18637/jss.v070.i01](https://doi.org/10.18637/jss.v070.i01) -
**PART** — Identifying clusters in genomics data by recursive
partitioning. doi:
[10.1515/sagmb-2013-0016](https://doi.org/10.1515/sagmb-2013-0016) -
**cluster** (Gower distance / PAM / silhouette) — Kaufman L, Rousseeuw
PJ. *Finding Groups in Data: An Introduction to Cluster Analysis*.
Wiley, 1990. ISBN 978-0-471-87876-6 (book, no DOI). - **grpreg** (group
lasso, discriminative analysis) — Group descent algorithms for nonconvex
penalized linear and logistic regression models with grouped predictors.
doi:
[10.1007/s11222-013-9424-2](https://doi.org/10.1007/s11222-013-9424-2) -
**nnet** (`multinom()`, multinomial logistic regression) — Venables WN,
Ripley BD. *Modern Applied Statistics with S*, 4th ed. Springer, 2002.
ISBN 0-387-95457-0 (book, no DOI).

**Deconvolution** - **CIBERSORT** *(source code must be obtained
separately)* — Profiling Tumor Infiltrating Immune Cells with CIBERSORT.
doi:
[10.1007/978-1-4939-7493-1_12](https://doi.org/10.1007/978-1-4939-7493-1_12)

**Annotation & enrichment** - **HOMER** — Simple combinations of
lineage-determining transcription factors prime cis-regulatory elements
required for macrophage and B cell identities. doi:
[10.1016/j.molcel.2010.05.004](https://doi.org/10.1016/j.molcel.2010.05.004) -
**ChIPseeker** — An R/Bioconductor package for ChIP peak annotation,
comparison and visualization. doi:
[10.1093/bioinformatics/btv145](https://doi.org/10.1093/bioinformatics/btv145) -
**clusterProfiler** — clusterProfiler 4.0: A universal enrichment tool
for interpreting omics data. doi:
[10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141) -
**gprofiler2** — An R package for gene list functional enrichment
analysis and namespace conversion toolset g:Profiler. doi:
[10.12688/f1000research.24956.2](https://doi.org/10.12688/f1000research.24956.2) -
**fgsea** — Fast gene set enrichment analysis. bioRxiv. doi:
[10.1101/060012](https://doi.org/10.1101/060012) - **MSigDB** (via
msigdbr) — Gene set enrichment analysis: a knowledge-based approach for
interpreting genome-wide expression profiles. doi:
[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102) -
**rrvgo** *(built on GOSemSim)* — GOSemSim: an R package for measuring
semantic similarity among GO terms and gene products. doi:
[10.1093/bioinformatics/btq064](https://doi.org/10.1093/bioinformatics/btq064) -
**EnrichedHeatmap** — EnrichedHeatmap: an R/Bioconductor package for
comprehensive visualization of genomic signal associations. doi:
[10.1186/s12864-018-4625-x](https://doi.org/10.1186/s12864-018-4625-x) -
**OmnipathR** — Integrated knowledgebase for multi-omics analysis. doi:
[10.1093/nar/gkaf1126](https://doi.org/10.1093/nar/gkaf1126) -
**DecoupleR** — Ensemble of computational methods to infer biological
activities from omics data. doi:
[10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)

**Normalization & imputation (proteomics)** - **vsn** — Variance
stabilization applied to microarray data calibration and to the
quantification of differential expression. doi:
[10.1093/bioinformatics/18.suppl_1.S96](https://doi.org/10.1093/bioinformatics/18.suppl_1.S96) -
**imputeLCMD** (QRILC) — Accounting for the Multiple Natures of Missing
Values in Label-Free Quantitative Proteomics Data Sets to Compare
Imputation Strategies. doi:
[10.1021/acs.jproteome.5b00981](https://doi.org/10.1021/acs.jproteome.5b00981)

**Isoform analysis** - **IsoformSwitchAnalyzeR** —
IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of
alternative splicing and its functional consequences. doi:
[10.1093/bioinformatics/bty1054](https://doi.org/10.1093/bioinformatics/bty1054)
