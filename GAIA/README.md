# GAIA

Multi-omics integration toolkit for CUT&TAG, ATAC-seq, and RNA-seq analysis. Covers quality control, preprocessing, annotation, enrichment, integration, machine learning, statistical analysis, visualization, and reporting.

Part of the [L-Atelier](../) repository. Its nine subordinate modules are each named after a figure from Greek mythology.

---

## Installation

GAIA is an R package installed directly from the repository. R >= 4.1.0 is required.

```r
install.packages("devtools")
devtools::install("path/to/L-Atelier/GAIA")
library(GAIA)
```

### Optional dependencies

Several modules have optional dependencies loaded on demand (listed in `Suggests`). Key ones:

- **Bioconductor** — `GenomicRanges`, `IRanges`, `ChIPseeker`, `clusterProfiler`, `Rsubread`, `rtracklayer`, etc. Install via `BiocManager::install()`.
- **Batch correction** — `sva` (ComBat), `limma` (removeBatchEffect).
- **Annotation** — `org.Hs.eg.db`, `AnnotationDbi`.
- **CIBERSORT** — Source code must be obtained separately from the authors. See [References](#references).
- **Quantile normalization** — `preprocessCore`.

If Bioconductor packages are unexpectedly unavailable after installation, run `BiocManager::valid()` to identify and fix any version mismatches.

---

## Modules

GAIA orchestrates nine subordinate modules. 
Individual modules can also be used independently.

### Hades — Quality Control & Filtering

Named after the deity of the underworld — Hades removes the "dead", filtering out low-quality data and problematic samples before they can contaminate downstream analysis.

### Poseidon — Preprocessing & Normalization

Like the god of the seas who controlled and purified waters, Poseidon cleanses and normalizes raw data, washing away technical noise and batch effects.

### Apollo — Annotation & Enrichment

Apollo, god of knowledge and prophecy, brings wisdom to the data through annotation — linking peaks to genes, identifying transcription factor motifs, and revealing biological pathways.

### Hephaestus — Integration

The master craftsman who forges disparate materials into powerful artifacts, Hephaestus combines multiple omics datasets into unified, coherent insights through multi-view integration methods.

### Minerva — Machine Learning & Validation

The Roman goddess of wisdom and strategic warfare, Minerva employs sophisticated validation strategies: cross-validation, hyperparameter tuning, and model evaluation.

### Artemis — Analysis & Statistics

Goddess of the hunt, Artemis precisely targets and captures statistical insights — differential analysis, statistical testing, correlation analysis, and deconvolution.

### Aether — Visualization

Named for the pure upper air breathed by gods, Aether clarifies and illuminates results, generating publication-quality plots from coverage tracks to circos diagrams.

### Eleuthia — Import, Export & Reporting

The cradle facility that brought new life into the world, Eleuthia handles the birth and delivery of data — importing raw inputs and delivering polished outputs and reports.

### Demeter — Utilities

Goddess of agriculture and harvest, Demeter cultivates the foundation — providing the essential tools and helper functions that nourish all other modules.

---

## Development

After modifying any source file in `GAIA/R/`:

```r
devtools::document("GAIA")   # Regenerate NAMESPACE and man/ from roxygen2
devtools::install("GAIA")    # Reinstall the updated package
library(GAIA)                # Reload
```

All source files are flat under `GAIA/R/` with module-prefixed names (e.g., `artemis_wgcna.R`, `aether_circos.R`). Function names follow the same convention: `MODULE_function_name` (e.g., `ARTEMIS_wgcna_detect_modules`, `AETHER_plot_circos`).

Module implementation status is tracked in `STATUS/` at the repository root (e.g., `STATUS/Artemis_STATUS.txt`).

---

## References

Methods implemented or wrapped within GAIA:

- **sCCA** — Integrating multi-OMICS data through sparse canonical correlation analysis for the prediction of complex traits: a comparison study. doi: [10.1093/bioinformatics/btaa530](https://doi.org/10.1093/bioinformatics/btaa530)
- **TiSA** — TiSA: TimeSeriesAnalysis — a pipeline for the analysis of longitudinal transcriptomics data. doi: [10.1093/nargab/lqad020](https://doi.org/10.1093/nargab/lqad020)
- **DESeq2** — Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. doi: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
- **FactoMineR** — An R Package for Multivariate Analysis. doi: [10.18637/jss.v025.i01](https://doi.org/10.18637/jss.v025.i01)
- **WGCNA** — An R package for weighted correlation network analysis. doi: [10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)
- **PART** — Identifying clusters in genomics data by recursive partitioning. doi: [10.1515/sagmb-2013-0016](https://doi.org/10.1515/sagmb-2013-0016)
- **gprofiler2** — An R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler. doi: [10.12688/f1000research.24956.2](https://doi.org/10.12688/f1000research.24956.2)
- **CIBERSORT** *(source code must be obtained separately)* — Profiling Tumor Infiltrating Immune Cells with CIBERSORT. doi: [10.1007/978-1-4939-7493-1_12](https://doi.org/10.1007/978-1-4939-7493-1_12)
- **HOMER** — Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. doi: [10.1016/j.molcel.2010.05.004](https://doi.org/10.1016/j.molcel.2010.05.004)
- **OmnipathR** — Integrated knowledgebase for multi-omics analysis. doi: [10.1093/nar/gkaf1126](https://doi.org/10.1093/nar/gkaf1126)
- **DecoupleR** — Ensemble of computational methods to infer biological activities from omics data. doi: [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)
- **ChIPseeker** — An R/Bioconductor package for ChIP peak annotation, comparison and visualization. doi: [10.1093/bioinformatics/btv145](https://doi.org/10.1093/bioinformatics/btv145)
