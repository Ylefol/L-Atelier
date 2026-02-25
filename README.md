# ZERO_DAWN

A private multi-omics integration toolkit a variety of omics - primarily focused on RNAseq, ATACseq, CHIPseq, and proteomics. Note that no single-cell sequencing is present in this repository.
This toolkit was developed with the aid of CLAUDE code as a means for me to test and familiarize myself with a new style of coding.

> **Naming note:** This repository is called ZERO_DAWN, but the installable R package within it is called **GAIA**. The repository name is temporary and will be updated before any public release; the package name GAIA is permanent.

---

## Installation

GAIA is an R package installed directly from the repository. R >= 4.1.0 is required.

```r
# Install devtools if not already available
install.packages("devtools")

# Install GAIA from the local repository
devtools::install("path/to/ZERO_DAWN/GAIA")

# Load the package
library(GAIA)
```

Loading `library(GAIA)` attaches `ggplot2` and loads the namespaces of `DESeq2`, `SummarizedExperiment`, and `WGCNA`.

### Optional dependencies

Several modules have optional dependencies that are loaded when needed (listed in `Suggests`). Key ones to be aware of:

- **CIBERSORT** — Source code must be obtained separately from the authors. See the [CIBERSORT reference](#references) below.
- **Bioconductor packages** — `GenomicRanges`, `IRanges`, `ChIPseeker`, `clusterProfiler`, `Rsubread`, etc. Install via `BiocManager::install()`.
- **Batch correction** — `sva` (ComBat), `limma` (removeBatchEffect).
- **Annotation** — `org.Hs.eg.db`, `AnnotationDbi`.
- **Quantile normalization** — `preprocessCore`.

---

## Project Structure

```
ZERO_DAWN/
├── GAIA/               # R package (install with devtools::install("GAIA"))
│   ├── R/              # All module source files (flat, MODULE_prefixed names)
│   ├── man/            # Auto-generated roxygen2 documentation
│   ├── DESCRIPTION
│   └── NAMESPACE
├── projects/           # Individual analysis projects (git-ignored)
├── data/               # Shared reference data (git-ignored)
└── to_process/         # Development notes and module STATUS files
    └── STATUS/         # Per-module implementation status tracking
```

All module source files live flat under `GAIA/R/` with module-prefixed names (e.g., `artemis_wgcna.R`, `aether_circos.R`). Function names follow the same convention: `MODULE_function_name` (e.g., `ARTEMIS_wgcna_detect_modules`, `AETHER_plot_circos`).

Note an YGGDRASIL package is in the works - this will be a package formatted similarly to GAIA intended for single cell/nuclei analyses.
---

## Issues
Feel free to open an issue for any found bugs or suggested improvements to the toolkit.

---

## GAIA Modules

Inspired by Horizon Zero Dawn, GAIA is the central system orchestrating specialized subordinate modules. Each module handles a specific stage of the analysis pipeline.

### Hades — Quality Control & Filtering

Named after the deity of the underworld, Hades removes the "dead" — filtering out low-quality data and problematic samples before they can contaminate downstream analysis.

### Poseidon — Preprocessing & Normalization

Like the god of the seas who controlled and purified waters, Poseidon cleanses and normalizes raw data, washing away technical noise and batch effects.

### Apollo — Annotation & Enrichment

Apollo, god of knowledge and prophecy, brings wisdom to the data through annotation — linking peaks to genes, identifying motifs, and revealing biological pathways.

### Hephaestus — Integration

The master craftsman who forges disparate materials into powerful artifacts, Hephaestus combines multiple omics datasets into unified, coherent insights.

### Minerva — Machine Learning & Validation

The Roman goddess of wisdom and strategic warfare, Minerva employs sophisticated validation strategies — cross-validation, hyperparameter tuning, and model evaluation.

### Artemis — Analysis & Statistics

Goddess of the hunt, Artemis precisely targets and captures statistical insights with the accuracy of her legendary arrows.

### Aether — Visualization

Named for the pure upper air breathed by gods, Aether clarifies and illuminates results, making the invisible visible through visualization.

### Eleuthia — Import, Export & Reporting

The cradle facility that brought new life into the world, Eleuthia handles the birth and delivery of data — importing raw inputs and delivering polished final reports.

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

Module implementation status is tracked in `to_process/STATUS/` (e.g., `Artemis_STATUS.txt`).

---

## References

The following methods are implemented or wrapped within GAIA:

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
