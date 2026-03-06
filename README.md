# ZERO_DAWN

A private multi-omics analysis toolkit for RNAseq, ATACseq, CUT&TAG, and proteomics data. The repository houses three companion R packages that together cover the full analysis lifecycle — from raw sequencing files to integrated multi-omics results.

This toolkit was developed with the aid of Claude Code as a means to test and familiarise with a new coding workflow.

> **Naming note:** The repository name ZERO_DAWN is temporary and will be updated before any public release. Package names (GAIA, HORIZON, CAULDRON) are permanent.

---

## Packages

### [GAIA](GAIA/) — Multi-omics Analysis

The primary analysis package. Covers the full downstream analytical workflow for bulk multi-omics data: quality control, preprocessing, annotation, enrichment, integration, machine learning, statistical analysis, visualization, and reporting. Named after the central AI system from *Horizon Zero Dawn*; its nine subordinate modules are named after figures from Greek mythology.

→ See [`GAIA/README.md`](GAIA/README.md) for installation, modules, and references.

### [CAULDRON](CAULDRON/) — Single-Cell Analysis *(early development)*

The single-cell and single-nuclei counterpart to GAIA, designed for scRNA-seq and snRNA-seq analysis. Built on `SingleCellExperiment` (Bioconductor). Follows the same modular architecture as GAIA, with subordinate modules named after Hephaestus's creations from Greek mythology. Currently in early development.

→ See [`CAULDRON/README.md`](CAULDRON/README.md) for the planned architecture.

### [HORIZON](HORIZON/) — Upstream Processing

Handles the computational work that sits upstream of GAIA and CAULDRON: QC, adapter trimming, alignment, BAM processing, read counting, and BigWig generation. Designed to produce outputs in formats that GAIA already knows how to import. A flat collection of wrappers with no sub-modules.

→ See [`HORIZON/README.md`](HORIZON/README.md) for installation and the processing workflow.

---

## Repository Structure

```
ZERO_DAWN/
├── GAIA/               # Bulk multi-omics analysis package
├── HORIZON/            # Upstream processing package
├── CAULDRON/           # Single-cell analysis package (early development)
├── STATUS/             # Per-module implementation tracking files
├── projects/           # Individual analysis projects (git-ignored)
├── data/               # Shared reference data (git-ignored)
└── to_process/         # Development notes and concept documents
```

---

## Issues

Feel free to open an issue for any found bugs or suggested improvements.

## Personal notes

Use of cat() instead of message() in order to enable console output and logging at the same time.
