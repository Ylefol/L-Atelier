# L-Atelier

L-Atelier (French for *the workshop*) is a multi-omics toolkit that I've been designing to simplify my own workflows. It is composed of four separate packages, each designed to fulfill some aspects of bioinformatic analyses ranging from early data processing and QC to data analysis and figures. Much of the naming convention is derived from the Horizon video game series, and the project has doubled as my way of learning to work with Claude Code.

> **Status:** This entire repository is under-active development.

---

## Packages

### [GAIA](GAIA/) — Multi-omics Analysis

Covers the full downstream analytical workflow for bulk multi-omics data: quality control, preprocessing, annotation, enrichment, integration, machine learning, statistical analysis, visualization, reporting ect... Its nine subordinate modules are named after figures from Greek mythology.

→ See [`GAIA/README.md`](GAIA/README.md) for installation, modules, and references.

### [CAULDRON](CAULDRON/) — Single-Cell Analysis *(early development)*

The single-cell and single-nuclei counterpart to GAIA, currently focused on scRNA-seq and snRNA-seq analysis. Built on `SingleCellExperiment` (Bioconductor). Follows the same modular architecture as GAIA, with subordinate modules named after Hephaestus's creations from Greek mythology.

→ See [`CAULDRON/README.md`](CAULDRON/README.md) for installation, modules, and references.

### [HORIZON](HORIZON/) — Upstream Processing

Handles the computational work that sits upstream of GAIA and CAULDRON: QC, adapter trimming, alignment, BAM processing, read counting, BigWig generation ect... Designed to produce outputs in formats that GAIA and CAULDRON already know how to import. A flat collection of wrappers with no sub-modules; all functions are prefixed `HORIZON_`

→ See [`HORIZON/README.md`](HORIZON/README.md) for installation and references.

### [CYAN](CYAN/) — Niche Genomic Analyses

A focused package for specialised genomic analyses that sit outside the scope of the other three packages. Currently implements RNA editing detection via REDItools2 (basilisk-managed Python environment, system2 invocation). Somatic variant calling via GATK MuTect2 is planned. Like HORIZON, CYAN is a flat collection of functions with no sub-modules; all functions are prefixed `CYAN_`.

→ See [`CYAN/README.md`](CYAN/README.md) for installation and references; [`CYAN/inst/status/CYAN_STATUS.txt`](CYAN/inst/status/CYAN_STATUS.txt) for current implementation status.

---

## Repository Structure

```
L-Atelier/
├── GAIA/               # Bulk multi-omics analysis package
├── HORIZON/            # Upstream processing package
├── CAULDRON/           # Single-cell analysis package (early development)
├── CYAN/               # Niche genomic analyses (RNA editing, variant calling)
├── STATUS/             # Per-module implementation tracking files
├── projects/           # Individual analysis projects (git-ignored)
├── data/               # Shared reference data (git-ignored)
└── to_process/         # Development notes and concept documents
```

---

## Design Philosophy

I want to emphasize that this is designed as my own personal toolkit, thus I've tailored the packages to suit my needs — i.e. exposing certain parameters while hiding others. Of course, eventual users may still open issues and I'll do my best to answer them quickly.

I've decided to split this into four packages to ease the burden of development and possible installation. In addition, this was mostly developed in R, mainly due to bioinformatics' dependency on R, but also due to R's friendly documentation practices. For this reason, certain elements which would require a Python script are still implemented in R and called via an R wrapper. The same applies to `system2` calls. This way, everything can be run through R.

---

## Known Issues

- **IDE compatibility with `system2` calls**: some IDEs don't play well with R scripts that use `system2` (Positron's R console has been confirmed to have issues, for example). If a script fails in a way that looks environment-related rather than a real code bug, try running it from the console/terminal directly instead before digging further.
- **IDE compatibility with some downloads**: I've observed a similar issue with certain downloads not working when called via an idea (such as an annotation database). I've yet to pin down exactly why some downloads fail this way and others do not. If you find yourself experience download related issues within the IDE there is a chance you can fix it by calling the script via the console (Rscript).

---

## Platform

This toolkit has only been developed and tested on Ubuntu 24.04 LTS (Noble Numbat). Cross-platform compatibility (macOS, Windows, other Linux distributions) hasn't been verified — conda/basilisk environment setup, `system2` calls, and file path handling in particular may behave differently elsewhere. If you're on a different OS, expect to do some of your own troubleshooting.

---

## License

All four packages (GAIA, HORIZON, CAULDRON, CYAN) are MIT-licensed — see each package's `LICENSE.md`. CYAN vendors one third-party file (`CYAN/inst/python/reditools.py`, REDItools2, unmodified) under its own CC BY 4.0 license; see the attribution notice at the top of that file and `CYAN/LICENSE.md` for details. This repository's own MIT license does not cover external tools GAIA/HORIZON/CYAN merely call or depend on (CIBERSORT, HOMER, GATK, MACS3, etc.) — those remain governed by their own licenses; CIBERSORT's source in particular must be obtained separately from its authors (see `GAIA/README.md`). The same is applicable to the spipe software (processing of parse biosciences single cell sequencing data).
