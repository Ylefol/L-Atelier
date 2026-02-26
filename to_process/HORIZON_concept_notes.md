# HORIZON — Concept Notes

## What is it?
A hypothetical affiliate R package for upstream bioinformatics processing.
Handles the heavy computational work (QC, trimming, alignment, counting) that
sits upstream of GAIA and CAULDRON. Not part of either package, but designed
to feed into both.

## Name Rationale
HORIZON references the Horizon Zero Dawn franchise directly. The horizon is
where the journey begins — raw sequencing data is the starting point before
it enters the analysis ecosystem. The data flow narrative is intentionally
clean: HORIZON → GAIA (bulk) / CAULDRON (single-cell). The name is simple,
immediately intuitive without HZD knowledge, and suits a foundational package
that does not warrant the complexity of named sub-modules.

## Structure
HORIZON has no sub-modules. The wrapped tools themselves (fastp, bowtie2,
HISAT2, featureCounts, etc.) serve that role implicitly. The package is a flat
collection of wrappers and format conversion functions — one layer of functions
that call the underlying Bioconductor tool packages and return outputs in
formats that GAIA and CAULDRON already expect.

## Scope
Raw sequencing data (FASTQ) → processed, analysis-ready files that GAIA's
ELEUTHIA module and CAULDRON's TALARIA module already know how to import.

Coverage:
- QC (read quality, adapter content, duplication)
- Trimming / adapter removal
- Alignment (RNA-seq, ChIP-seq, ATAC-seq, CUT&TAG)
- BAM processing (sorting, indexing, filtering)
- Read counting / peak-level quantification
- BED / BigWig generation

## R Package Dependencies (candidates)
| Task                        | Package                 | Notes                                   |
|-----------------------------|-------------------------|-----------------------------------------|
| QC                          | ShortRead / Rqc         | Native                                  |
| Trimming                    | Rfastp                  | Bundles fastp binary, cross-platform    |
| RNA-seq alignment           | Rhisat2 / Rsubread      | Both splice-aware; no R wrapper for STAR|
| ChIP/ATAC/CUT&TAG alignment | Rbowtie2                | Bundles bowtie2 binary, cross-platform  |
| BAM handling                | Rsamtools               | Mature, widely used                     |
| Counting                    | Rsubread::featureCounts | Same tool as command-line featureCounts |
| BED/BigWig I/O              | rtracklayer             | Native                                  |

## Cross-Platform Feasibility
Bioconductor packages bundle their compiled binaries within the package itself,
so installation via BiocManager::install() requires no external tools on PATH.
This makes the pipeline genuinely cross-platform:
- Linux / macOS: fully feasible
- Windows: possible but some packages have historical friction (worth verifying)

## The STAR Gap
There is no R package that wraps STAR. For RNA-seq alignment, HORIZON uses
HISAT2 (Rhisat2) or Subread (Rsubread::align()) instead. Both are splice-aware
and produce comparable results for standard use cases. The gap matters more for
edge cases (novel junction discovery, fusion detection).

## Interface with GAIA and CAULDRON
The interface is file format compatibility, not a shared object model.
HORIZON's outputs are the file formats that downstream packages already import:
- Count matrices (CSV / RDS) → GAIA ELEUTHIA / CAULDRON TALARIA
- BED files → GAIA ELEUTHIA
- BigWig files → GAIA ELEUTHIA / AETHER

HORIZON has no hard dependency on GAIA or CAULDRON and vice versa.

## Status
Hypothetical — no implementation planned yet. Notes captured for future reference.
