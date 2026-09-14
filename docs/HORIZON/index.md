# HORIZON

Upstream bioinformatics processing package for bulk sequencing data.
Handles QC, adapter trimming, alignment, BAM processing, read counting,
BigWig generation, and chromatin accessibility / ChIP-seq / CUT&TAG peak
calling. This upstream processing is also used in regards to single cell
data - Making the output of HORIZON compatable with both GAIA and
CAULDRON.

Part of the [L-Atelier](https://ylefol.github.io/L-Atelier/) repository.
Unlike GAIA and CAULDRON, HORIZON has no named sub-modules; the wrapped
tools serve that role implicitly. HORIZON has no hard dependency on GAIA
or CAULDRON, and vice versa — the interface is file format
compatibility.

> **Status:** This package, like the others, is under active
> development. Some elements are unlikely to change while others are
> very prone to change.

------------------------------------------------------------------------

## Installation

See the
[Installation](https://ylefol.github.io/L-Atelier/HORIZON/articles/installation.md)
guide for the R package, the `horizon_cli` conda environment, and the
optional PARSE Biosciences environment.

------------------------------------------------------------------------

Output files are in formats GAIA’s import module already understands:

| HORIZON output               | GAIA import                              |
|------------------------------|------------------------------------------|
| `count_matrix.csv` / `.rds`  | `ELEUTHIA_load_counts_from_dir()`        |
| `*_fragments.bed`            | `ELEUTHIA_call_regions_from_fragments()` |
| `*_peaks.narrowPeak`         | `ELEUTHIA_load_peaks()`                  |
| `*_RPKM.bw` / `*_spikein.bw` | `AETHER_plot_coverage_tracks()`          |

------------------------------------------------------------------------

## Output Directory Structure

All functions write into a consistent tree under `output_dir`:

    <output_dir>/
    └── <sample_id>/
        ├── qc/          # Trimmed FASTQs + fastp HTML/JSON report
        ├── aligned/     # BAM files, BAI indexes, BigWig, fragment size plot
        ├── aligned_star/ # STAR-aligned BAM (multi-mapping reported; TE quantification)
        ├── counts/      # featureCounts gene-level counts (RNA-seq)
        ├── salmon/      # Salmon quant.sf (isoform-level RNA-seq)
        ├── tecount/     # TEcount .cntTable (transposable element quantification)
        └── peaks/       # MACS3 peak files (narrowPeak / broadPeak)
    <output_dir>/aggregated/counts/   # Combined gene count matrix + metadata (RNA-seq)
    <output_dir>/aggregated/salmon/   # Combined transcript TPM/count matrices + metadata (RNA-seq)
    <output_dir>/aggregated/tecount/  # Combined TE subfamily count matrix + metadata (RNA-seq)

------------------------------------------------------------------------

## Development

Implementation status is tracked in `STATUS/HORIZON_STATUS.txt` at the
repository root.

------------------------------------------------------------------------

## References

Methods/tools wrapped or reimplemented within HORIZON:

- **fastp** (via Rfastp) — fastp: an ultra-fast all-in-one FASTQ
  preprocessor. doi:
  [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
- **Rsubread / featureCounts** — The R package Rsubread is easier,
  faster, cheaper and better for alignment and quantification of RNA
  sequencing reads. doi:
  [10.1093/nar/gkz114](https://doi.org/10.1093/nar/gkz114)
- **Bowtie2** — Fast gapped-read alignment with Bowtie 2. doi:
  [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)
- **samtools** — Twelve years of SAMtools and BCFtools. doi:
  [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
- **bedtools** — BEDTools: a flexible suite of utilities for comparing
  genomic features. doi:
  [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)
- **MACS3/MACS2** — Model-based Analysis of ChIP-Seq (MACS). doi:
  [10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)
- **deeptools** — deepTools2: a next generation web server for
  deep-sequencing data analysis. doi:
  [10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257)
- **Salmon** — Salmon provides fast and bias-aware quantification of
  transcript expression. doi:
  [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)
- **STAR** — STAR: ultrafast universal RNA-seq aligner. doi:
  [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)
- **TEtranscripts/TEcount** — TEtranscripts: a package for including
  transposable elements in differential expression analysis of RNA-seq
  datasets. doi:
  [10.1093/bioinformatics/btv422](https://doi.org/10.1093/bioinformatics/btv422)
- **SEACR** *(method reimplemented in R by
  [`HORIZON_call_peaks_seacr()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_call_peaks_seacr.md),
  not called as an external tool)* — Peak calling by Sparse Enrichment
  Analysis for CUT&RUN chromatin profiling. doi:
  [10.1186/s13072-019-0287-4](https://doi.org/10.1186/s13072-019-0287-4).
  Original tool:
  [github.com/FredHutch/SEACR](https://github.com/FredHutch/SEACR).
- **scVelo** *(target format for
  [`HORIZON_run_parse_velocity()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_parse_velocity.md)’s
  spliced/unspliced AnnData output; the counting step itself is a custom
  script, not a scVelo wrapper)* — Generalizing RNA velocity to
  transient cell states through dynamical modeling. doi:
  [10.1038/s41587-020-0591-5](https://doi.org/10.1038/s41587-020-0591-5)
- **PARSE Biosciences `split-pipe`** — proprietary vendor software; see
  parsebiosciences website for code/software.
