# HORIZON

Upstream bioinformatics processing package for bulk sequencing data. Handles QC, adapter trimming, alignment, BAM processing, read counting, BigWig generation, and chromatin accessibility / ChIP-seq / CUT&TAG peak calling — the computational work that sits upstream of GAIA.

Part of the [L-Atelier](../) repository. Unlike GAIA and CAULDRON, HORIZON has no named sub-modules; the wrapped tools serve that role implicitly. HORIZON has no hard dependency on GAIA or CAULDRON, and vice versa — the interface is file format compatibility.

---

Output files are in formats GAIA's import module already understands:

| HORIZON output | GAIA import |
|---|---|
| `count_matrix.csv` / `.rds` | `ELEUTHIA_load_counts_from_dir()` |
| `*_fragments.bed` | `ELEUTHIA_call_regions_from_fragments()` |
| `*_peaks.narrowPeak` | `ELEUTHIA_load_peaks()` |
| `*_RPKM.bw` / `*_spikein.bw` | `AETHER_plot_coverage_tracks()` |

---

## Installation

### 1. Install the R package

R >= 4.1.0 is required.

```r
install.packages("devtools")
devtools::install("path/to/L-Atelier/HORIZON")
library(HORIZON)
```

Core R dependencies (`Rsubread`, `Rfastp`, `Rsamtools`) are in `Imports` and installed automatically. `bowtie2` is called directly via the conda environment (see step 2).

### 2. Create the HORIZON conda environment

The chromatin pipeline (ATAC-seq, ChIP-seq, CUT&TAG/CUT&RUN) calls external CLI tools that must be available in a conda environment. Create it once:

```bash
conda create -n horizon_cli -c bioconda -c conda-forge \
  samtools=1.23.1 bedtools=2.31.1 macs3=3.0.2 bowtie2
```

> **Note:** deeptools may need to be installed via pip from within the environment if the bioconda build is unavailable or conflicts:
> ```bash
> conda activate horizon_cli
> pip install deeptools==3.5.5
> ```

### 3. Register the environment at the start of each session

```r
HORIZON_set_conda_env("~/miniconda3/envs/horizon_cli")
```

This stores the path for the session. Call it once at the top of any project script that uses the chromatin pipeline. RNA-seq functions (QC, align, count) do not require the conda environment.

---

## Output Directory Structure

All functions write into a consistent tree under `output_dir`:

```
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
```

---

## Functions

### Sample sheet
- `HORIZON_create_sample_sheet()` — generate a CSV template with all required columns
- `HORIZON_validate_sample_sheet()` — validate paths and column values

### Index building
- `HORIZON_build_index()` — build a Subread genome index (RNA-seq)
- `HORIZON_build_bowtie2_index()` — build a Bowtie2 index (chromatin assays)
- `HORIZON_build_combined_index()` — concatenate host + spike-in FASTAs and build a combined Bowtie2 index

### QC and trimming
- `HORIZON_run_qc_trim()` — QC + adapter trimming (wraps fastp via Rfastp)

### Alignment
- `HORIZON_run_align()` — splice-aware alignment for RNA-seq (wraps Rsubread)
- `HORIZON_run_bowtie2()` — paired-end alignment for chromatin assays (bowtie2 piped to samtools view; no SAM intermediate); MAPQ filtering included

### BAM processing
- `HORIZON_sort_index_bam()` — sort and index BAM (RNA-seq; wraps Rsamtools)
- `HORIZON_process_bam()` — collate → fixmate → sort → markdup (duplicate removal) → optional chrM removal; writes flagstat
- `HORIZON_filter_blacklist()` — remove blacklist-overlapping reads with bedtools
- `HORIZON_downsample_bam()` — subsample BAM by a spike-in scale factor (samtools)
- `HORIZON_rename_chromosomes()` — rename chromosome identifiers in BAM / BED / BigWig files; built-in T2T CHM13 map

### Spike-in normalisation
- `HORIZON_separate_spike_in()` — split a combined-genome BAM into host and spike-in BAMs by chromosome prefix
- `HORIZON_compute_spike_in_factors()` — count spike-in read pairs and compute scale factors for BigWig (`scale_factor_bw`) and DESeq2 (`size_factor_deseq2`)

### Format conversion and BigWig
- `HORIZON_bam_to_bed()` — BAM → fragment BED via bedtools; optional Tn5 shift for ATAC-seq
- `HORIZON_bam_to_bigwig()` — RPKM or spike-in normalised BigWig via deeptools `bamCoverage`; built-in effective genome sizes for GRCh38, GRCm38, T2TCHM13, WBcel235
- `HORIZON_run_fragment_size()` — insert size histogram via deeptools `bamPEFragmentSize`; nucleosomal periodicity diagnostic for ATAC-seq

### Peak calling
- `HORIZON_call_peaks()` — wraps `macs3 callpeak`; supports ATAC-seq, ChIP-seq (TF and histone), CUT&TAG, CUT&RUN via explicit parameter control

### Counting and aggregation (RNA-seq)
- `HORIZON_run_count()` — gene-level feature counting (wraps featureCounts via Rsubread); exon-union counting, not suitable for isoform-level resolution (overlapping-isoform reads are ambiguous/discarded)
- `HORIZON_aggregate_counts()` — merge per-sample count files into a single matrix

### Transcript quantification (isoform-level RNA-seq)
- `HORIZON_build_salmon_index()` — decoy-aware Salmon transcriptome index from a GTF + genome FASTA (`GenomicFeatures::extractTranscriptSeqs()`); one-time per reference
- `HORIZON_run_salmon()` — per-sample `salmon quant` (selective-alignment mode, direct from FASTQ — does not use the genome-aligned BAM); resolves isoform-ambiguous reads via Salmon's EM algorithm, unlike `HORIZON_run_count()`
- `HORIZON_aggregate_salmon()` — merge per-sample `quant.sf` files into transcript x sample TPM and count matrices

### Transposable element quantification
Reads from repetitive TE loci mostly multi-map, so `HORIZON_run_align()` + `HORIZON_run_count()` (tuned to keep one best hit per read, correct for gene-level counting) systematically undercounts TEs. This path uses STAR with a relaxed multi-mapping filter, then TEcount's EM algorithm to redistribute those reads across candidate TE loci, rather than discarding them.

- `HORIZON_build_star_index()` — build a STAR genome index (separate from the Rsubread index used by `HORIZON_run_align()`); one-time per reference/read-length
- `HORIZON_run_star_align()` — STAR alignment with multi-mapping alignments reported (`--outFilterMultimapNmax`/`--winAnchorMultimapNmax` relaxed from STAR's defaults); writes to `aligned_star/`, distinct from Rsubread's `aligned/`
- `HORIZON_run_tecount()` — per-sample `TEcount` (bioconda `tetranscripts` package's single-BAM counting tool, *not* the 2-group `TEtranscripts` DE binary) quantifying genes + TE subfamilies from a STAR BAM
- `HORIZON_aggregate_tecount()` — merge per-sample `.cntTable` files into a TE subfamily x sample count matrix (TE rows only; gene rows are dropped — canonical gene counts should keep coming from `HORIZON_aggregate_counts()`)
- `HORIZON_merge_gene_te_counts()` — combine a gene-level count matrix and a TE count matrix into one feature x sample matrix for differential analysis (e.g. `ARTEMIS_normalize_counts()` in GAIA)

STAR and TEcount share the standard `horizon_cli` conda environment — no dependency conflict with the tools already there (confirmed via `conda install --dry-run`), unlike PARSE split-pipe which genuinely needs its own isolated environment:
```bash
conda install -n horizon_cli -c bioconda -c conda-forge star tetranscripts
```
Both functions call `.horizon_run_cli()` internally, so `HORIZON_set_conda_env()` must be called first (same as any other `horizon_cli` tool) — no separate `conda_env` argument is needed.

---

## Typical Chromatin Workflow

```
HORIZON_build_combined_index()          # once per reference pair
  ↓ per sample
HORIZON_run_qc_trim()
HORIZON_run_bowtie2()
HORIZON_separate_spike_in()             # spike-in workflows only
HORIZON_process_bam()
HORIZON_filter_blacklist()
HORIZON_run_fragment_size()             # QC diagnostic
  ↓ spike-in only: compute factors, then downsample
HORIZON_compute_spike_in_factors()
HORIZON_downsample_bam()
  ↓ all samples
HORIZON_bam_to_bed()
HORIZON_bam_to_bigwig()
HORIZON_call_peaks()
```

### Recommended MACS3 parameters by assay type

| Assay | Key parameters |
|---|---|
| ATAC-seq | `nomodel=TRUE, shift=-100, extsize=200` |
| ChIP-seq (TF) | defaults; provide `control_bed` |
| ChIP-seq (histone) | `broad=TRUE`; provide `control_bed` |
| CUT&TAG / CUT&RUN | `nomodel=TRUE`; optionally `extra_flags=c("--keep-dup","all")` |

---

## Development

```r
devtools::document("HORIZON")
devtools::install("HORIZON")
```

Implementation status is tracked in `STATUS/HORIZON_STATUS.txt` at the repository root.
