# HORIZON

Upstream bioinformatics processing package for bulk sequencing data. Handles QC, adapter trimming, alignment, BAM processing, read counting, and BigWig generation — the computational work that sits upstream of GAIA.

Part of the [ZERO_DAWN](../) repository. Unlike GAIA and CAULDRON, HORIZON has no named sub-modules; the wrapped tools (fastp, Subread, Rsamtools, etc.) serve that role implicitly. HORIZON has no hard dependency on GAIA or CAULDRON, and vice versa — the interface is file format compatibility.

---

Output files are in formats GAIA's import module already understands:

| HORIZON output | GAIA import |
|---|---|
| `count_matrix.csv` / `.rds` | `ELEUTHIA_load_counts_from_dir()` |
| `*_sorted.bam` | `ELEUTHIA_call_regions_from_fragments()` |
| `*_RPKM.bw` | `AETHER_plot_coverage_tracks()` |

---

## Installation

HORIZON is an R package installed directly from the repository. R >= 4.1.0 is required.

```r
install.packages("devtools")
devtools::install("path/to/ZERO_DAWN/HORIZON")
library(HORIZON)
```

Core dependencies (`Rsubread`, `Rfastp`, `Rsamtools`) are in `Imports` and installed automatically. BigWig generation (`HORIZON_bam_to_bigwig`) additionally requires deeptools on `PATH`:

```bash
conda install -c bioconda deeptools
# or
pip install deeptools
```

---

## Output Directory Structure

All functions write into a consistent tree under `output_dir`:

```
<output_dir>/
└── <sample_id>/
    ├── qc/          # Trimmed FASTQs + fastp HTML/JSON report
    ├── aligned/     # Sorted BAM, BAI index, optional RPKM BigWig
    └── counts/      # Per-sample counts.txt + featurecounts.rds
<output_dir>/aggregated/counts/   # Combined count matrix + metadata
```

---

## Functions

**Sample sheet**
- `HORIZON_create_sample_sheet()` — generate a CSV template with all required columns
- `HORIZON_validate_sample_sheet()` — validate paths and column values

**Index**
- `HORIZON_build_index()` — build a genome index (once per reference genome)

**Per-sample pipeline**
- `HORIZON_run_qc_trim()` — QC + adapter trimming (wraps fastp via Rfastp)
- `HORIZON_run_align()` — splice-aware alignment (wraps Rsubread)
- `HORIZON_sort_index_bam()` — sort and index BAM (wraps Rsamtools)
- `HORIZON_run_count()` — feature counting (wraps featureCounts via Rsubread)
- `HORIZON_bam_to_bigwig()` — RPKM-normalised BigWig via deeptools `bamCoverage`; requires `effective_genome_size`; default bin size 10 bp

**Aggregation**
- `HORIZON_aggregate_counts()` — merge per-sample count files into a single matrix

---

## Development

```r
devtools::document("HORIZON")
devtools::install("HORIZON")
```

Implementation status is tracked in `STATUS/HORIZON_STATUS.txt` at the repository root.
