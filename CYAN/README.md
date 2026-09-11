# CYAN

Focused package for niche genomic analyses that sit outside the scope of the main GAIA/HORIZON/CAULDRON pipeline. Currently implements RNA editing detection via REDItools2. Somatic variant calling via GATK MuTect2 is planned but not yet implemented. This package has seen very little use and should be treated as unreliable for the time being. As the need arises, I'll develop this further.

Part of the [L-Atelier](../) repository. Like HORIZON, CYAN is a flat collection of functions with no sub-modules; all public functions are prefixed `CYAN_`. Design principles: one sample per call (cohort looping left to the user), large intermediates written to `work_dir`, final outputs to `output_dir`.

---

## Installation

```r
install.packages("devtools")
devtools::install("path/to/L-Atelier/CYAN")
library(CYAN)
```

R >= 4.3 is required. `basilisk` and `reticulate` are in `Imports` and installed automatically.

No manual conda setup is needed for RNA editing: `CYAN_run_reditools()` runs REDItools2 through a `basilisk`-managed Python environment (Python 3.10.14, `pysam`, `sortedcontainers`, `psutil`, `netifaces`, pinned versions) that `basilisk` creates automatically on first use. REDItools2 itself is vendored in `inst/python/reditools.py` (the original script, unmodified).

---

## Functions

### Validation utilities
- `CYAN_check_bam()` — validate a BAM file exists and has a `.bai` index
- `CYAN_check_references()` — validate a reference FASTA has a `.fai` index; optionally check additional VCF/BED paths exist
- `CYAN_resolve_gatk()` — resolve the GATK executable (explicit path → `GATK_PATH` env var → `PATH`); not yet called by any function since MuTect2 isn't implemented

### RNA editing detection
- `CYAN_run_reditools()` — detect RNA editing events (e.g. A-to-I) in a single BAM via REDItools2. Checkpoints raw output before parsing (re-running skips REDItools2 if the checkpoint exists); defaults match REDItools2's own native defaults (`min_coverage=1`, `min_edits=0`, `strict=FALSE`, unfiltered output for post-hoc filtering). Guards against a common REDItools2 pitfall: its duplicate filter only excludes reads that already carry the SAM duplicate flag (0x400) — on a BAM that was never `markdup`'d, that filter is a silent no-op, so `require_dup_marked=TRUE` (default) hard-stops unless at least one duplicate-flagged read is found first.
- `CYAN_run_reditools_batch()` — convenience wrapper iterating `CYAN_run_reditools()` over a vector of BAMs; per-sample error handling (a failed sample warns and continues rather than aborting the batch), inherits the same checkpoint behaviour.

---

## Output

`CYAN_run_reditools()` writes to `output_dir`:
- `{sample}_{date}_reditools.txt` — full raw REDItools2 table (all covered positions)
- `{sample}_{date}_reditools_subs.txt` — substitution sites only
- `{sample}_{date}_reditools.rds` — parsed `cyan_editing_result` object (opt-out with `save_rds = FALSE`)

The parsed `cyan_editing_result` is a `data.frame` (one row per genomic position: `region`, `position`, `reference`, `strand`, `coverage`, `mean_quality`, per-base counts, `all_subs`, `frequency`), with a `params` attribute recording the call's arguments for reproducibility.

---

## References

- **REDItools2** — Investigating RNA editing in deep transcriptome datasets with REDItools and REDIportal. doi: [10.1038/s41596-019-0279-7](https://doi.org/10.1038/s41596-019-0279-7)

---

## Development

Implementation status is tracked in `inst/status/CYAN_STATUS.txt`, following the same format as GAIA's `STATUS/` files at the repository root.
