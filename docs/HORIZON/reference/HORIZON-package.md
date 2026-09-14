# HORIZON: Upstream Bioinformatics Processing

Upstream preprocessing pipeline for bulk RNA-seq data. Wraps Rfastp,
Rsubread, and Rsamtools into a consistent, RAM-aware interface designed
to feed downstream analysis in GAIA and CAULDRON.

Typical workflow:

1.  [`HORIZON_create_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_create_sample_sheet.md)
    — generate a sample sheet template

2.  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md)
    — validate paths and columns

3.  [`HORIZON_build_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_index.md)
    — build genome index (once per reference)

4.  [`HORIZON_run_qc_trim`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md)
    — QC + adapter trimming (per sample)

5.  [`HORIZON_run_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_align.md)
    — splice-aware alignment (per sample)

6.  [`HORIZON_sort_index_bam`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_sort_index_bam.md)
    — sort and index BAM (per sample)

7.  [`HORIZON_run_count`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_count.md)
    — feature counting (per sample)

8.  [`HORIZON_aggregate_counts`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_counts.md)
    — combine into one count matrix

## See also

Useful links:

- <https://ylefol.github.io/L-Atelier/HORIZON/>

- <https://github.com/Ylefol/L-Atelier/tree/master/HORIZON>

- Report bugs at <https://github.com/Ylefol/L-Atelier/issues>

## Author

**Maintainer**: Yohan Lefol <yohanlefol@proton.me>

Authors:

- Yohan Lefol <yohanlefol@proton.me>
