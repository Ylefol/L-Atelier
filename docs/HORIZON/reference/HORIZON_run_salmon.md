# Quantify transcript abundance with Salmon

Runs `salmon quant` in selective-alignment (mapping-based) mode directly
against FASTQ reads and a pre-built decoy-aware index (see
[`HORIZON_build_salmon_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_salmon_index.md)).
Unlike
[`HORIZON_run_count`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_count.md),
this does not use the genome-aligned BAM produced by
[`HORIZON_run_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_align.md)
at all – Salmon performs its own lightweight mapping against the
transcriptome, which is what allows it to resolve reads shared between
overlapping isoforms via its EM algorithm.

**Two usage modes:**

1.  **Sample-sheet mode** (default): provide `sample_sheet` and
    `sample_id`. Reads are taken from the trimmed FASTQs produced by
    [`HORIZON_run_qc_trim`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md)
    (`<output_dir>/<sample_id>/qc/<sample_id>_trimmed_R1/R2.fastq.gz`).
    Output is written to `<output_dir>/<sample_id>/salmon/`.

2.  **Direct FASTQ mode**: provide `fastq_r1` (and `fastq_r2` for
    paired-end data) directly. `sample_id` is inferred from the R1
    filename if omitted. Output goes to `output_dir` (defaults to the
    directory containing `fastq_r1`).

## Usage

``` r
HORIZON_run_salmon(
  sample_sheet = NULL,
  sample_id = NULL,
  fastq_r1 = NULL,
  fastq_r2 = NULL,
  output_dir = NULL,
  index,
  lib_type = "A",
  threads = 4L,
  extra_flags = character(0),
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame. Required in sample-sheet mode; must
  be `NULL` when `fastq_r1` is supplied.

- sample_id:

  Character. Sample ID (sample-sheet mode). Optional in direct FASTQ
  mode – inferred from the first element of `fastq_r1` if omitted.

- fastq_r1:

  Character vector or `NULL`. Path(s) to (trimmed) R1 FASTQ. Multiple
  files (e.g. multi-lane runs) are passed straight through to
  `salmon quant -1`, which treats them as one combined library – no
  pre-concatenation needed. When supplied, `sample_sheet` must be
  `NULL`. Default `NULL`.

- fastq_r2:

  Character vector or `NULL`. Path(s) to R2 FASTQ for paired-end data
  (direct FASTQ mode only), positionally matched to `fastq_r1` (same
  length, same lane order). `NULL` for single-end.

- output_dir:

  Character or `NULL`. Output directory for direct FASTQ mode. Ignored
  in sample-sheet mode. Defaults to the directory containing `fastq_r1`.

- index:

  Character. Path to a Salmon index directory (from
  [`HORIZON_build_salmon_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_salmon_index.md)).

- lib_type:

  Character. Salmon library type (`-l`). Default `"A"` (automatic
  detection – recommended unless you have a specific reason to pin it;
  see the Salmon documentation for the ISR/ISF/IU codes corresponding to
  the sample sheet's `strandedness` values).

- threads:

  Integer. Threads for `salmon quant`. Default 4.

- extra_flags:

  Character vector of additional `salmon quant` flags appended verbatim
  (e.g. `c("--seqBias", "--gcBias")`). Default `character(0)` –
  `--validateMappings` used to be the default here, but current Salmon
  versions warn that it has no effect (selective alignment, what
  `--validateMappings` used to enable, is now always the default mapping
  mode; `--sketch` is the flag to opt back into pseudoalignment).

- force:

  Logical. If `FALSE` (default), skip quantification when `quant.sf`
  already exists. Set `TRUE` to rerun.

## Value

Character. Path to the `quant.sf` file, invisibly.
