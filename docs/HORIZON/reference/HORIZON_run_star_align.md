# Align reads with STAR, reporting multi-mapping alignments

Wraps `STAR --runMode alignReads` with a relaxed multi-mapping filter
suitable for downstream transposable-element (TE) quantification via
[`HORIZON_run_tecount`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_tecount.md).
This is a separate alignment path from
[`HORIZON_run_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_align.md)
(Rsubread) – output is written to a distinct `aligned_star/` subfolder
so both alignment methods can coexist per sample.

## Usage

``` r
HORIZON_run_star_align(
  sample_sheet,
  sample_id,
  index,
  use_trimmed = TRUE,
  threads = 8L,
  outfilter_multimap_nmax = 100L,
  win_anchor_multimap_nmax = 100L,
  tmp_dir = NULL,
  force = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).

- sample_id:

  Character. Sample ID to process.

- index:

  Character. Path to a STAR genome index directory (from
  [`HORIZON_build_star_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_star_index.md)).

- use_trimmed:

  Logical. Use trimmed FASTQs from the `qc/` subfolder (output of
  [`HORIZON_run_qc_trim`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md)).
  Default TRUE. Set FALSE to align the raw FASTQs listed in the sample
  sheet.

- threads:

  Integer. Number of alignment threads (`--runThreadN`). Default 8.

- outfilter_multimap_nmax:

  Integer. STAR's `--outFilterMultimapNmax` – maximum number of loci a
  read is allowed to map to before being discarded entirely. Default
  100, the value cited in the TEtranscripts tutorial/paper as standard
  practice for TE quantification; worth a second look if a TE family in
  this genome build has a pathologically higher copy number than that.

- win_anchor_multimap_nmax:

  Integer. STAR's `--winAnchorMultimapNmax` – maximum number of loci an
  "anchor" seed is allowed to map to during seed search. Default 100,
  paired with `outfilter_multimap_nmax` per the same
  TEtranscripts-tutorial recipe.

- tmp_dir:

  Character or `NULL`. Parent directory for STAR's own working temp
  directory (`--outTmpDir`), where it creates a FIFO to feed
  `--readFilesCommand zcat` output through. STAR's own default
  (`<outFileNamePrefix>_STARtmp`, alongside the output BAM) breaks when
  the sample sheet's `output_dir` points at an exFAT/NTFS-mounted drive,
  since those filesystems can't hold POSIX FIFOs – the alignment output
  BAM itself has no such restriction. Default `NULL` uses
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) (the R session's
  own temp directory), which is always on the local filesystem
  regardless of where `output_dir` points. Pass an explicit path if
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) doesn't have
  enough free space for STAR's temp files.

- force:

  Logical. If `FALSE` (default), skip alignment when the sorted output
  BAM already exists. Set `TRUE` to realign.

- verbose:

  Logical. If `TRUE` (default), prints the full STAR command before
  executing it.

- ...:

  Additional STAR flags passed verbatim (e.g.
  `c("--outSAMattributes", "NH", "HI", "AS", "nM")` to add fields beyond
  STAR's own `Standard` set). Note that `Standard` already includes the
  `NH` (number-of-hits) tag
  [`HORIZON_run_tecount()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_tecount.md)
  needs for multi-mapper redistribution – do not strip it via a custom
  `--outSAMattributes` list here.

## Value

Character. Path to the sorted, indexed BAM file, invisibly.

## Details

By default, uses the trimmed FASTQ files produced by
[`HORIZON_run_qc_trim`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md).
Raw FASTQs from the sample sheet can be used instead by setting
`use_trimmed = FALSE`.

**Why the multi-mapping filter matters here:** STAR's own default
(`--outFilterMultimapNmax 10`) silently discards any read mapping to
more than 10 genomic loci. For gene-level RNA-seq this is a reasonable
behavior, but it systematically discards reads from exactly the
high-copy-number TE families this workflow is meant to quantify.
Relaxing this filter (default here: 100) keeps those reads in the
alignment so
[`HORIZON_run_tecount`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_tecount.md)
can redistribute them across TE loci via its EM algorithm rather than
losing them entirely.
