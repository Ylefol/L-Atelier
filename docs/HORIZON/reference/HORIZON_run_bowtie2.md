# Align FASTQ reads with Bowtie2

Aligns paired-end FASTQ files against a Bowtie2 index, piping output
directly through `samtools view` (MAPQ filtering) and `samtools sort`
(coordinate sort) without writing an intermediate SAM file to disk.
Output BAM is coordinate-sorted and ready for indexing. Applies
`--no-mixed --no-discordant` by default (appropriate for all paired-end
chromatin assays).

## Usage

``` r
HORIZON_run_bowtie2(
  sample_sheet,
  sample_id,
  index,
  max_insert = 2000L,
  use_trimmed = TRUE,
  threads = 4L,
  min_mapq = 30L,
  force = FALSE,
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

  Character. Bowtie2 index basename from
  [`HORIZON_build_bowtie2_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_bowtie2_index.md).

- max_insert:

  Integer. Maximum fragment insert size (`-X`). Default `2000L`.

- use_trimmed:

  Logical. If `TRUE` (default), use Rfastp-trimmed FASTQs from
  `<output_dir>/<sample_id>/qc/`. If `FALSE`, use the raw FASTQs from
  the sample sheet (`fastq_r1` / `fastq_r2`).

- threads:

  Integer. Number of threads passed to bowtie2 (`--threads`). Default 4.

- min_mapq:

  Integer. Minimum mapping quality. Reads below this threshold are
  discarded by `samtools view -q`. Default 30.

- force:

  Logical. If `FALSE` (default), skip alignment when the output
  `_mapq_filtered.bam` already exists. Set `TRUE` to re-align and
  overwrite.

- ...:

  Additional bowtie2 flags as character strings (e.g.
  `"--very-sensitive"`). Appended verbatim to the bowtie2 command.

## Value

Character. Path to the MAPQ-filtered BAM file, invisibly. Pass to
[`HORIZON_process_bam`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_process_bam.md).

## Details

Requires `bowtie2` and `samtools` to be present in the HORIZON conda
environment registered via
[`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md).
Install via:
`conda install -n horizon_cli -c bioconda bowtie2 samtools`.

All chromatin assay processing assumes **paired-end** sequencing.
Single-end data is not supported.

**Recommended `max_insert` values by assay:**

- ATAC-seq — `2000` (fragments span up to tri-nucleosome arrays)

- ChIP-seq / CUT&RUN / CUT&TAG — `1000` (shorter sonication/tagmentation
  fragments)
