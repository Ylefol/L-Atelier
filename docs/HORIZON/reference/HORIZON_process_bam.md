# Process a BAM file: collate, fixmate, sort, mark/remove duplicates

Runs the standard samtools duplicate-removal pipeline on the
MAPQ-filtered BAM produced by
[`HORIZON_run_bowtie2`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_bowtie2.md).
Steps performed in order:

1.  **collate** — name-sort for fixmate

2.  **fixmate** — fill mate coordinates (required by markdup)

3.  **sort** — coordinate sort

4.  **markdup** — mark and remove PCR/optical duplicates

5.  **chrM removal** (optional) — remove mitochondrial reads

6.  **index** — index the final BAM

7.  **flagstat** — write alignment statistics

## Usage

``` r
HORIZON_process_bam(
  sample_sheet,
  sample_id,
  threads = 4L,
  memory_per_thread = 2000L,
  remove_chrM = TRUE,
  remove_tmp = TRUE,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame.

- sample_id:

  Character. Sample ID to process.

- threads:

  Integer. Threads for all samtools commands. Default 4.

- memory_per_thread:

  Integer. Memory in MB per thread for `samtools sort` (`-m`). Default
  2000.

- remove_chrM:

  Logical. Remove mitochondrial reads from the final BAM. Default
  `TRUE`.

- remove_tmp:

  Logical. Remove intermediate BAM files after processing. Default
  `TRUE`.

- force:

  Logical. If `FALSE` (default), skip processing when `_processed.bam`
  already exists. Set `TRUE` to reprocess.

## Value

Character. Path to the processed BAM file, invisibly. Pass to
[`HORIZON_filter_blacklist`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_filter_blacklist.md)
or
[`HORIZON_bam_to_bigwig`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_bam_to_bigwig.md).
