# Split a combined BAM into host and spike-in BAMs

After aligning to a combined host + spike-in genome with
[`HORIZON_run_bowtie2`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_bowtie2.md),
this function splits the `_mapq_filtered.bam` into two BAM files by
chromosome prefix:

- `<sample_id>_host.bam` — host reads (feed to
  [`HORIZON_process_bam`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_process_bam.md)
  and downstream)

- `<sample_id>_spikein.bam` — spike-in reads (used only for counting in
  [`HORIZON_compute_spike_in_factors`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_compute_spike_in_factors.md))

The host BAM is automatically picked up by `.latest_bam()` and all
subsequent pipeline steps without any changes to those functions.

## Usage

``` r
HORIZON_separate_spike_in(
  sample_sheet,
  sample_id,
  spikein_prefix = "spikein_",
  threads = 4L,
  remove_combined = TRUE,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame.

- sample_id:

  Character. Sample ID to process.

- spikein_prefix:

  Character. Prefix that identifies spike-in chromosomes in the BAM
  header. Must match the value used in
  [`HORIZON_build_combined_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_combined_index.md).
  Default `"spikein_"`.

- threads:

  Integer. Threads for samtools commands. Default 4.

- remove_combined:

  Logical. Remove the combined `_mapq_filtered.bam` after separation to
  save disk space. Default `TRUE`.

- force:

  Logical. If `FALSE` (default), skip separation when both `_host.bam`
  and `_spikein.bam` already exist. Set `TRUE` to re-separate and
  overwrite.

## Value

Named list with elements `host` and `spikein` giving the paths to the
two output BAMs, invisibly.
