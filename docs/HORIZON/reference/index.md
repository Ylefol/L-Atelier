# Package index

## Sample sheets & setup

Building/validating sample sheets and registering the conda environment.

- [`HORIZON_create_sample_sheet()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_create_sample_sheet.md)
  : Create a sample sheet template
- [`HORIZON_validate_sample_sheet()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md)
  : Validate and load a sample sheet
- [`HORIZON_set_conda_env()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md)
  : Set the conda environment used by HORIZON for all external tools
- [`HORIZON_create_parse_sheet()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_create_parse_sheet.md)
  : Create a PARSE Biosciences sample sheet template
- [`HORIZON_validate_parse_sheet()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_parse_sheet.md)
  : Validate a PARSE Biosciences sample sheet

## QC & trimming

Read QC and adapter trimming.

- [`HORIZON_run_qc_trim()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_qc_trim.md)
  : Run QC and adapter trimming with fastp

## Reference indices & chromosome handling

Building aligner indices and chromosome-name/size utilities.

- [`HORIZON_build_bowtie2_index()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_bowtie2_index.md)
  : Build a Bowtie2 genome index
- [`HORIZON_build_combined_index()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_combined_index.md)
  : Build a Bowtie2 index from a combined host + spike-in genome
- [`HORIZON_build_index()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_index.md)
  : Build a Rsubread genome index
- [`HORIZON_build_salmon_index()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_salmon_index.md)
  : Build a decoy-aware Salmon transcriptome index
- [`HORIZON_build_star_index()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_star_index.md)
  : Build a STAR genome index
- [`HORIZON_get_chrom_sizes()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_get_chrom_sizes.md)
  : Get chromosome sizes
- [`HORIZON_rename_chromosomes()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_rename_chromosomes.md)
  : Rename chromosome identifiers in BAM, BED, or BigWig files

## Alignment

Genome/transcriptome alignment and BAM sorting/indexing.

- [`HORIZON_run_align()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_align.md)
  : Align reads to the genome using Rsubread
- [`HORIZON_run_bowtie2()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_bowtie2.md)
  : Align FASTQ reads with Bowtie2
- [`HORIZON_run_star_align()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_star_align.md)
  : Align reads with STAR, reporting multi-mapping alignments
- [`HORIZON_run_salmon()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_salmon.md)
  : Quantify transcript abundance with Salmon
- [`HORIZON_sort_index_bam()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_sort_index_bam.md)
  : Sort and index a BAM file

## Spike-in normalization

Host/spike-in read separation and scale-factor computation.

- [`HORIZON_separate_spike_in()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_separate_spike_in.md)
  : Split a combined BAM into host and spike-in BAMs
- [`HORIZON_compute_spike_in_factors()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_compute_spike_in_factors.md)
  : Compute per-sample scale factors from spike-in read counts

## BAM processing & filtering

Duplicate marking, blacklist filtering, downsampling, BED conversion.

- [`HORIZON_process_bam()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_process_bam.md)
  : Process a BAM file: collate, fixmate, sort, mark/remove duplicates
- [`HORIZON_filter_blacklist()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_filter_blacklist.md)
  : Remove reads overlapping a genomic blacklist
- [`HORIZON_downsample_bam()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_downsample_bam.md)
  : Downsample a BAM file by a spike-in derived scale factor
- [`HORIZON_bam_to_bed()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_bam_to_bed.md)
  : Convert a processed BAM to a fragment-level BED file
- [`HORIZON_concat_lanes()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_concat_lanes.md)
  : Concatenate multi-lane FASTQ files into one file per read direction

## Peak calling

MACS3 and SEACR peak calling.

- [`HORIZON_call_peaks()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_call_peaks.md)
  : Call peaks with MACS3
- [`HORIZON_call_peaks_seacr()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_call_peaks_seacr.md)
  : Call peaks with SEACR (pure-R reimplementation)

## Quantification

Gene, transposable element (TE), and transcript-level
quantification/aggregation.

- [`HORIZON_run_count()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_count.md)
  : Count reads per feature using featureCounts
- [`HORIZON_aggregate_counts()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_counts.md)
  : Aggregate per-sample count files into a single count matrix
- [`HORIZON_run_tecount()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_tecount.md)
  : Quantify genes + TE subfamilies for a single BAM with TEcount
- [`HORIZON_aggregate_tecount()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_tecount.md)
  : Aggregate per-sample TEcount output into a TE count matrix
- [`HORIZON_merge_gene_te_counts()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_merge_gene_te_counts.md)
  : Merge gene-level and TE-level count matrices
- [`HORIZON_aggregate_salmon()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_salmon.md)
  : Aggregate per-sample Salmon quant.sf files into TPM / count matrices

## Tracks & fragment size

BigWig generation and fragment-size distribution.

- [`HORIZON_bam_to_bigwig()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_bam_to_bigwig.md)
  : Convert a processed BAM to a BigWig using bamCoverage
- [`HORIZON_run_fragment_size()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_fragment_size.md)
  : Compute insert size distribution with bamPEFragmentSize

## Single-cell library prep (PARSE Biosciences)

split-pipe sheet creation/validation, running split-pipe, RNA velocity,
DGE filtering.

- [`HORIZON_run_splitpipe()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_splitpipe.md)
  : Run PARSE Biosciences split-pipe
- [`HORIZON_combine_splitpipe()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_combine_splitpipe.md)
  : Combine multiple PARSE split-pipe runs
- [`HORIZON_run_parse_velocity()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_parse_velocity.md)
  : Generate spliced/unspliced velocity matrices from PARSE split-pipe
  output
- [`HORIZON_parse_DGE_filter()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_parse_DGE_filter.md)
  : Filter velocity AnnData using PARSE DGE output and add metadata

## Logging

Run logging utilities shared across pipeline steps.

- [`HORIZON_generate_log_filename()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_generate_log_filename.md)
  : Generate Timestamp-Based Log Filename
- [`HORIZON_get_log_file()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_get_log_file.md)
  : Get Current Log File Path
- [`HORIZON_is_logging()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_is_logging.md)
  : Check if Logging is Active
- [`HORIZON_start_log()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_start_log.md)
  : Start Console Logging
- [`HORIZON_stop_log()`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_stop_log.md)
  : Stop Console Logging
- [`.log_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/dot-log_env.md)
  : HORIZON - Logging Functions
