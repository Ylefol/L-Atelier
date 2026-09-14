# Eleuthia - Quantification Functions

Functions for quantifying reads against genomic features. Includes both
BAM-based (featureCounts) and BED-based quantification. Quantify Reads
in Peak Regions (BAM-based)

Counts reads from BAM files in consensus peak regions using
Rsubread::featureCounts().

## Usage

``` r
ELEUTHIA_quantify_peaks(
  sample_sheet,
  consensus_peaks,
  omics,
  paired_end = TRUE,
  count_fragments = TRUE,
  min_mapq = 10,
  nthreads = 4,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame with bam_loc column. Note: The
  standard sample sheet uses bed_loc; you may need to add bam_loc
  manually if using this function.

- consensus_peaks:

  A consensus peak data.frame from ELEUTHIA_create_consensus_peaks().

- omics:

  Character string. Omics type to quantify ("ATACseq" or "CHIPseq").

- paired_end:

  Logical. Are the reads paired-end? (default = TRUE for ATAC-seq).

- count_fragments:

  Logical. If paired-end, count fragments instead of reads (default =
  TRUE).

- min_mapq:

  Integer. Minimum mapping quality to count a read (default = 10).

- nthreads:

  Integer. Number of threads for parallel processing (default = 4).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A list containing:

- counts:

  Matrix of read counts (peaks x samples)

- annotation:

  Data.frame with peak annotations

- targets:

  Data.frame with sample information

- stat:

  Data.frame with counting statistics

## Details

This function wraps Rsubread::featureCounts() for ease of use with the
GAIA pipeline. It:

1.  Extracts BAM file paths from the sample sheet

2.  Converts consensus peaks to SAF format

3.  Runs featureCounts with appropriate parameters for ATAC/ChIP-seq

4.  Returns a structured list with counts and metadata

For ATAC-seq data with paired-end reads, the function counts fragments
(read pairs) by default rather than individual reads.

## Note

This function requires BAM files and a sample sheet with a bam_loc
column. For most use cases, consider using ELEUTHIA_quantify_bed()
instead, which quantifies directly from BED fragment files and is more
lightweight.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load and create consensus peaks
peak_list <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
consensus <- ELEUTHIA_create_consensus_peaks(peak_list, min_overlap = 2)

# Quantify reads
counts <- ELEUTHIA_quantify_peaks(sample_sheet, consensus, "ATACseq")

# Access count matrix
count_matrix <- counts$counts

} # }
```
