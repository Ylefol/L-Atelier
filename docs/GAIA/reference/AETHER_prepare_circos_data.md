# Prepare Data for Circos Plots

Bins fragment counts from BED files into genomic tiles for circos plot
visualization. Processes files one at a time for memory efficiency, caps
outliers using mean + 3\*SD threshold.

## Usage

``` r
AETHER_prepare_circos_data(
  sample_sheet,
  chr_sizes,
  omics,
  bin_size = 1e+06,
  aggregate = FALSE,
  cap_outliers = TRUE,
  verbose = TRUE
)
```

## Arguments

- sample_sheet:

  Data frame with sample metadata. Must contain columns: sample_id,
  bed_loc, omics, group.

- chr_sizes:

  Data frame with chromosome sizes. Must contain columns: chr, size.
  Typically from APOLLO_get_chromosome_sizes().

- omics:

  Character string or vector. Which omics types to include (e.g.,
  "ATACseq", "CHIPseq", or c("ATACseq", "CHIPseq")).

- bin_size:

  Numeric. Bin size in base pairs (default = 1e6 = 1Mb).

- aggregate:

  Logical. If TRUE, aggregate counts by group (one value per group). If
  FALSE (default), keep individual sample values.

- cap_outliers:

  Logical. If TRUE (default), cap outliers at mean + 3\*SD per
  sample/group.

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A list with:

- bins:

  Data frame with chr, start, end for each genomic bin

- signal:

  Matrix of binned signal values (bins x samples/groups)

- sample_info:

  Data frame with sample/group metadata and colors

- chr_sizes:

  The input chr_sizes (passed through for plotting)

- bin_size:

  The bin size used

- omics:

  The omics type(s) processed

## Details

The function operates directly from the sample_sheet without loading all
data into memory. BED files are processed one at a time:

1.  Load BED file

2.  Count fragments in each genomic bin

3.  Discard BED data

4.  Move to next sample

This approach minimizes memory usage for large datasets.

For aggregate = TRUE, samples within each group are summed, then
normalized by the number of samples in that group.

## Note

If more than 6 samples are processed in non-aggregated mode, a message
suggests using aggregate = TRUE for cleaner visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load chromosome sizes
chr_sizes <- APOLLO_get_chromosome_sizes(
  "annotation.gtf",
  name_mapping = "T2T",
  chromosomes = c(paste0("chr", 1:22), "chrX", "chrY")
)

# Prepare circos data for ATAC-seq
circos_data <- AETHER_prepare_circos_data(
  sample_sheet = my_sample_sheet,
  chr_sizes = chr_sizes,
  omics = "ATACseq",
  bin_size = 1e6,
  aggregate = FALSE
)

# Aggregated by group
circos_data_agg <- AETHER_prepare_circos_data(
  sample_sheet = my_sample_sheet,
  chr_sizes = chr_sizes,
  omics = c("ATACseq", "CHIPseq"),
  aggregate = TRUE
)
} # }
```
