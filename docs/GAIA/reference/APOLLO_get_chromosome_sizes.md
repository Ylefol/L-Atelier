# Get Chromosome Sizes from GTF/GFF Annotation File

Extracts chromosome names and sizes from a GTF or GFF annotation file.
Useful for circos plots and other visualizations that need chromosome
dimensions, especially for non-standard genomes (e.g., T2T instead of
hg38). Supports caching to avoid re-parsing large annotation files.

## Usage

``` r
APOLLO_get_chromosome_sizes(
  annotation_path,
  chromosomes = NULL,
  name_mapping = NULL,
  add_chr_prefix = FALSE,
  cache_dir = NULL,
  cache_name = NULL,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- annotation_path:

  Character string. Path to GTF or GFF file. Supports gzipped files (.gz
  extension).

- chromosomes:

  Optional character vector of chromosome names to include. If NULL
  (default), includes all chromosomes found. Use this to filter to
  standard chromosomes, e.g., c(paste0("chr", 1:22), "chrX", "chrY").
  Note: filtering happens AFTER name_mapping is applied.

- name_mapping:

  Optional. Rename chromosomes from NCBI accessions to UCSC-style names.
  Can be:

  - A character string for built-in mappings: "T2T" for T2T-CHM13v2.0

  - A named character vector: c("NC_060925.1" = "chr1", ...)

- add_chr_prefix:

  Logical. If TRUE and chromosome names don't start with "chr", adds the
  prefix. Applied AFTER name_mapping (default = FALSE).

- cache_dir:

  Character string. Directory to store cached chromosome sizes. Default
  is "data/" relative to ZERO_DAWN root. Set to NULL to disable caching.

- cache_name:

  Character string or NULL. Name for the cached RDS file (without .rds
  extension). If NULL, derives from annotation filename.

- force:

  Logical. If TRUE, re-extract sizes even if cache exists (default =
  FALSE).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A data.frame with columns:

- chr:

  Chromosome name

- size:

  Chromosome size (max coordinate found in annotation)

## Details

The function determines chromosome sizes by finding the maximum
coordinate (end position) for each chromosome in the annotation file.
This is reliable for well-annotated genomes where annotations extend
near chromosome ends.

For GFF3 files with \##sequence-region pragmas, those values are used
directly as they provide exact chromosome lengths.

**Caching**: Results are cached as RDS files in cache_dir. The cache
filename includes the name_mapping used, so different mapping
configurations get separate cache files. Use force=TRUE to regenerate
the cache.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get all chromosome sizes (first run slow, subsequent runs fast)
chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf.gz")

# Filter to standard human chromosomes
standard_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf", chromosomes = standard_chrs)

# For T2T genome with built-in mapping
chr_sizes <- APOLLO_get_chromosome_sizes(
  "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf",
  name_mapping = "T2T",
  chromosomes = c(paste0("chr", 1:22), "chrX", "chrY")
)

# Force regeneration of cache
chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf", force = TRUE)

} # }
```
