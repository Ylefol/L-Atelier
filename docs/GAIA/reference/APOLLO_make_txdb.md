# Apollo - Peak Annotation & Enrichment Functions

Functions for annotating genomic peaks with gene information and
performing pathway/GO enrichment analysis. Create or Load TxDb from
GTF/GFF Annotation

Creates a TxDb object from a GTF/GFF annotation file, with automatic
caching to avoid recreating it on subsequent runs.

## Usage

``` r
APOLLO_make_txdb(
  gtf_path,
  cache_dir = NULL,
  cache_name = NULL,
  chr_mapping = NULL,
  organism = "Homo sapiens",
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- gtf_path:

  Character string. Path to GTF or GFF annotation file.

- cache_dir:

  Character string. Directory to store cached TxDb SQLite file. Default
  is "data/" relative to ZERO_DAWN root.

- cache_name:

  Character string or NULL. Name for the cached SQLite file (without
  .sqlite extension). If NULL, derives from GTF filename.

- chr_mapping:

  Character string or named vector. Chromosome name mapping to apply.
  Can be:

  - "T2T" - Use built-in T2T-CHM13 NCBI-\>UCSC mapping

  - A named vector: c("NC_060925.1" = "chr1", ...)

  - NULL - No renaming (default)

- organism:

  Character string. Organism name for TxDb metadata (default = "Homo
  sapiens").

- force:

  Logical. If TRUE, recreate TxDb even if cache exists (default =
  FALSE).

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

A TxDb object.

## Details

The function checks for a cached TxDb in cache_dir:

- If found, loads and returns it (fast)

- If not found, creates TxDb from GTF and saves to cache (slow, but
  one-time)

TxDb objects are saved using AnnotationDbi::saveDb() as SQLite
databases, which is the proper serialization method for these objects.

This is particularly useful for non-standard genomes like T2T-CHM13 that
don't have pre-built TxDb packages on Bioconductor.

## Examples

``` r
if (FALSE) { # \dontrun{
# First run: creates TxDb from GTF (slow)
txdb <- APOLLO_make_txdb("path/to/T2T_annotation.gtf", chr_mapping = "T2T")

# Subsequent runs: loads from cache (fast)
txdb <- APOLLO_make_txdb("path/to/T2T_annotation.gtf", chr_mapping = "T2T")

# Force recreation
txdb <- APOLLO_make_txdb("path/to/T2T_annotation.gtf", chr_mapping = "T2T", force = TRUE)

} # }
```
