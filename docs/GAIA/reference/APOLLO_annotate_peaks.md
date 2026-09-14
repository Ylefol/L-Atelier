# Annotate Peaks with Genomic Features and Nearest Genes

Annotates genomic peaks with their genomic context (promoter, exon,
intron, intergenic, etc.) and nearest gene using ChIPseeker.

## Usage

``` r
APOLLO_annotate_peaks(
  regions,
  txdb,
  tss_region = c(-3000, 3000),
  level = "transcript",
  verbose = TRUE
)
```

## Arguments

- regions:

  Data.frame with at minimum: chr, start, end. Can also have peak_id or
  name column for identification.

- txdb:

  A TxDb object (from APOLLO_make_txdb or Bioconductor).

- tss_region:

  Numeric vector of length 2. Region around TSS to define as promoter,
  e.g., c(-3000, 3000) means 3kb upstream to 3kb downstream.

- level:

  Character. Annotation level - "transcript" or "gene" (default =
  "transcript").

- verbose:

  Logical. Print summary (default = TRUE).

## Value

A data.frame with original region info plus annotation columns:

- annotation:

  Genomic feature (Promoter, Exon, Intron, etc.)

- annotation_simple:

  Simplified category

- gene_id:

  Nearest gene ID

- gene_name:

  Gene symbol (if available in TxDb)

- distance_to_tss:

  Distance to nearest TSS

- transcript_id:

  Associated transcript ID

## Details

Uses ChIPseeker::annotatePeak() for annotation. The TSS region defines
what counts as a "promoter" - regions within this window of any TSS are
labeled as promoter regions.

## Examples

``` r
if (FALSE) { # \dontrun{
txdb <- APOLLO_make_txdb("annotation.gtf")
annotated <- APOLLO_annotate_peaks(significant_peaks, txdb)

# See distribution of genomic features
table(annotated$annotation_simple)

# Get genes for pathway analysis
genes <- unique(annotated$gene_id)

} # }
```
