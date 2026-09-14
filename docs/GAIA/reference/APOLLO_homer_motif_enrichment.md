# Run HOMER motif enrichment on a single BED file

Wraps findMotifsGenome.pl to perform known motif enrichment and optional
de novo motif discovery on a set of genomic regions.

## Usage

``` r
APOLLO_homer_motif_enrichment(
  bed_file,
  genome,
  output_dir,
  size = 200,
  mask = TRUE,
  denovo = FALSE,
  denovo_n = 10,
  bg = NULL,
  nproc = 4,
  extra_args = NULL,
  verbose = TRUE
)
```

## Arguments

- bed_file:

  Path to BED file with peak coordinates.

- genome:

  HOMER genome string (e.g., "hg38", "mm10", "hg19").

- output_dir:

  Directory for HOMER output. Created if it doesn't exist.

- size:

  Region size around peak center for motif scanning. Default 200. Use
  "given" to use actual peak widths from the BED file.

- mask:

  Logical. Mask repeat sequences. Default TRUE (recommended for
  ATAC-seq).

- denovo:

  Logical. Enable de novo motif discovery. Default FALSE (known only).
  When FALSE, uses -nomotif flag for faster execution.

- denovo_n:

  Integer. Number of de novo motifs to find. Default 10. Maps to -S.

- bg:

  Optional path to background BED file for custom background comparison.
  E.g., all peaks as background when testing enrichment in differential
  subset.

- nproc:

  Number of processors. Default 4. Maps to -p.

- extra_args:

  Character vector of additional HOMER arguments passed directly.

- verbose:

  Logical. Print progress messages. Default TRUE.

## Value

S3 object of class "homer_motif" with components:

- known:

  Data.frame of known motif enrichment results

- denovo:

  Data.frame of de novo motifs (NULL if denovo=FALSE)

- output_dir:

  Path to HOMER output directory

- metadata:

  List with genome, parameters, run info

## Examples

``` r
if (FALSE) { # \dontrun{
# Known motifs only (fast)
result <- APOLLO_homer_motif_enrichment(
  "peaks.bed", "hg38", "results/motifs/"
)

# With custom background and de novo discovery
result <- APOLLO_homer_motif_enrichment(
  "diff_peaks.bed", "hg38", "results/motifs/",
  bg = "all_peaks.bed", denovo = TRUE
)

} # }
```
