# Aggregate Transcript-Level Counts to Gene-Level Counts

Collapses a transcript x sample count matrix (e.g. from
Nanopore/long-read quantification) into a gene x sample matrix by
summing transcript counts per parent gene, using the transcript-to-gene
mapping from a GTF annotation.

## Usage

``` r
POSEIDON_aggregate_transcript_to_gene(counts, gtf_file, verbose = TRUE)
```

## Arguments

- counts:

  Numeric matrix, transcripts (rows) x samples (cols). Row names must
  match the `transcript_id`(.`transcript_version`) values in `gtf_file`.

- gtf_file:

  Path to the GTF annotation used to quantify `counts`.

- verbose:

  Logical. Print a summary. Default: TRUE.

## Value

A numeric matrix, genes (rows) x samples (cols), with row names equal to
`gene_id`(.`gene_version`) as found in `gtf_file`.

## Details

Transcript and gene IDs are versioned (e.g. `"ENST00000511072.5"`) when
the GTF provides `transcript_version`/`gene_version` attributes, and
left unversioned otherwise – matching whichever convention the counts
were generated with. Transcripts in `counts` with no matching entry in
the GTF are dropped (reported when `verbose = TRUE`); this is expected
for e.g. spike-ins or annotation-version mismatches and should be
checked if the dropped fraction is large.

## Examples

``` r
if (FALSE) { # \dontrun{
gene_counts <- POSEIDON_aggregate_transcript_to_gene(
  counts   = rna_data_Nano$counts,
  gtf_file = "data/Nanopore/Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf"
)
} # }
```
