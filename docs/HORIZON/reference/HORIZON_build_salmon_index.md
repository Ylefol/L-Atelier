# Build a decoy-aware Salmon transcriptome index

Extracts transcript sequences from a genome FASTA + GTF annotation (via
[`extractTranscriptSeqs`](https://rdrr.io/pkg/GenomicFeatures/man/extractTranscriptSeqs.html)),
builds a decoy-aware Salmon index following Salmon's documented best
practice (the whole genome is included as "decoy" sequence, which
prevents reads originating from unspliced/intergenic genomic sequence
from being spuriously assigned to a similar transcript), and writes the
index to `index_dir`.

This is a one-time operation per reference genome/annotation pair. The
`salmon` binary must be available in the conda environment registered
with
[`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md).

## Usage

``` r
HORIZON_build_salmon_index(
  gtf,
  genome_fasta,
  index_dir,
  kmer_length = 31L,
  keep_duplicates = FALSE,
  threads = 4L,
  force = FALSE
)
```

## Arguments

- gtf:

  Character. Path to the GTF annotation used to define transcripts.

- genome_fasta:

  Character. Path to the reference genome FASTA. Indexed in place
  (`.fai`) if not already indexed.

- index_dir:

  Character. Output directory for the Salmon index.

- kmer_length:

  Integer. Salmon `-k` (minimum acceptable match length for
  quasi-mapping). Default 31, appropriate for reads \>= 75bp; reduce
  (e.g. to 23) for shorter reads. See the Salmon documentation for
  guidance.

- keep_duplicates:

  Logical. Passes `--keepDuplicates` to `salmon index` when `TRUE`.
  Without it, Salmon collapses transcripts with identical sequence into
  a single representative ID at index time, silently dropping the others
  from quantification – this matters for isoform-level analyses (e.g.
  downstream isoform-switch testing) where every transcript ID needs its
  own entry, and is generally recommended by Salmon/tximport's own
  tutorials regardless of downstream use. Set `FALSE` (default) to match
  Salmon's own out-of-the-box default.

- threads:

  Integer. Threads for `salmon index`. Default 4.

- force:

  Logical. Rebuild even if an index already exists at `index_dir`.
  Default `FALSE`.

## Value

Character. `index_dir`, invisibly.

## Details

Intermediate files (the extracted transcript FASTA, the
genome+transcript "gentrome" FASTA used as Salmon's direct input, and
the decoy sequence name list) are written to `<index_dir>_prep/`
alongside the index, for provenance/debugging. They are not required
after the index is built and may be deleted to reclaim disk space.

Requires: `GenomicFeatures`, `Biostrings`, `Rsamtools`.
