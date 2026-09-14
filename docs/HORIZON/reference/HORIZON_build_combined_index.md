# Build a Bowtie2 index from a combined host + spike-in genome

Concatenates a host reference FASTA with a spike-in reference FASTA,
prefixing all spike-in chromosome names so that host and spike-in reads
can be cleanly separated after alignment. The combined FASTA is written
to disk and reused on subsequent calls (set `overwrite = TRUE` to
regenerate). A Bowtie2 index is then built from the combined FASTA via
[`HORIZON_build_bowtie2_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_bowtie2_index.md).

## Usage

``` r
HORIZON_build_combined_index(
  host_fasta,
  spikein_fasta,
  combined_dir,
  index_dir = combined_dir,
  index_name = "combined",
  spikein_prefix = "spikein_",
  threads = 4L,
  memory = 8000L,
  overwrite = FALSE
)
```

## Arguments

- host_fasta:

  Character. Path to the host reference FASTA (e.g. `GRCh38.fa`).

- spikein_fasta:

  Character. Path to the spike-in reference FASTA (e.g. `sacCer3.fa` for
  yeast).

- combined_dir:

  Character. Directory where the combined FASTA is written.

- index_dir:

  Character. Directory for the Bowtie2 index files. Defaults to
  `combined_dir`.

- index_name:

  Character. Basename prefix for the Bowtie2 index and the combined
  FASTA file. Default `"combined"`.

- spikein_prefix:

  Character. String prepended to every spike-in chromosome name
  (`>header` lines in the FASTA). Must not be a prefix of any host
  chromosome name. Default `"spikein_"`. Pass the same value to
  [`HORIZON_separate_spike_in`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_separate_spike_in.md).

- threads:

  Integer. Threads for `bowtie2-build`. Default 4.

- memory:

  Integer. Memory in MB for `bowtie2-build`. Default 8000.

- overwrite:

  Logical. Regenerate the combined FASTA even if it already exists.
  Default `FALSE`.

## Value

Character. Bowtie2 index basename (pass to
[`HORIZON_run_bowtie2`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_bowtie2.md)),
invisibly.

## Details

**Why a combined index?** Aligning to a single combined reference forces
reads to compete for the better mapping location. Reads that originate
from the spike-in genome map uniquely to spike-in chromosomes, while
host reads align to host chromosomes. This is more accurate than
aligning to each genome separately.
