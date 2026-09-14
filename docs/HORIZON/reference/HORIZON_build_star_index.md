# Build a STAR genome index

Wraps `STAR --runMode genomeGenerate` to build a splice-aware genome
index for use with
[`HORIZON_run_star_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_star_align.md).
This is a one-time operation per reference genome/annotation/read-length
combination.

## Usage

``` r
HORIZON_build_star_index(
  genome_fasta,
  gtf,
  index_dir,
  sjdb_overhang = 100L,
  threads = 8L,
  extra_flags = character(0),
  force = FALSE
)
```

## Arguments

- genome_fasta:

  Character. Path to the reference genome FASTA.

- gtf:

  Character. Path to the gene GTF annotation (used to build the splice
  junction database via `--sjdbGTFfile`).

- index_dir:

  Character. Output directory for the STAR index (STAR requires a real
  directory, unlike Rsubread/Salmon's basename-prefix index convention).

- sjdb_overhang:

  Integer. STAR's `--sjdbOverhang`, ideally `(read_length - 1)` for the
  actual FASTQ read length of the run this index will be used for –
  **not** a value to leave at the default without checking. Default
  `100L` is STAR's own documented general-purpose value (works
  reasonably for most Illumina read lengths via a few bp of
  junction-overhang slack) but should be verified against the real read
  length, e.g. via
  `zcat <fastq> | head -400 | awk 'NR%4==2\{print length($0)\}' | sort -u`.

- threads:

  Integer. Threads for `--runThreadN`. Default 8.

- extra_flags:

  Character vector of additional STAR `genomeGenerate` flags appended
  verbatim (e.g. `c("--genomeSAindexNbases", "13")` for unusually small
  genomes). Default `character(0)` – STAR's own default of 14 is already
  appropriate for a standard human/mouse-sized genome and does not need
  overriding here.

- force:

  Logical. Rebuild even if an index already exists at `index_dir`.
  Default `FALSE`.

## Value

Character. `index_dir`, invisibly.

## Details

`STAR` must be available in the conda environment registered with
[`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md).
Install it alongside the existing `horizon_cli` tools:

      conda install -n horizon_cli -c bioconda -c conda-forge star tetranscripts
