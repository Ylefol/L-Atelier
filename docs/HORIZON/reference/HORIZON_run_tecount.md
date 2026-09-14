# Quantify genes + TE subfamilies for a single BAM with TEcount

Wraps the `TEcount` CLI (part of the bioconda `tetranscripts` package)
to quantify both genes and transposable element (TE) subfamilies from a
single BAM, using an EM algorithm to redistribute multi-mapping reads
across candidate TE loci rather than discarding them.

## Usage

``` r
HORIZON_run_tecount(
  sample_sheet = NULL,
  sample_id = NULL,
  bam_file = NULL,
  output_dir = NULL,
  strandedness = NULL,
  gene_gtf,
  te_gtf,
  mode = "multi",
  force = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).
  Required in sample-sheet mode; must be `NULL` when `bam_file` is
  supplied.

- sample_id:

  Character. Sample ID to process (sample-sheet mode). Optional in
  direct BAM mode – inferred from the filename if omitted.

- bam_file:

  Character or `NULL`. Path to a STAR-aligned, coordinate-sorted BAM
  with multi-mapping alignments reported (see
  [`HORIZON_run_star_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_star_align.md)).
  When supplied, `sample_sheet` must be `NULL`. Default `NULL`.

- output_dir:

  Character or `NULL`. Output directory for direct BAM mode. Ignored in
  sample-sheet mode. Defaults to the directory containing `bam_file`.

- strandedness:

  Character or `NULL`. One of `"unstranded"`, `"forward"`, `"reverse"`.
  Required in direct BAM mode; read from the sample sheet in
  sample-sheet mode. Mapped internally to `TEcount`'s own `--stranded`
  vocabulary (`"no"`/`"forward"`/`"reverse"`).

- gene_gtf:

  Character. Path to the gene GTF annotation – should be the same GTF
  already used for the project's canonical gene-level counting, for
  consistency (even though its gene-count output is discarded here).

- te_gtf:

  Character. Path to the RepeatMasker-derived TE GTF
  (`gene_id`/`family_id`/`class_id` composite format).

- mode:

  Character. `TEcount`'s `--mode`: `"multi"` (default) applies EM-based
  redistribution of multi-mapping reads across TE loci – the actual
  point of using TEcount over featureCounts. `"uniq"` counts
  unique-mapping reads only, useful as a robustness comparison against
  the `"multi"` result but not the primary mode.

- force:

  Logical. If `FALSE` (default), skip quantification when the output
  `.cntTable` already exists. Set `TRUE` to rerun.

- verbose:

  Logical. If `TRUE` (default), prints the full `TEcount` command before
  executing it.

- ...:

  Additional `TEcount` flags passed verbatim.

## Value

Character. Path to the `.cntTable` file, invisibly.

## Details

`TEcount` must be available in the conda environment registered with
[`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md).
Install it alongside `STAR` (see
[`HORIZON_build_star_index`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_build_star_index.md))
and the existing `horizon_cli` tools:

      conda install -n horizon_cli -c bioconda -c conda-forge star tetranscripts

**Only the TE rows of `TEcount`'s output are meant to be used
downstream** (see
[`HORIZON_aggregate_tecount`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_tecount.md))
– its gene-count half comes from a different aligner/counting method
than HORIZON's canonical gene-level pipeline
([`HORIZON_run_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_align.md) +
[`HORIZON_run_count`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_count.md))
and would be inconsistent with already-computed differential expression
results if mixed in. Canonical gene counts should keep coming from the
existing Rsubread/featureCounts path; merge the two feature sets with
[`HORIZON_merge_gene_te_counts`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_merge_gene_te_counts.md).

**Two usage modes:**

1.  **Sample-sheet mode** (default): provide `sample_sheet` and
    `sample_id`. `strandedness` is read from the sample sheet row, and
    the BAM is located at the standard HORIZON path
    `<output_dir>/<sample_id>/aligned_star/<sample_id>_sorted.bam` (i.e.
    [`HORIZON_run_star_align`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_run_star_align.md)'s
    output – **not** the Rsubread BAM). Output is written to
    `<output_dir>/<sample_id>/tecount/`.

2.  **Direct BAM mode**: provide `bam_file`, plus `strandedness`
    explicitly, since it cannot be reliably inferred without a sample
    sheet row. `sample_id` is inferred from the filename if omitted.
    Output goes to `output_dir` (defaults to the directory containing
    the BAM).
