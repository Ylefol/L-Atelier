# Call peaks with MACS3

Wraps the `macs3 callpeak` command via the user-managed conda
environment registered with
[`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md).

## Usage

``` r
HORIZON_call_peaks(
  treatment_bed,
  control_bed = NULL,
  output_dir,
  sample_name,
  genome_size = "hs",
  q_value = 0.05,
  broad = FALSE,
  nomodel = FALSE,
  shift = NULL,
  extsize = NULL,
  extra_flags = character(0),
  force = FALSE
)
```

## Arguments

- treatment_bed:

  Character. Path to the treatment fragment BED file (from
  [`HORIZON_bam_to_bed`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_bam_to_bed.md)).

- control_bed:

  Character or `NULL`. Path to a control (input/IgG) BED file.
  Recommended for ChIP-seq; optional for CUT&TAG/CUT&RUN; not used for
  ATAC-seq.

- output_dir:

  Character. Directory for MACS3 output files.

- sample_name:

  Character. Prefix applied to all output filenames.

- genome_size:

  Character or numeric. Effective genome size for MACS3 (`-g`).
  Shortcuts: `"hs"`, `"mm"`, `"ce"`, `"dm"`. Pass a numeric value for
  other genomes.

- q_value:

  Numeric. FDR threshold (`-q`). Default 0.05.

- broad:

  Logical. Call broad peaks (`--broad`). Recommended for histone mark
  ChIP-seq. Default `FALSE`.

- nomodel:

  Logical. Skip the MACS3 shift-model building step (`--nomodel`).
  Required when `shift` and `extsize` are set manually. Default `FALSE`.

- shift:

  Integer or `NULL`. Read shift in bp (`--shift`). Only used when
  `nomodel=TRUE`. Use `-100` for ATAC-seq. Default `NULL` (not passed).

- extsize:

  Integer or `NULL`. Fragment extension size in bp (`--extsize`). Only
  used when `nomodel=TRUE`. Use `200` for ATAC-seq. Default `NULL` (not
  passed).

- extra_flags:

  Character vector of additional MACS3 flags to append verbatim (e.g.
  `c("--keep-dup", "all")` for CUT&TAG/CUT&RUN). Default `character(0)`.

- force:

  Logical. If `FALSE` (default), skip peak calling when the output peak
  file already exists. Set `TRUE` to rerun.

## Value

Character. Path to the narrowPeak (or broadPeak) file, invisibly.

## Details

**Recommended settings by assay:**

- ATAC-seq — `nomodel=TRUE, shift=-100, extsize=200` (shift model is not
  appropriate for open-chromatin data)

- ChIP-seq (transcription factor) — defaults (MACS3 builds a shift model
  from the data; provide `control_bed`)

- ChIP-seq (histone mark) — `broad=TRUE`

- CUT&TAG / CUT&RUN — `nomodel=TRUE`, optionally
  `extra_flags="--keep-dup all"` if duplicates were not removed upstream
  (CUT&TAG can generate high local duplicate rates from genuine signal
  rather than PCR artefact)
