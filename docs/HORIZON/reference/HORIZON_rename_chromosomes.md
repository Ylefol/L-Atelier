# Rename chromosome identifiers in BAM, BED, or BigWig files

Scans the sample output directory for files whose names end with
`file_suffix`, then renames chromosome identifiers according to
`chr_map`. Useful when aligning to a T2T or NCBI-style reference that
uses accession-style names (e.g. `NC_060925.1`) rather than conventional
UCSC-style names (e.g. `chr1`).

## Usage

``` r
HORIZON_rename_chromosomes(
  sample_sheet,
  sample_id,
  chr_map,
  file_suffix,
  overwrite = FALSE,
  rename_suffix = "_chrrenamed"
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).

- sample_id:

  Character. Sample ID to process.

- chr_map:

  Named character vector (names = old chromosome identifiers, values =
  new names), a path to a two-column tab-separated file (no header;
  column 1 = old, column 2 = new), or `"T2TCHM13"` for the built-in
  T2T-CHM13v2.0 NCBI-to-UCSC map (`NC_060925.1` \\\to\\ `chr1`, ...).

- file_suffix:

  Character. Files whose basenames end with this string (case-sensitive)
  are processed. Use a specific trailing fragment to target one file
  (e.g. `"RPKM.bw"`, `"_sorted.bam"`) or a broader ending to match
  several files of the same type. The file type (`bam`, `bed`, or `bw`)
  is derived from the extension of `file_suffix` and determines the
  renaming strategy.

- overwrite:

  Logical. If `TRUE`, the renamed file replaces the original (via a safe
  temp-file swap). Default `FALSE`: a new file is written with
  `rename_suffix` inserted before the extension.

- rename_suffix:

  Character. String appended to the base name before the extension when
  `overwrite = FALSE`. Default `"_chrrenamed"`.

## Value

Character vector of output file paths, invisibly.

## Details

**BAM**: Requires `samtools` in the conda environment registered via
[`HORIZON_set_conda_env`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_set_conda_env.md)
(only the header `@SQ SN:` fields are rewritten; the binary read data is
unchanged). The output BAM is automatically re-indexed.

**BigWig**: Requires the Bioconductor package `rtracklayer`. Install via
`BiocManager::install("rtracklayer")`.

**BED**: Handled in pure R (no extra dependencies).
