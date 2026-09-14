# Pre-load BigWig Coverage for Fast Multi-Region Extraction

Reads all BigWig files in a sample sheet once and stores their coverage
as RLE objects. The result can be passed to
[`ARTEMIS_prepare_bw_comparison`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_prepare_bw_comparison.md)
via the `coverage` argument to skip redundant file I/O when extracting
signal over many region sets from the same samples.

## Usage

``` r
ARTEMIS_load_bw_coverage(
  sample_sheet,
  sample_col = "sample_id",
  bw_col = "bw_path"
)
```

## Arguments

- sample_sheet:

  data.frame. Must include columns `sample_col` and `bw_col`. Duplicate
  sample IDs are deduplicated automatically.

- sample_col:

  Character. Column holding sample IDs. Default `"sample_id"`.

- bw_col:

  Character. Column holding BigWig file paths. Default `"bw_path"`.

## Value

An `artemis_bw_coverage` object (a list) with:

- `coverage`:

  Named list of RleList objects, one per sample.

- `bw_chrnames`:

  Character vector of chromosome names from the first BigWig, used for
  chromosome name harmonization in `ARTEMIS_prepare_bw_comparison`.
