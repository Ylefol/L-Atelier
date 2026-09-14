# Plot Gene Length Distribution by Differential Expression Status

Visualises whether up-regulated, down-regulated, and non-significant
genes differ in their total exonic length. Gene lengths are derived from
the supplied GTF file (sum of exon lengths per gene) and plotted on a
`log10` scale as smoothed density curves, one per significance class.

## Usage

``` r
AETHER_plot_gene_length_distribution(
  de_result,
  gtf_file,
  padj_thresh = 0.05,
  l2fc_thresh = 1,
  title = NULL
)
```

## Arguments

- de_result:

  Data.frame with columns `gene_id`, `padj`, and `log2FoldChange`. Also
  accepts a list with a `$results` data.frame (e.g. the direct output of
  [`ARTEMIS_differential_counts`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)).

- gtf_file:

  Character. Path to a GTF/GFF annotation file used to compute exonic
  gene lengths.

- padj_thresh:

  Numeric. Adjusted p-value threshold for significance. Default `0.05`.

- l2fc_thresh:

  Numeric. Absolute log2 fold-change threshold. Default `1`.

- title:

  Character or `NULL`. Plot title. If `NULL` a default title
  incorporating the thresholds is generated.

## Value

A `ggplot` object.

## Details

Gene length is computed internally via `.demeter_add_gene_length()`.
Genes with no matching GTF entry or a computed length of zero are
silently dropped before plotting.

Colours follow the GAIA convention: red for up-regulated, blue for
down-regulated, grey for non-significant.

## See also

[`ARTEMIS_differential_counts`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
