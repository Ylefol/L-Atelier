# Pairwise ranked gene list from marker scores (for pairwise GSEA)

Extracts a per-gene ranking statistic for one specific pair of clusters
from a
[`KERAUNOS_score_markers`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_score_markers.md)
result run with `full_stats = TRUE`. Unlike `$top`/`$markers`' default
summary columns (`median.AUC`, etc., aggregated across *every* other
cluster), this pulls the effect size against exactly `cluster2`, giving
a proper pairwise ranking for two clusters without requiring pseudobulk
replicate samples — effect sizes are computed at single-cell resolution
by
[`scran::scoreMarkers()`](https://rdrr.io/pkg/scran/man/scoreMarkers.html).
Feed the result directly into
[`KERAUNOS_gsea`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md)'s
`ranked_genes` argument.

## Usage

``` r
KERAUNOS_rank_pairwise_markers(
  markers_result,
  cluster1,
  cluster2,
  stat = c("AUC", "logFC.cohen", "logFC.detected"),
  verbose = TRUE
)
```

## Arguments

- markers_result:

  A `keraunos_markers` object from
  `KERAUNOS_score_markers(..., full_stats = TRUE)`.

- cluster1:

  Character. The focal cluster (positive ranking values mean higher in
  this cluster).

- cluster2:

  Character. The comparison cluster (negative ranking values, for
  `stat = "AUC"`, mean higher in this cluster).

- stat:

  Character. Which pairwise effect size to use: `"AUC"` (default; scaled
  cell-detection-rank statistic, shifted by `-0.5` here so 0 = no
  difference, matching `scran`'s convention that 0.5 = random/no
  difference — the shift makes the vector suitable as a signed GSEA
  ranking statistic), `"logFC.cohen"` (Cohen's d, already signed and
  centred at 0), or `"logFC.detected"` (log-fold change in detection
  rate, already signed).

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

A named numeric vector (names = gene symbols), NA-removed and sorted
descending — ready to pass to `KERAUNOS_gsea(ranked_genes = )`.
