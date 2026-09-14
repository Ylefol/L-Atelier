# Filter motif enrichment results

Filter known motif results by significance thresholds and/or top N.
Works with both single (homer_motif) and batch (homer_motif_batch)
objects.

## Usage

``` r
APOLLO_filter_motifs(
  homer_result,
  q_thresh = 0.05,
  p_thresh = NULL,
  min_target_pct = NULL,
  top_n = NULL,
  verbose = TRUE
)
```

## Arguments

- homer_result:

  A homer_motif or homer_motif_batch object.

- q_thresh:

  Maximum q-value (Benjamini). Default 0.05.

- p_thresh:

  Maximum raw p-value. If set, overrides q_thresh.

- min_target_pct:

  Minimum percentage of target sequences with motif.

- top_n:

  Keep top N motifs by q-value (per peak set for batch).

- verbose:

  Logical. Print filter statistics. Default TRUE.

## Value

Same class as input, with filtered results.
