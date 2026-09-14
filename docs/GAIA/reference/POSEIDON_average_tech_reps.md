# Average Technical Replicates

Aggregates technical replicates by averaging counts within each
biological replicate. Produces one sample per unique combination of
grouping variables (excluding tech_rep).

## Usage

``` r
POSEIDON_average_tech_reps(
  quant_result,
  bio_rep_col = "bio_rep",
  tech_rep_col = "tech_rep",
  group_by = NULL,
  method = "mean",
  verbose = TRUE
)
```

## Arguments

- quant_result:

  A list containing:

  - counts: Matrix of counts (features x samples)

  - targets: Data.frame with sample metadata

- bio_rep_col:

  Character string. Column identifying biological replicates (default =
  "bio_rep").

- tech_rep_col:

  Character string. Column identifying technical replicates to average
  over (default = "tech_rep").

- group_by:

  Character vector. Additional columns to group by when averaging.
  Samples are averaged within each unique combination of bio_rep_col +
  group_by columns. Default includes "group" and "batch" if present in
  targets.

- method:

  Character string. Aggregation method: "mean" (default) or "median".

- verbose:

  Logical. Print summary (default = TRUE).

## Value

A quant_result with:

- counts: Averaged count matrix (features x biological samples)

- targets: Updated metadata with one row per biological sample

- n_tech_reps: Number of technical replicates averaged per sample

## Details

This function averages technical replicates to produce one observation
per biological sample. This is appropriate when:

- Preparing data for integration methods requiring matched samples

- Reducing technical noise while preserving biological variation

- Creating a cleaner dataset for correlation-based analyses

The function automatically detects which columns to use for grouping:

- Always groups by bio_rep_col

- Includes "group" if present (preserves experimental groups)

- Includes "batch" if present (keeps batches separate)

- Additional columns can be specified via group_by parameter

Can be used before or after filtering by group - the function respects
whatever samples are present in the input.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic averaging
averaged <- POSEIDON_average_tech_reps(atac_counts)

# Average after filtering for WT
wt_data <- POSEIDON_filter_by_group(atac_counts, groups = "WT")
wt_averaged <- POSEIDON_average_tech_reps(wt_data)

# Average before filtering (groups will be preserved)
averaged <- POSEIDON_average_tech_reps(atac_counts)
wt_averaged <- POSEIDON_filter_by_group(averaged, groups = "WT")

# Custom grouping
averaged <- POSEIDON_average_tech_reps(data, group_by = c("group", "batch", "treatment"))

} # }
```
