# Plot Distribution of Fragment Counts in Candidate Regions

Visualizes the distribution of fragment counts across candidate regions
to help choose an appropriate selection threshold.

## Usage

``` r
ELEUTHIA_plot_region_distribution(
  candidate_regions,
  show_quantiles = c(0.5, 0.75, 0.9, 0.95, 0.99),
  log_scale = TRUE,
  title = "Distribution of Fragment Counts"
)
```

## Arguments

- candidate_regions:

  Data.frame from ELEUTHIA_merge_fragments() with n_fragments column.

- show_quantiles:

  Numeric vector of quantiles to mark on the plot (default = c(0.50,
  0.75, 0.90, 0.95, 0.99)).

- log_scale:

  Logical. Use log10 scale for x-axis (default = TRUE). Useful when
  fragment counts span several orders of magnitude.

- title:

  Character. Plot title.

## Value

A ggplot object. Also prints quantile summary to console.

## Details

The plot shows:

- Histogram of fragment counts across all candidate regions

- Vertical lines marking key quantiles

- A table showing how many regions would be kept at each quantile

Use this to decide on a threshold before calling
ELEUTHIA_select_regions().

## Examples

``` r
if (FALSE) { # \dontrun{
candidates <- ELEUTHIA_merge_fragments(bed_files)
p <- ELEUTHIA_plot_region_distribution(candidates)
print(p)

# Then select based on what you see
regions <- ELEUTHIA_select_regions(candidates, quantile_threshold = 0.90)

} # }
```
