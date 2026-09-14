# Calculate Odds Ratio from Enrichment Results

Helper function to calculate odds ratio from GeneRatio and BgRatio
columns in clusterProfiler enrichment results.

## Usage

``` r
.calculate_odds_ratio(gene_ratio, bg_ratio)
```

## Arguments

- gene_ratio:

  Character vector. GeneRatio values (e.g., "10/50").

- bg_ratio:

  Character vector. BgRatio values (e.g., "100/10000").

## Value

Numeric vector of odds ratios.

## Details

Odds Ratio = (k/(n-k)) / ((M-k)/(N-M-n+k)) Where:

- k = genes in query that hit the pathway (from GeneRatio numerator)

- n = total genes in query (from GeneRatio denominator)

- M = total genes in pathway (from BgRatio numerator)

- N = total background genes (from BgRatio denominator)
