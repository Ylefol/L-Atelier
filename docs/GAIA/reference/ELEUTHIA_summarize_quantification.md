# Summarize Quantification Statistics

Prints a summary of quantification results. Works with both BAM-based
(featureCounts) and BED-based quantification outputs.

## Usage

``` r
ELEUTHIA_summarize_quantification(quant_result)
```

## Arguments

- quant_result:

  Result from ELEUTHIA_quantify_peaks() or ELEUTHIA_quantify_bed().

## Value

Invisibly returns a summary data.frame with per-sample statistics.
