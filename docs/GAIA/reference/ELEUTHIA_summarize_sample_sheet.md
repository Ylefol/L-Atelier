# Summarize Sample Sheet Structure

Prints a summary of the sample sheet structure showing the distribution
of samples across omics types, groups, batches, etc.

## Usage

``` r
ELEUTHIA_summarize_sample_sheet(sample_sheet)
```

## Arguments

- sample_sheet:

  A validated sample sheet data.frame.

## Value

Invisibly returns a list with summary statistics.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ELEUTHIA_validate_sample_sheet("sample_sheet.csv")
ELEUTHIA_summarize_sample_sheet(result$sample_sheet)

} # }
```
