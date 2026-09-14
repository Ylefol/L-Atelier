# Parse REDItools2 tab-delimited output into a data.frame

REDItools2 emits a header line beginning with "Region" and then one row
per position. The BaseCount\[ACGT\] column (e.g. `"[0, 5, 100, 3]"`) is
split into four integer columns.

## Usage

``` r
.cyan_parse_reditools_output(path)
```
