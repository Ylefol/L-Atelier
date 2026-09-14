# Eleuthia - Peak Loading Functions

Functions for loading and processing peak files from ATAC-seq and
ChIP-seq analyses (narrowPeak, broadPeak, BED formats). Load a
narrowPeak File

Reads a MACS2 narrowPeak file into a data.frame with standardized column
names.

## Usage

``` r
ELEUTHIA_load_narrowpeak(file_path, min_score = 0, min_qvalue = 0)
```

## Arguments

- file_path:

  Character string. Path to the narrowPeak file.

- min_score:

  Numeric. Minimum peak score to retain (default = 0, keep all).

- min_qvalue:

  Numeric. Minimum -log10(qvalue) to retain (default = 0).

## Value

A data.frame with columns: chr, start, end, name, score, strand,
signal_value, pvalue, qvalue, peak_summit, and summit_pos (absolute
position).

## Details

narrowPeak format (BED6+4):

1.  chr - Chromosome

2.  start - Start position (0-based)

3.  end - End position

4.  name - Peak name

5.  score - Peak score (0-1000)

6.  strand - Strand (usually ".")

7.  signalValue - Fold enrichment

8.  pValue - -log10(p-value)

9.  qValue - -log10(q-value)

10. peak - Offset from start to peak summit

## Examples

``` r
if (FALSE) { # \dontrun{
peaks <- ELEUTHIA_load_narrowpeak("sample_peaks.narrowPeak")
peaks <- ELEUTHIA_load_narrowpeak("sample_peaks.narrowPeak", min_qvalue = 2)

} # }
```
