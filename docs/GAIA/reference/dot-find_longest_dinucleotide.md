# Find Longest Dinucleotide Repeat

Internal helper function to find the longest dinucleotide repeat in a
sequence (e.g., ATATAT, CGCGCG).

## Usage

``` r
.find_longest_dinucleotide(seq)
```

## Arguments

- seq:

  Character string. DNA sequence (uppercase).

## Value

List with components:

- length:

  Length in base pairs of the longest dinucleotide repeat

- motif:

  The two-base motif (e.g., "AT", "CG")
