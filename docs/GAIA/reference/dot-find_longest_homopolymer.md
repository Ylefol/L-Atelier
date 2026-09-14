# Find Longest Homopolymer Run

Internal helper function to find the longest run of consecutive
identical bases in a sequence.

## Usage

``` r
.find_longest_homopolymer(seq)
```

## Arguments

- seq:

  Character string. DNA sequence (uppercase).

## Value

List with components:

- length:

  Length of longest homopolymer

- base:

  The base forming the longest run (A, T, G, C, or N)
