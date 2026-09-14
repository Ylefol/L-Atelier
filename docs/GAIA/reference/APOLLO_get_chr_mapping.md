# Built-in Chromosome Name Mappings

Returns a named vector mapping NCBI accession numbers to UCSC-style
chromosome names for supported genomes.

## Usage

``` r
APOLLO_get_chr_mapping(genome)
```

## Arguments

- genome:

  Character string. Genome identifier. Currently supported: "T2T" or
  "T2T-CHM13v2.0" for the Telomere-to-Telomere human genome.

## Value

A named character vector where names are NCBI accessions and values are
UCSC-style chromosome names.

## Examples

``` r
if (FALSE) { # \dontrun{
mapping <- APOLLO_get_chr_mapping("T2T")
# Use with APOLLO_get_chromosome_sizes:
chr_sizes <- APOLLO_get_chromosome_sizes("annotation.gtf", name_mapping = "T2T")

} # }
```
