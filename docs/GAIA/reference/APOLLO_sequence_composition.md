# Calculate Sequence Composition for Genomic Regions

Extracts DNA sequences for genomic regions and calculates composition
metrics including GC/AT content and repeat patterns.

## Usage

``` r
APOLLO_sequence_composition(
  regions,
  fasta_path,
  chr_mapping = NULL,
  extend = 0,
  min_width = 100,
  include_repeats = TRUE,
  include_sequence = TRUE,
  verbose = TRUE
)
```

## Arguments

- regions:

  A data.frame with columns: chr, start, end. Additional columns are
  preserved. Coordinates should be 0-based (BED format).

- fasta_path:

  Path to a FASTA file. Must have an accompanying .fai index (create
  with `samtools faidx`).

- chr_mapping:

  Either:

  - NULL (default) - no chromosome name translation

  - A genome name (e.g., "T2T") - uses built-in mapping

  - A named character vector mapping region chr names to FASTA chr names

- extend:

  Integer. Extend regions by this many bp on each side before extracting
  sequence (default = 0). Useful for standardizing region sizes.

- min_width:

  Integer. Minimum region width (after extension) to calculate
  composition (default = 100). Regions below this return NA values.

- include_repeats:

  Logical. Calculate homopolymer and dinucleotide repeat metrics
  (default = TRUE).

- include_sequence:

  Logical. Include the extracted sequence in output (default = TRUE).
  Set to FALSE to reduce memory usage for large datasets.

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

When `regions` is a data.frame: a data.frame with one row per input
region. When `regions` is a named list of data.frames: a named list of
such data.frames, one per input set.

Each data.frame contains:

- peak_id:

  Region identifier (from input or generated)

- width:

  Final region width after extension

- gc_percent:

  Percentage of G+C bases

- at_percent:

  Percentage of A+T bases

- n_Nbases:

  Number of N (ambiguous) bases

- longest_homopolymer:

  Length of longest single-base repeat (if include_repeats)

- homopolymer_base:

  Base forming the longest homopolymer (if include_repeats)

- longest_dinucleotide:

  Length of longest dinucleotide repeat in bp (if include_repeats)

- dinucleotide_motif:

  Motif of longest dinucleotide repeat (if include_repeats)

- sequence:

  The actual nucleotide sequence extracted (if include_sequence)

## Details

**Chromosome mapping**: When using genomes like T2T-CHM13 where your
regions use UCSC-style names (chr1, chr2) but the FASTA uses NCBI
accessions (NC_060925.1, NC_060926.1), use chr_mapping = "T2T" to
automatically translate. The mapping is applied internally - your output
will retain the original region chromosome names.

**Extension behavior**: Consistent with ELEUTHIA_expand_regions(), the
extend parameter adds the specified bp to EACH side. So extend = 50 adds
100bp total. Start coordinates are clamped at 0.

**Minimum width threshold**: Regions smaller than min_width (after
extension) return NA for composition values. This avoids unreliable
percentages from very short sequences (e.g., 10bp regions where each
base is 10%).

**Repeat detection**:

- Homopolymers: Consecutive identical bases (e.g., AAAA, TTTTTT)

- Dinucleotide repeats: Alternating two-base patterns (e.g., ATATAT,
  CGCGCG)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with genomic regions
comp <- APOLLO_sequence_composition(regions, "reference.fa")

# Named list of peak sets (returns named list of data.frames)
comp_list <- APOLLO_sequence_composition(
  list(WT = wt_peaks, KO = ko_peaks),
  "T2T.fna", chr_mapping = "T2T"
)

# With T2T genome (regions have chr1, FASTA has NC_060925.1)
comp <- APOLLO_sequence_composition(regions, "T2T.fna", chr_mapping = "T2T")

# With extension to standardize region size
comp <- APOLLO_sequence_composition(regions, "reference.fa", extend = 100)

# Without sequence column (saves memory)
comp <- APOLLO_sequence_composition(regions, "ref.fa", include_sequence = FALSE)

# Merge back to original data
regions_with_comp <- cbind(regions, comp[, c("gc_percent", "at_percent")])

} # }
```
