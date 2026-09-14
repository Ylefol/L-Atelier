# Export Sequence Composition Results

Saves the data.frame returned by
[`APOLLO_sequence_composition()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_sequence_composition.md)
to CSV, with an optional separate file for the raw DNA sequences and a
plain-text summary of composition statistics.

## Usage

``` r
ELEUTHIA_export_sequence_composition(
  composition,
  output_dir,
  prefix = "composition",
  save_sequences = FALSE,
  save_rds = FALSE,
  verbose = TRUE
)
```

## Arguments

- composition:

  A data.frame from
  [`APOLLO_sequence_composition()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_sequence_composition.md),
  or a named list of such data.frames (one per peak set). When a list is
  provided, each element is exported to its own subdirectory under
  `output_dir` named after the list element. Expected columns:
  `peak_id`, `width`, `gc_percent`, `at_percent`, `n_Nbases`. Optional:
  `longest_homopolymer`, `homopolymer_base`, `longest_dinucleotide`,
  `dinucleotide_motif`, `sequence`.

- output_dir:

  Character. Directory to save results. Created if needed.

- prefix:

  Character. Prefix for output filenames. Default: "composition".

- save_sequences:

  Logical. If `TRUE` and a `sequence` column exists, save sequences to a
  separate CSV and drop the column from the main export. Default: FALSE
  (sequences omitted from main CSV but not saved separately).

- save_rds:

  Logical. Save the full data.frame as RDS. Default: FALSE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

Invisible character vector of file paths created.

## Details

Exports the following files:

- **Composition CSV** — all numeric columns without the sequence string
  (`{prefix}_composition.csv`)

- **Sequences CSV** — two-column table (`peak_id`, `sequence`) if
  `save_sequences = TRUE` and a `sequence` column is present
  (`{prefix}_sequences.csv`)

- **Summary TXT** — descriptive statistics for GC%, AT%, region width,
  and repeat metrics (`{prefix}_summary.txt`)

- **RDS** — full data.frame if `save_rds = TRUE`

## Examples

``` r
if (FALSE) { # \dontrun{
comp <- APOLLO_sequence_composition(peaks, "genome.fa")
ELEUTHIA_export_sequence_composition(comp, "results/composition/")

# Named list — each set gets its own subdirectory
comp_list <- APOLLO_sequence_composition(list(WT = wt, KO = ko), "genome.fa")
ELEUTHIA_export_sequence_composition(comp_list, "results/composition/")
# -> results/composition/WT/composition_composition.csv
# -> results/composition/KO/composition_composition.csv

# Save sequences separately
ELEUTHIA_export_sequence_composition(comp, "results/composition/",
                                      save_sequences = TRUE)
} # }
```
