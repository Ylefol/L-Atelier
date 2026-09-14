# Get DoRothEA TF regulons

Retrieves DoRothEA transcription factor regulons with confidence level
filtering. DoRothEA provides curated TF-target interactions with
assigned confidence scores.

## Usage

``` r
APOLLO_get_dorothea(
  organism = "human",
  levels = c("A", "B", "C", "D", "E"),
  verbose = TRUE,
  ...
)
```

## Arguments

- organism:

  Character. Organism to use: "human" or "mouse". Default: "human".

- levels:

  Character vector. Confidence levels to include. Options: "A"
  (highest), "B", "C", "D", "E" (lowest). Default: c("A", "B", "C", "D",
  "E") (all levels).

- verbose:

  Logical. Print information about the network. Default: TRUE.

- ...:

  Additional arguments passed to decoupleR::get_dorothea().

## Value

A tibble with columns:

- source:

  Transcription factor name

- target:

  Target gene name

- mor:

  Mode of regulation: 1 (activation) or -1 (repression)

- confidence:

  Confidence level (A-E)

## Details

DoRothEA confidence levels:

- A: Highest confidence (multiple evidence types)

- B: High confidence (literature + ChIP-seq)

- C: Medium confidence (ChIP-seq or inference)

- D: Low confidence (inference only)

- E: Lowest confidence (predicted)

For most analyses, levels A-C are recommended. All levels are returned
by default to allow downstream filtering (e.g., at visualization).

## Examples

``` r
if (FALSE) { # \dontrun{
# Get all DoRothEA regulons
network <- APOLLO_get_dorothea()

# Get only high-confidence regulons
network_hc <- APOLLO_get_dorothea(levels = c("A", "B"))

} # }
```
