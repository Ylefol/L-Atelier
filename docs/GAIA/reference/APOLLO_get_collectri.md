# Get CollecTRI TF-target network

Retrieves the CollecTRI transcription factor-target gene network.
CollecTRI is a comprehensive collection of TF-target interactions
compiled from multiple sources.

## Usage

``` r
APOLLO_get_collectri(
  organism = "human",
  split_complexes = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- organism:

  Character. Organism to use: "human" or "mouse". Default: "human".

- split_complexes:

  Logical. If TRUE, splits TF complexes into individual TFs. Default:
  FALSE.

- verbose:

  Logical. Print information about the network. Default: TRUE.

- ...:

  Additional arguments passed to decoupleR::get_collectri().

## Value

A tibble with columns:

- source:

  Transcription factor name

- target:

  Target gene name

- mor:

  Mode of regulation: 1 (activation) or -1 (repression)

## Details

CollecTRI integrates TF-target interactions from:

- Literature curation

- ChIP-seq experiments

- Inference from gene expression

This is the recommended default for TF activity inference as it provides
good coverage while maintaining reasonable accuracy.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get human TF-target network
network <- APOLLO_get_collectri()

# Get mouse network
network_mouse <- APOLLO_get_collectri(organism = "mouse")

} # }
```
