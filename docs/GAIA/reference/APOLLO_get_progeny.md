# Get PROGENy pathway signatures

Retrieves PROGENy pathway-responsive gene signatures. PROGENy provides
genes that are responsive to pathway activity, derived from perturbation
experiments.

## Usage

``` r
APOLLO_get_progeny(organism = "human", top = 500, verbose = TRUE, ...)
```

## Arguments

- organism:

  Character. Organism to use: "human" or "mouse". Default: "human".

- top:

  Integer. Number of top responsive genes per pathway to include.
  Default: 500. Use higher values for more coverage, lower for more
  specificity.

- verbose:

  Logical. Print information about the network. Default: TRUE.

- ...:

  Additional arguments passed to decoupleR::get_progeny().

## Value

A tibble with columns:

- source:

  Pathway name

- target:

  Responsive gene name

- weight:

  Gene weight/importance for pathway activity

## Details

PROGENy pathways include:

- Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K,
  TGFb, TNFa, Trail, VEGF, WNT

The `top` parameter controls how many genes per pathway are included.
More genes = better coverage but potentially more noise. Fewer genes =
more specific but may miss relevant signals.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get PROGENy signatures (default 500 genes per pathway)
network <- APOLLO_get_progeny()

# Get more specific signatures (100 genes per pathway)
network_specific <- APOLLO_get_progeny(top = 100)

} # }
```
