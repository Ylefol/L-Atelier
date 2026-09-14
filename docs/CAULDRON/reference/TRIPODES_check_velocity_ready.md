# Check whether an SCE object is ready for RNA velocity analysis

Validates that a SingleCellExperiment contains the spliced and unspliced
count assays required by `TRIPODES_run_velocity`. Detects alternate or
incompatible naming conventions and reports actionable guidance. Also
checks data type (must be raw integer counts) and reports summary
splicing metrics.

## Usage

``` r
TRIPODES_check_velocity_ready(
  sce,
  spliced_assay = "spliced",
  unspliced_assay = "unspliced",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- spliced_assay:

  Name of the spliced counts assay (default `"spliced"`).

- unspliced_assay:

  Name of the unspliced counts assay (default `"unspliced"`).

- verbose:

  Logical; print diagnostic report (default `TRUE`).

## Value

Invisibly returns a named list:

- ready:

  Logical — `TRUE` only if both assays are present, have matching
  dimensions, contain integer-like counts, and no threshold warnings
  were triggered.

- spliced_assay, unspliced_assay:

  Confirmed assay names, or `NA` if not found.

- n_cells, n_genes:

  Dimensions of the SCE.

- median_spliced, median_unspliced:

  Median total UMI per cell for each modality (`NA` if assays absent).

- splicing_ratio:

  Global splicing ratio: sum(spliced) / (sum(spliced) + sum(unspliced));
  `NA` if assays absent.

- issues:

  Character vector of all warnings/errors found (empty if
  `ready = TRUE`).

## Details

This function is non-destructive — it does **not** modify the SCE. Run
it before `TRIPODES_assess_splicing` or `TRIPODES_run_velocity` to
confirm your object is usable.
