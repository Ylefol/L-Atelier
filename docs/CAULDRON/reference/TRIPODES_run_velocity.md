# Estimate RNA velocity via scVelo

Runs RNA velocity estimation using the Python scVelo library, called
directly through a basilisk-managed environment and a bundled script
(`inst/python/scvelo_run.py`) — not via the velociraptor Bioconductor
wrapper (CAULDRON has no dependency on velociraptor). Spliced and
unspliced count assays must be present in the SCE. Results are stored
back into the input SCE and returned.

## Usage

``` r
TRIPODES_run_velocity(
  sce,
  spliced_assay = "spliced",
  unspliced_assay = "unspliced",
  mode = c("deterministic", "stochastic", "dynamical"),
  use_dimred = "PCA",
  n_pcs = 30L,
  n_neighbors = 30L,
  store_moments = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- sce:

  A `SingleCellExperiment`; must contain spliced and unspliced count
  assays.

- spliced_assay:

  Name of the spliced counts assay (default `"spliced"`).

- unspliced_assay:

  Name of the unspliced counts assay (default `"unspliced"`).

- mode:

  scVelo fitting mode. One of `"deterministic"` (default, fast
  first-order OLS fit), `"stochastic"` (accounts for transcriptional
  noise via second-order moments; requires more memory), or
  `"dynamical"` (most accurate, substantially slower — fits full kinetic
  model per gene via `scv.tl.recover_dynamics()`, run automatically
  first when this mode is selected).

- use_dimred:

  Name of the PCA-like dimensionality reduction in the SCE to pass to
  scVelo for neighbour graph construction (default `"PCA"`). This should
  be the same reduction used to build your UMAP, ensuring that velocity
  neighbours are consistent with the embedding. `n_pcs` is automatically
  capped to the number of components available. Set to `NULL` to let
  scVelo compute its own PCA internally (not recommended if you intend
  to project velocity onto an existing UMAP).

- n_pcs:

  Number of PCs to use for neighbour graph construction within scVelo
  (default `30L`).

- n_neighbors:

  Number of neighbours for the scVelo graph (default `30L`).

- store_moments:

  Logical; also store the moments-smoothed spliced/ unspliced matrices
  as `assay(sce, "Ms")`/`assay(sce, "Mu")` (default `FALSE`). These are
  computed internally by scVelo regardless (needed to fit velocity
  itself); this only controls whether they're additionally returned and
  stored, since doing so roughly doubles the SCE's assay footprint.
  Enable for phase-portrait-style diagnostics (unspliced vs. spliced per
  gene).

- verbose:

  Logical; print progress messages (default `TRUE`).

- ...:

  Currently unused — kept for backward-compatible signature stability;
  nothing in the current implementation forwards these arguments
  anywhere.

## Value

The input SCE with velocity results added (see Details).
`metadata(sce)$velocity_run` is set to `TRUE` and
`metadata(sce)$velocity_mode` records the mode used.

## Details

If `TRIPODES_assess_splicing` has not been run, a warning is issued but
execution continues. If the SCE fails the hard readiness checks (missing
assays, dimension mismatch, continuous data), execution stops with a
message directing the user to `TRIPODES_check_velocity_ready`.

Stored results:

- `assay(sce, "velocity")`:

  Gene × cell velocity matrix.

- `assay(sce, "Ms")`, `assay(sce, "Mu")`:

  Moments-smoothed spliced/unspliced counts (gene × cell), as fit
  against by `scv.pp.moments()`. Only stored when `store_moments = TRUE`
  (default `FALSE`) — these are full dense gene × cell matrices and
  substantially increase the SCE's size.

- `colData(sce)$velocity_confidence`:

  Per-cell velocity confidence score from scVelo.

- `colData(sce)$velocity_length`:

  Per-cell velocity vector magnitude.

- `reducedDim(sce, "velocity_<use_dimred>")`:

  Projected velocity coordinates on the chosen embedding (if
  `use_dimred` is present in the SCE).
