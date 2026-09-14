# Plot RNA velocity as smooth streamlines over a KDE cluster background

Replicates the `scv.pl.velocity_embedding_stream` style. Per-cluster KDE
density blobs provide the colour context; streamlines show the inferred
directional flow of RNA velocity across the embedding.

## Usage

``` r
ASPIS_plot_velocity_stream(
  sce,
  dimred = "UMAP",
  colour_by = "cluster",
  grid_res = 30,
  n_streams = 200,
  stream_length = 50,
  stream_colour = "black",
  stream_size = 0.4,
  stream_alpha = 0.8,
  bandwidth = NULL,
  min_density = 0.05,
  blob_alpha = 0.25,
  blob_bins = 4,
  show_points = FALSE,
  point_size = 0.3,
  point_alpha = 0.3,
  min_vel_frac = 0.05,
  palette = NULL,
  assay_name = "logcounts",
  title = NULL
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with velocity results from
  [`TRIPODES_run_velocity`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_velocity.md).

- dimred:

  Character. Embedding to use. Default `"UMAP"`.

- colour_by:

  Character. `colData` column used for blob and (optionally) streamline
  colour. Must be discrete. Default `"cluster"`.

- grid_res:

  Integer. Resolution of the velocity grid (grid_res × grid_res cells).
  Higher values give smoother fields but are slower to compute. Default
  `30`.

- n_streams:

  Integer. Number of streamline seed points. Default `200`.

- stream_length:

  Integer. Maximum integration steps per streamline. Default `50`.

- stream_colour:

  Character colour string, or `NULL`. A fixed string (e.g. `"black"`)
  renders all streamlines in that colour. `NULL` colours each streamline
  by the nearest cell's `colour_by` value, matching the blob palette.
  Default `"black"`.

- stream_size:

  Numeric. Streamline line width. Default `0.4`.

- stream_alpha:

  Numeric. Streamline opacity. Default `0.8`.

- bandwidth:

  Numeric or `NULL`. Gaussian kernel bandwidth for velocity field
  interpolation. `NULL` auto-scales to grid spacing.

- min_density:

  Numeric in \[0, 1\]. Grid points whose local cell density falls below
  this fraction of the grid's peak density are masked out of the
  velocity field entirely (no seeding, streamlines cannot cross them).
  Without this mask, the Gaussian kernel's infinite support means
  distant grid points (e.g. in empty corners of a non-convex embedding)
  still get a spuriously "valid" direction dominated by whichever cells
  are least-far-away, producing streamlines that wander through empty
  space with no relation to the actual data. Raise this value if
  streamlines still appear disconnected from the point cloud; lower it
  if streamlines are being cut off too early near the edges of dense
  clusters. Default `0.05`.

- blob_alpha:

  Numeric. Maximum opacity of KDE polygon fills. Default `0.25`.

- blob_bins:

  Integer. Number of density contour levels per cluster. Fewer bins
  produce softer blobs. Default `4`.

- show_points:

  Logical. Whether to overlay individual cell points. Default `FALSE`.

- point_size:

  Numeric. Cell point size (if `show_points=TRUE`). Default `0.3`.

- point_alpha:

  Numeric. Cell point opacity. Default `0.3`.

- min_vel_frac:

  Numeric. Fraction of mean velocity magnitude below which a streamline
  step is considered stalled and tracing stops. Default `0.05`.

- palette:

  Character vector or `NULL` for defaults.

- assay_name:

  Character. Assay for gene expression. Default `"logcounts"`.

- title:

  Character or `NULL`. Plot title.

## Value

A `ggplot` object.

## Details

Streamlines are computed by interpolating per-cell velocity vectors onto
a regular grid (Gaussian kernel weighting) and then integrating forward
along unit-normalised direction vectors using Euler steps.
