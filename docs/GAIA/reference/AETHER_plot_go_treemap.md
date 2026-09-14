# GO Term Treemap via Semantic Similarity Reduction

Reduces a set of enriched GO terms to representative clusters using
semantic similarity (via rrvgo), then renders a treemap where each tile
represents one GO term coloured by the module/gene-list that found it
most significantly. Tiles are grouped into larger semantic clusters
(e.g., "immune response") based on GO hierarchy.

## Usage

``` r
AETHER_plot_go_treemap(
  enrichment_result,
  orgdb,
  ont = "BP",
  top_n = 20L,
  threshold = 0.7,
  method = "Rel",
  show_shared = FALSE,
  show_group_labels = TRUE,
  min_label_scale = 0.1,
  palette = NULL,
  title = NULL,
  verbose = TRUE
)
```

## Arguments

- enrichment_result:

  A `gost_enrichment` object from
  [`APOLLO_enrich_gost`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)
  or
  [`DEMETER_load_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_load_enrichment.md),
  or a plain data.frame with at least the columns `term_id`, `source`,
  `p_value`, and `module`.

- orgdb:

  An OrgDb object (e.g. `org.Hs.eg.db`) or a package name string (e.g.
  `"org.Hs.eg.db"`) passed to `rrvgo::calculateSimMatrix()` and
  `rrvgo::reduceSimMatrix()`. Required — load the package first with
  [`library(org.Hs.eg.db)`](https://rdrr.io/r/base/library.html) (human)
  or [`library(org.Mm.eg.db)`](https://rdrr.io/r/base/library.html)
  (mouse), then pass the object.

- ont:

  Character vector. GO ontology/ontologies to plot. One or more of
  `"BP"`, `"MF"`, `"CC"`. Multiple values produce one panel per ontology
  in a single figure. Default: `"BP"`.

- top_n:

  Integer. Top N GO terms per module (by p-value) pooled before semantic
  reduction. Higher values include more terms but may slow down
  `rrvgo::calculateSimMatrix()`. Default: `20L`.

- threshold:

  Numeric (0–1). Similarity threshold for `rrvgo::reduceSimMatrix()`.
  Higher values → fewer, broader clusters; lower values → more, finer
  clusters. Default: `0.7`.

- method:

  Character. Semantic similarity measure passed to
  `rrvgo::calculateSimMatrix()`. One of `"Rel"` (default), `"Wang"`,
  `"Lin"`, `"Resnik"`, `"Jiang"`.

- show_shared:

  Logical. If `FALSE` (default), each GO term is assigned to the single
  module with the best (lowest) p-value for that term ("winner takes
  all"). If `TRUE`, terms found in multiple modules are assigned to a
  combined group labelled `"C1-C2"` (module names joined by `"-"`,
  sorted alphabetically) and given a distinct color automatically;
  useful for spotting convergent biology across modules.

- show_group_labels:

  Logical. If `TRUE` (default), each large parent-cluster tile shows the
  module group attribution beneath the semantic cluster name, e.g.
  `"immune response\n(C7)"`. The group shown is the one assigned to the
  cluster's representative (highest-scoring) term. Set to `FALSE` to
  display the semantic cluster name only.

- min_label_scale:

  Numeric (0–1). Minimum text scaling fraction passed to `treemap`'s
  `lowerbound.cex.labels`. When a tile is too small to fit text at the
  full `fontsize.labels`, treemap scales the text down; if the required
  scale falls below this value the label is suppressed entirely. Default
  `0.1` is more permissive than `treemap`'s own default of `0.4`, so
  more labels appear in small tiles. Set to `0` to force all labels
  regardless of tile size (may produce very small text). Opening a
  larger graphics device before calling this function also helps, as
  larger tiles require less scaling.

- palette:

  Named character vector. Module name → hex color mapping. If `NULL`
  (default), colors are auto-assigned from an internal 20-color
  categorical palette. Unknown modules fall back to auto-generated
  colors. The combined/shared group always uses `"#AAAAAA"` regardless
  of the palette.

- title:

  Character. Main plot title. If `NULL` (default), the ontology label is
  used (e.g., "GO Biological Process"). When multiple ontologies are
  requested, the ontology label is appended automatically.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

Invisibly returns a named list of `reducedTerms` data.frames (one per
ontology, named by ontology code: "BP", "MF", "CC"). The primary output
is the treemap rendered to the current graphics device. Save with
[`png()`](https://rdrr.io/r/grDevices/png.html) /
[`pdf()`](https://rdrr.io/r/grDevices/pdf.html) wrappers before calling.

## Details

**Workflow:**

1.  Filters the combined enrichment table to the requested GO source(s).

2.  Takes the top `top_n` terms per module by p-value and pools them.

3.  For each unique GO term, uses the best (lowest) p-value across
    modules as the rrvgo score (`-log10(p)`).

4.  Computes a semantic similarity matrix via
    `rrvgo::calculateSimMatrix()` (requires the relevant OrgDb package).

5.  Reduces to representative clusters via `rrvgo::reduceSimMatrix()`.

6.  Assigns each term a module group (winner or shared, per
    `show_shared`).

7.  Renders via `treemap::treemap()` with parent clusters as outer
    groupings and individual terms as inner tiles coloured by module
    group.

**OrgDb auto-detection:** The organism is read from
`enrichment_result$metadata$organism` (set automatically by
[`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)).
Supported organisms and their OrgDb packages: `hsapiens` →
`org.Hs.eg.db`; `mmusculus` → `org.Mm.eg.db`; `rnorvegicus` →
`org.Rn.eg.db`; `drerio` → `org.Dr.eg.db`; `dmelanogaster` →
`org.Dm.eg.db`; `celegans` → `org.Ce.eg.db`; `scerevisiae` →
`org.Sc.sgd.db`; `sscrofa` → `org.Ss.eg.db`; `btaurus` → `org.Bt.eg.db`;
`cfamiliaris` → `org.Cf.eg.db`; `ggallus` → `org.Gg.eg.db`. The
corresponding OrgDb package must be installed.

**Required packages:** rrvgo (Bioconductor) and treemap (CRAN), plus the
appropriate OrgDb package for the organism.

## See also

[`APOLLO_enrich_gost`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md),
[`DEMETER_load_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_load_enrichment.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(org.Hs.eg.db)  # human; use org.Mm.eg.db for mouse

# Basic GO:BP treemap
AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db)

# Multiple ontologies side by side
AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db, ont = c("BP", "MF"))

# Show shared terms across modules in grey
AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db, show_shared = TRUE)

# Custom module palette + save to PNG
my_pal <- c(C1 = "#E41A1C", C2 = "#377EB8", C3 = "#4DAF4A")
png("treemap_BP.png", width = 10, height = 8, units = "in", res = 300)
AETHER_plot_go_treemap(gost_res, orgdb = org.Hs.eg.db, ont = "BP", palette = my_pal)
dev.off()
} # }
```
