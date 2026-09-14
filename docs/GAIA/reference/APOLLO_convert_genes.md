# Convert gene IDs with fallback

Converts gene IDs from one type to another (e.g., ENSEMBL to SYMBOL),
preserving the original ID (optionally stripped of version numbers) when
conversion fails. Useful for converting count matrix rownames while
keeping all genes.

## Usage

``` r
APOLLO_convert_genes(
  genes,
  org_db,
  from_type = "ENSEMBL",
  to_type = "SYMBOL",
  strip_version = TRUE,
  verbose = TRUE
)
```

## Arguments

- genes:

  Character vector of gene IDs (e.g., rownames of a count matrix).

- org_db:

  OrgDb object for ID conversion (e.g., org.Hs.eg.db).

- from_type:

  Type of input gene IDs. Default: "ENSEMBL".

- to_type:

  Type of output gene IDs. Default: "SYMBOL".

- strip_version:

  Logical. Remove version numbers from Ensembl-style IDs before
  conversion (e.g., ENSG00000141510.16 -\> ENSG00000141510). Default:
  TRUE.

- verbose:

  Logical. Print conversion statistics. Default: TRUE.

## Value

Character vector of converted IDs, same length as input. Genes that
couldn't be converted retain their original ID (stripped if requested).
Has attribute "conversion_stats" with success/failure counts.

## Details

Unlike APOLLO_prepare_genelist() which drops unmapped genes (suitable
for enrichment analysis), this function preserves all genes - useful
when you need to maintain the same number of rows in a count matrix.

When duplicates occur after conversion (multiple ENSEMBL IDs mapping to
the same SYMBOL), a suffix is added to make names unique (e.g., "TP53",
"TP53.1").

## Examples

``` r
if (FALSE) { # \dontrun{
library(org.Hs.eg.db)

# Convert ENSEMBL to SYMBOL for count matrix rownames
ensembl_ids <- c("ENSG00000141510.16", "ENSG00000012048.23", "ENSG00000000003.15")
symbols <- APOLLO_convert_genes(ensembl_ids, org.Hs.eg.db)

# Apply to count matrix
rownames(counts) <- APOLLO_convert_genes(rownames(counts), org.Hs.eg.db)

} # }
```
