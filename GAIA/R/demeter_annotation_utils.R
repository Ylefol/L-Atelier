# GAIA/Demeter/annotation_utils.R
# General-purpose annotation helpers that attach external data to result data.frames.
#
# These utilities are module-agnostic: they operate on any data.frame that
# carries the expected columns and do not depend on peak or omics-specific
# infrastructure.


# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Add total exonic gene length to a result data.frame
#'
#' Reads a GTF file, sums exon lengths per gene, and appends a \code{gene_length}
#' column to \code{df} by matching on \code{gene_id}.  Rows whose gene ID is
#' absent from the GTF receive \code{NA}.
#'
#' Used internally by \code{\link{AETHER_plot_gene_length_distribution}}.
#'
#' @param df Data.frame with a \code{gene_id} column.
#' @param gtf_file Path to a GTF/GFF annotation file.
#'
#' @return \code{df} with an additional \code{gene_length} column (integer,
#'   total exonic bp per gene).
#' @noRd
.demeter_add_gene_length <- function(df, gtf_file) {

  if (!file.exists(gtf_file))
    stop("GTF file not found: ", gtf_file, call. = FALSE)

  if (!"gene_id" %in% colnames(df))
    stop("'df' must contain a 'gene_id' column.", call. = FALSE)

  g <- rtracklayer::readGFF(gtf_file)
  g <- g[g$type == "exon", ]
  g <- g[g$gene_id %in% df$gene_id, ]

  if (nrow(g) == 0)
    stop("No exon records found in GTF matching gene IDs in 'df'. ",
         "Check that gene ID format matches (e.g. ENSG vs symbol).",
         call. = FALSE)

  g$exon_length <- abs(g$end - g$start)
  gene_lengths  <- aggregate(exon_length ~ gene_id, data = g, FUN = sum)
  colnames(gene_lengths)[colnames(gene_lengths) == "exon_length"] <- "gene_length"

  merge(df, gene_lengths, by = "gene_id", all.x = TRUE)
}
