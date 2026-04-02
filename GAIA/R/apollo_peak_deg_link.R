# ==============================================================================
# APOLLO - Peak-DEG Linking
# ==============================================================================
# Identifies ATAC-seq peaks located within a defined distance of the TSS of
# differentially expressed genes. Window-based overlap captures all peak-gene
# pairs within range (many-to-many), not just the nearest gene per peak.
# ==============================================================================


#' Link Differentially Expressed Genes to Nearby ATAC Peaks
#'
#' @description Identifies ATAC-seq peaks located within a defined distance
#' of the TSS of differentially expressed genes (DEGs). Uses a window-based
#' approach that captures all peak-gene pairs within range (many-to-many),
#' not just the single nearest gene per peak.
#'
#' Significance filtering should be applied to \code{deg_result} by the caller
#' before passing it to this function.
#'
#' @param deg_result DEA result. Either an \code{artemis_dea} S3 object from
#'   \code{ARTEMIS_perform_dea()}, or a plain data.frame. The data.frame must
#'   contain at minimum a gene name column (see \code{gene_col}).
#' @param peaks Peak data. Either a data.frame with at least three columns
#'   (chr, start, end — BED-style 0-based starts) or a \code{GRanges} object.
#' @param gtf_file Character. Path to GTF/GFF annotation file. Used to extract
#'   TSS positions. Either \code{gtf_file} or \code{txdb} must be provided.
#'   When both are given, \code{gtf_file} takes precedence.
#' @param txdb A TxDb object (from \code{APOLLO_make_txdb()}). Alternative to
#'   \code{gtf_file}. Gene IDs in the TxDb must match the values in
#'   \code{gene_col} of \code{deg_result} (e.g. both use gene symbols).
#' @param distance Numeric. Maximum distance in bp from a peak (any edge) to
#'   a TSS for the pair to be included. Default \code{50000} (50 kb).
#' @param gene_col Character. Column in \code{deg_result} holding gene symbols.
#'   Default \code{"gene_name"}.
#' @param chr_mapping Character or named vector. Chromosome renaming to apply
#'   to gene coordinates from the GTF so they match peak coordinates. Pass
#'   \code{"T2T"} for the T2T-CHM13 NCBI→UCSC mapping. Ignored when using
#'   \code{txdb} (mapping is applied at TxDb creation time).
#'   Default \code{NULL} (no renaming).
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An S3 object of class \code{"apollo_deg_peaks"} (a named list):
#' \describe{
#'   \item{\code{$linked}}{data.frame. One row per peak × DEG pair within
#'     \code{distance}. Contains all columns from \code{deg_result}, plus peak
#'     coordinates, TSS position, and \code{distance_to_tss} (signed: negative
#'     = peak is upstream of TSS relative to transcription direction, positive
#'     = downstream).}
#'   \item{\code{$deg_summary}}{data.frame. One row per input DEG. Includes
#'     all \code{deg_result} columns plus \code{has_peak} (logical),
#'     \code{n_peaks} (integer), and \code{nearest_peak_distance} (bp).}
#'   \item{\code{$peak_summary}}{data.frame. One row per input peak. Includes
#'     peak coordinates plus \code{near_deg} (logical), \code{n_degs},
#'     \code{nearest_deg} (gene name), and \code{nearest_deg_distance}.}
#'   \item{\code{$stats}}{Named list of summary counts.}
#'   \item{\code{$params}}{Named list of parameters used.}
#' }
#'
#' @details
#' **Signed distance convention** (strand-aware):
#' \itemize{
#'   \item + strand: \code{distance_to_tss = peak_midpoint - TSS_position}
#'   \item - strand: \code{distance_to_tss = TSS_position - peak_midpoint}
#' }
#' Negative values = peak is upstream (regulatory region); positive = downstream
#' (gene body or beyond).
#'
#' **TSS extraction:** When using \code{gtf_file}, only gene-level records are
#' read, making import fast even for large annotation files. The most 5' TSS
#' per gene symbol is used (in transcription direction). Gene name is read from
#' the \code{gene_name} or \code{gene} GTF attribute (tried in that order).
#'
#' @export
APOLLO_link_degs_to_peaks <- function(deg_result,
                                       peaks,
                                       gtf_file    = NULL,
                                       txdb        = NULL,
                                       distance    = 50000,
                                       gene_col    = "gene_name",
                                       chr_mapping = NULL,
                                       verbose     = TRUE) {

  # --- Package checks ---------------------------------------------------------
  for (pkg in c("GenomicRanges", "IRanges")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' required: BiocManager::install('", pkg, "')",
           call. = FALSE)
  }

  if (is.null(gtf_file) && is.null(txdb))
    stop("Provide either gtf_file or txdb to determine TSS positions.",
         call. = FALSE)

  # --- Extract results data.frame ---------------------------------------------
  deg_df <- .apollo_extract_deg_df(deg_result)

  if (!gene_col %in% colnames(deg_df))
    stop("gene_col '", gene_col, "' not found in deg_result.", call. = FALSE)

  if (nrow(deg_df) == 0)
    stop("deg_result contains no rows.", call. = FALSE)

  deg_genes <- deg_df[[gene_col]]
  deg_genes <- deg_genes[!is.na(deg_genes) & nchar(deg_genes) > 0]

  if (verbose) message("[APOLLO] Input genes: ", length(deg_genes))

  # --- Build peaks GRanges (0-based BED → 1-based GRanges) -------------------
  if (inherits(peaks, "GRanges")) {
    peaks_gr  <- peaks
    peaks_df  <- as.data.frame(peaks)[, c("seqnames", "start", "end"),
                                       drop = FALSE]
    colnames(peaks_df) <- c("chr", "start", "end")
    peaks_df$start <- peaks_df$start - 1L  # back to 0-based for output
  } else {
    peaks_df <- as.data.frame(peaks)
    if (!all(c("chr", "start", "end") %in% colnames(peaks_df)))
      stop("peaks must have columns: chr, start, end", call. = FALSE)
    peaks_gr <- GenomicRanges::GRanges(
      seqnames = peaks_df$chr,
      ranges   = IRanges::IRanges(
        start = as.integer(peaks_df$start) + 1L,  # 0-based → 1-based
        end   = as.integer(peaks_df$end)
      )
    )
  }
  peaks_df$.peak_idx <- seq_len(nrow(peaks_df))
  n_peaks <- nrow(peaks_df)

  if (verbose) message("[APOLLO] Input peaks: ", n_peaks)

  # --- Get TSS positions for input genes --------------------------------------
  tss_df <- if (!is.null(gtf_file)) {
    .apollo_tss_from_gtf(gtf_file, deg_genes, gene_col, chr_mapping, verbose)
  } else {
    .apollo_tss_from_txdb(txdb, deg_genes, gene_col, verbose)
  }

  tss_df <- tss_df[tss_df[[gene_col]] %in% deg_genes, , drop = FALSE]

  n_matched <- length(intersect(deg_genes, tss_df[[gene_col]]))
  n_missing <- length(setdiff(deg_genes, tss_df[[gene_col]]))
  if (verbose) {
    message("[APOLLO] Genes with TSS found: ", n_matched,
            if (n_missing > 0) paste0(" (", n_missing, " not in annotation)") else "")
  }

  if (nrow(tss_df) == 0)
    stop("No genes could be matched to coordinates in the annotation. ",
         "Check gene names and chr_mapping.", call. = FALSE)

  # --- Build TSS windows GRanges (±distance) ----------------------------------
  tss_gr <- GenomicRanges::GRanges(
    seqnames = tss_df$chr,
    ranges   = IRanges::IRanges(
      start = pmax(1L, as.integer(tss_df$tss) - as.integer(distance)),
      end   = as.integer(tss_df$tss) + as.integer(distance)
    ),
    strand = tss_df$strand
  )

  # --- Find overlaps ----------------------------------------------------------
  hits <- GenomicRanges::findOverlaps(peaks_gr, tss_gr, ignore.strand = TRUE)

  if (length(hits) == 0) {
    if (verbose) message("[APOLLO] No peaks found within ", distance,
                         " bp of any gene TSS.")
    return(.apollo_empty_result(deg_df, peaks_df, gene_col, distance,
                                chr_mapping))
  }

  if (verbose) message("[APOLLO] Peak-gene pairs within ", distance,
                       " bp: ", length(hits))

  # --- Compute signed distance (peak midpoint → TSS) --------------------------
  peak_idx <- S4Vectors::queryHits(hits)
  gene_idx <- S4Vectors::subjectHits(hits)

  peak_mid    <- (GenomicRanges::start(peaks_gr[peak_idx]) +
                    GenomicRanges::end(peaks_gr[peak_idx])) / 2
  tss_pos     <- tss_df$tss[gene_idx]
  strand      <- tss_df$strand[gene_idx]

  # Signed: negative = upstream, positive = downstream (transcription-relative)
  signed_dist <- ifelse(strand == "-",
                        tss_pos - peak_mid,
                        peak_mid - tss_pos)

  # --- Build linked data.frame ------------------------------------------------
  linked <- data.frame(
    tss_df[gene_idx, , drop = FALSE],
    chr_peak        = peaks_df$chr[peak_idx],
    start_peak      = peaks_df$start[peak_idx],
    end_peak        = peaks_df$end[peak_idx],
    distance_to_tss = round(signed_dist),
    abs_distance    = abs(round(signed_dist)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Merge in all DEG columns
  deg_stat_cols <- setdiff(colnames(deg_df), gene_col)
  linked <- merge(linked, deg_df[, c(gene_col, deg_stat_cols), drop = FALSE],
                  by = gene_col, all.x = TRUE)
  linked <- linked[order(linked$abs_distance), , drop = FALSE]
  rownames(linked) <- NULL

  # --- DEG summary (one row per DEG) ------------------------------------------
  linked_genes <- split(linked$abs_distance, linked[[gene_col]])

  deg_summary <- deg_df
  deg_summary$has_peak <- deg_summary[[gene_col]] %in% names(linked_genes)
  deg_summary$n_peaks  <- vapply(deg_summary[[gene_col]], function(g) {
    if (g %in% names(linked_genes)) length(linked_genes[[g]]) else 0L
  }, integer(1))
  deg_summary$nearest_peak_distance <- vapply(deg_summary[[gene_col]], function(g) {
    if (g %in% names(linked_genes)) min(linked_genes[[g]]) else NA_real_
  }, numeric(1))
  deg_summary <- deg_summary[order(deg_summary$nearest_peak_distance,
                                   na.last = TRUE), , drop = FALSE]
  rownames(deg_summary) <- NULL

  # --- Peak summary (one row per peak) ----------------------------------------
  linked$.peak_key   <- paste0(linked$chr_peak, ":", linked$start_peak, "-",
                                linked$end_peak)
  peaks_df$.peak_key <- paste0(peaks_df$chr, ":", peaks_df$start, "-",
                                peaks_df$end)

  peak_link_by_key <- split(
    linked[, c(gene_col, "abs_distance"), drop = FALSE],
    linked$.peak_key
  )

  peak_summary <- peaks_df[, c("chr", "start", "end"), drop = FALSE]
  peak_summary$near_deg <- peaks_df$.peak_key %in% names(peak_link_by_key)

  peak_summary$n_degs <- vapply(peaks_df$.peak_key, function(k) {
    if (k %in% names(peak_link_by_key)) nrow(peak_link_by_key[[k]]) else 0L
  }, integer(1))

  peak_summary$nearest_deg <- vapply(peaks_df$.peak_key, function(k) {
    if (k %in% names(peak_link_by_key)) {
      df <- peak_link_by_key[[k]]
      df[[gene_col]][which.min(df$abs_distance)]
    } else NA_character_
  }, character(1))

  peak_summary$nearest_deg_distance <- vapply(peaks_df$.peak_key, function(k) {
    if (k %in% names(peak_link_by_key))
      min(peak_link_by_key[[k]]$abs_distance)
    else NA_real_
  }, numeric(1))

  linked$.peak_key   <- NULL
  peaks_df$.peak_key <- NULL
  rownames(peak_summary) <- NULL

  # --- Stats ------------------------------------------------------------------
  stats <- list(
    n_genes_input     = nrow(deg_df),
    n_genes_with_peaks = sum(deg_summary$has_peak),
    n_genes_no_peaks  = sum(!deg_summary$has_peak),
    n_peaks_input     = n_peaks,
    n_peaks_near_genes = sum(peak_summary$near_deg),
    n_pairs           = nrow(linked),
    distance_used     = distance
  )

  if (verbose) {
    message("[APOLLO] Genes with ≥1 nearby peak: ", stats$n_genes_with_peaks,
            " / ", stats$n_genes_input,
            " (", round(100 * stats$n_genes_with_peaks / stats$n_genes_input, 1), "%)")
    message("[APOLLO] Peaks near ≥1 gene: ", stats$n_peaks_near_genes,
            " / ", stats$n_peaks_input)
    message("[APOLLO] Total peak-gene pairs: ", stats$n_pairs)
  }

  # --- Return -----------------------------------------------------------------
  result <- list(
    linked       = linked,
    deg_summary  = deg_summary,
    peak_summary = peak_summary,
    stats        = stats,
    params       = list(
      distance    = distance,
      gene_col    = gene_col,
      chr_mapping = chr_mapping
    )
  )
  class(result) <- c("apollo_deg_peaks", "list")
  result
}


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

# Extract a flat data.frame from artemis_dea S3 or plain data.frame
.apollo_extract_deg_df <- function(deg_result) {
  if (inherits(deg_result, "artemis_dea") || inherits(deg_result, "list")) {
    if ("results" %in% names(deg_result) && is.data.frame(deg_result$results))
      return(deg_result$results)
  }
  if (is.data.frame(deg_result)) return(deg_result)
  stop("deg_result must be an artemis_dea object or a data.frame.", call. = FALSE)
}


# Extract TSS positions from GTF file (only reads gene-level records — fast)
.apollo_tss_from_gtf <- function(gtf_file, deg_genes, gene_col,
                                  chr_mapping, verbose) {
  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("Package 'rtracklayer' required: BiocManager::install('rtracklayer')",
         call. = FALSE)

  if (!file.exists(gtf_file))
    stop("gtf_file not found: ", gtf_file, call. = FALSE)

  if (verbose) message("[APOLLO] Reading gene coordinates from GTF...")

  genes_gr <- tryCatch(
    rtracklayer::import(gtf_file, feature.type = "gene"),
    error = function(e) stop("Failed to read GTF: ", conditionMessage(e),
                             call. = FALSE)
  )

  # Determine which mcols column holds the gene symbol
  mc      <- as.data.frame(GenomicRanges::mcols(genes_gr))
  sym_col <- intersect(c("gene_name", "gene"), colnames(mc))[1]
  if (is.na(sym_col))
    stop("GTF gene records have neither 'gene_name' nor 'gene' attribute. ",
         "Use txdb with matching gene IDs instead.", call. = FALSE)

  gene_syms <- as.character(mc[[sym_col]])

  # Apply chromosome name mapping if requested
  if (!is.null(chr_mapping)) {
    if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
      stop("GenomeInfoDb required for chr_mapping: BiocManager::install('GenomeInfoDb')",
           call. = FALSE)
    map_vec <- if (is.character(chr_mapping) && length(chr_mapping) == 1 &&
                   is.null(names(chr_mapping))) {
      APOLLO_get_chr_mapping(chr_mapping)
    } else {
      chr_mapping
    }
    cur        <- GenomeInfoDb::seqlevels(genes_gr)
    to_rename  <- intersect(cur, names(map_vec))
    if (length(to_rename) > 0) {
      rv <- map_vec[to_rename]; names(rv) <- to_rename
      genes_gr <- GenomeInfoDb::renameSeqlevels(genes_gr, rv)
      if (verbose) message("[APOLLO]   Renamed ", length(to_rename),
                           " chromosomes in GTF gene coords")
    }
  }

  seqnm  <- as.character(GenomicRanges::seqnames(genes_gr))
  starts <- GenomicRanges::start(genes_gr)
  ends   <- GenomicRanges::end(genes_gr)
  strnds <- as.character(GenomicRanges::strand(genes_gr))
  tss_pos <- ifelse(strnds == "-", ends, starts)

  raw_df <- data.frame(
    gene_name = gene_syms,
    chr       = seqnm,
    tss       = tss_pos,
    strand    = strnds,
    stringsAsFactors = FALSE
  )
  colnames(raw_df)[1] <- gene_col

  raw_df <- raw_df[raw_df[[gene_col]] %in% deg_genes, , drop = FALSE]
  if (nrow(raw_df) == 0) return(raw_df)

  # For genes with multiple records, keep the most 5' TSS
  most5prime <- function(df) {
    if (nrow(df) == 1) return(df)
    if (df$strand[1] == "-") df[which.max(df$tss), ] else df[which.min(df$tss), ]
  }
  result <- do.call(rbind, lapply(split(raw_df, raw_df[[gene_col]]), most5prime))
  rownames(result) <- NULL
  result
}


# Extract TSS positions from a TxDb object
.apollo_tss_from_txdb <- function(txdb, deg_genes, gene_col, verbose) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
    stop("GenomicFeatures required: BiocManager::install('GenomicFeatures')",
         call. = FALSE)

  if (verbose) message("[APOLLO] Extracting TSS positions from TxDb...")

  genes_gr <- tryCatch(
    GenomicFeatures::genes(txdb),
    error = function(e) stop("Failed to extract genes from TxDb: ",
                             conditionMessage(e), call. = FALSE)
  )

  strnds  <- as.character(GenomicRanges::strand(genes_gr))
  tss_pos <- ifelse(strnds == "-",
                    GenomicRanges::end(genes_gr),
                    GenomicRanges::start(genes_gr))

  raw_df <- data.frame(
    gene_id = names(genes_gr),
    chr     = as.character(GenomicRanges::seqnames(genes_gr)),
    tss     = tss_pos,
    strand  = strnds,
    stringsAsFactors = FALSE
  )
  colnames(raw_df)[1] <- gene_col

  raw_df <- raw_df[raw_df[[gene_col]] %in% deg_genes, , drop = FALSE]

  if (nrow(raw_df) == 0)
    warning("No genes matched IDs in TxDb. TxDb gene IDs may not be gene ",
            "symbols — consider using gtf_file instead.", call. = FALSE)

  rownames(raw_df) <- NULL
  raw_df
}


# Build an empty apollo_deg_peaks result when no overlaps are found
.apollo_empty_result <- function(deg_df, peaks_df, gene_col, distance,
                                  chr_mapping) {
  deg_df$has_peak              <- FALSE
  deg_df$n_peaks               <- 0L
  deg_df$nearest_peak_distance <- NA_real_

  peaks_df$near_deg             <- FALSE
  peaks_df$n_degs               <- 0L
  peaks_df$nearest_deg          <- NA_character_
  peaks_df$nearest_deg_distance <- NA_real_

  linked_cols <- c(gene_col, "chr", "tss", "strand",
                   "chr_peak", "start_peak", "end_peak",
                   "distance_to_tss", "abs_distance")
  linked <- as.data.frame(matrix(nrow = 0, ncol = length(linked_cols),
                                 dimnames = list(NULL, linked_cols)))

  result <- list(
    linked       = linked,
    deg_summary  = deg_df,
    peak_summary = peaks_df[, c("chr", "start", "end", "near_deg", "n_degs",
                                "nearest_deg", "nearest_deg_distance"),
                             drop = FALSE],
    stats = list(
      n_genes_input      = nrow(deg_df),
      n_genes_with_peaks = 0L,
      n_genes_no_peaks   = nrow(deg_df),
      n_peaks_input      = nrow(peaks_df),
      n_peaks_near_genes = 0L,
      n_pairs            = 0L,
      distance_used      = distance
    ),
    params = list(distance = distance, gene_col = gene_col,
                  chr_mapping = chr_mapping)
  )
  class(result) <- c("apollo_deg_peaks", "list")
  result
}


#' Print method for apollo_deg_peaks
#' @export
print.apollo_deg_peaks <- function(x, ...) {
  s <- x$stats
  cat("DEG peaks\n")
  cat("------------------------------\n")
  cat("Distance window   : \u00B1", format(s$distance_used, big.mark = ","),
      " bp\n", sep = "")
  cat("Genes input       : ", s$n_genes_input, "\n", sep = "")
  cat("Genes with peaks  : ", s$n_genes_with_peaks,
      " (", round(100 * s$n_genes_with_peaks / max(s$n_genes_input, 1), 1),
      "%)\n", sep = "")
  cat("Peaks near genes  : ", s$n_peaks_near_genes, " / ",
      s$n_peaks_input, "\n", sep = "")
  cat("Total pairs       : ", s$n_pairs, "\n", sep = "")
  invisible(x)
}
