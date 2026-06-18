# GAIA/Apollo/activity_networks.R
# Prior knowledge networks for activity inference
#
# Wrapper functions for accessing TF-target and pathway databases
# used by decoupleR for activity inference.
#
# Databases:
# - CollecTRI: Comprehensive TF-target network (default for TF activity)
# - DoRothEA: Curated TF regulons with confidence levels
# - PROGENy: Pathway responsive genes (for pathway activity)


# ==============================================================================
# TF-TARGET NETWORKS
# ==============================================================================

#' Get CollecTRI TF-target network
#'
#' Retrieves the CollecTRI transcription factor-target gene network.
#' CollecTRI is a comprehensive collection of TF-target interactions
#' compiled from multiple sources.
#'
#' @param organism Character. Organism to use: "human" or "mouse". Default: "human".
#' @param split_complexes Logical. If TRUE, splits TF complexes into individual
#'   TFs. Default: FALSE.
#' @param verbose Logical. Print information about the network. Default: TRUE.
#' @param ... Additional arguments passed to decoupleR::get_collectri().
#'
#' @return A tibble with columns:
#'   \item{source}{Transcription factor name}
#'   \item{target}{Target gene name}
#'   \item{mor}{Mode of regulation: 1 (activation) or -1 (repression)}
#'
#' @details
#' CollecTRI integrates TF-target interactions from:
#' - Literature curation
#' - ChIP-seq experiments
#' - Inference from gene expression
#'
#' This is the recommended default for TF activity inference as it provides
#' good coverage while maintaining reasonable accuracy.
#'
#' @examples
#' \dontrun{
#' # Get human TF-target network
#' network <- APOLLO_get_collectri()
#'
#' # Get mouse network
#' network_mouse <- APOLLO_get_collectri(organism = "mouse")
#'
#' }
#' @export
APOLLO_get_collectri <- function(organism = "human",
                                  split_complexes = FALSE,
                                  verbose = TRUE,
                                  ...) {

  # Validate organism

organism <- tolower(organism)
  if (!organism %in% c("human", "mouse")) {
    stop("organism must be 'human' or 'mouse'")
  }

  if (verbose) {
    cat("[APOLLO] Fetching CollecTRI network for", organism, "...\n")
  }

  # Get network
  network <- decoupleR::get_collectri(
    organism = organism,
    split_complexes = split_complexes,
    ...
  )

  if (verbose) {
    n_tfs <- length(unique(network$source))
    n_targets <- length(unique(network$target))
    n_interactions <- nrow(network)
    cat("   TFs:", n_tfs, "\n")
    cat("   Targets:", n_targets, "\n")
    cat("   Interactions:", n_interactions, "\n")
  }

  # Add metadata as attributes
attr(network, "database") <- "collectri"
  attr(network, "organism") <- organism
  attr(network, "network_type") <- "tf_target"

  return(network)
}


#' Get DoRothEA TF regulons
#'
#' Retrieves DoRothEA transcription factor regulons with confidence level
#' filtering. DoRothEA provides curated TF-target interactions with
#' assigned confidence scores.
#'
#' @param organism Character. Organism to use: "human" or "mouse". Default: "human".
#' @param levels Character vector. Confidence levels to include.
#'   Options: "A" (highest), "B", "C", "D", "E" (lowest).
#'   Default: c("A", "B", "C", "D", "E") (all levels).
#' @param verbose Logical. Print information about the network. Default: TRUE.
#' @param ... Additional arguments passed to decoupleR::get_dorothea().
#'
#' @return A tibble with columns:
#'   \item{source}{Transcription factor name}
#'   \item{target}{Target gene name}
#'   \item{mor}{Mode of regulation: 1 (activation) or -1 (repression)}
#'   \item{confidence}{Confidence level (A-E)}
#'
#' @details
#' DoRothEA confidence levels:
#' - A: Highest confidence (multiple evidence types)
#' - B: High confidence (literature + ChIP-seq)
#' - C: Medium confidence (ChIP-seq or inference)
#' - D: Low confidence (inference only)
#' - E: Lowest confidence (predicted)
#'
#' For most analyses, levels A-C are recommended. All levels are returned
#' by default to allow downstream filtering (e.g., at visualization).
#'
#' @examples
#' \dontrun{
#' # Get all DoRothEA regulons
#' network <- APOLLO_get_dorothea()
#'
#' # Get only high-confidence regulons
#' network_hc <- APOLLO_get_dorothea(levels = c("A", "B"))
#'
#' }
#' @export
APOLLO_get_dorothea <- function(organism = "human",
                                 levels = c("A", "B", "C", "D", "E"),
                                 verbose = TRUE,
                                 ...) {

  # Validate organism
  organism <- tolower(organism)
  if (!organism %in% c("human", "mouse")) {
    stop("organism must be 'human' or 'mouse'")
  }

  # Validate levels
  valid_levels <- c("A", "B", "C", "D", "E")
  levels <- toupper(levels)
  if (!all(levels %in% valid_levels)) {
    stop("levels must be one or more of: ", paste(valid_levels, collapse = ", "))
  }

  if (verbose) {
    cat("[APOLLO] Fetching DoRothEA regulons for", organism, "...\n")
    cat("   Confidence levels:", paste(levels, collapse = ", "), "\n")
  }

  # Get network
  network <- decoupleR::get_dorothea(
    organism = organism,
    levels = levels,
    ...
  )

  if (verbose) {
    n_tfs <- length(unique(network$source))
    n_targets <- length(unique(network$target))
    n_interactions <- nrow(network)
    cat("   TFs:", n_tfs, "\n")
    cat("   Targets:", n_targets, "\n")
    cat("   Interactions:", n_interactions, "\n")

    # Show breakdown by confidence
    if ("confidence" %in% colnames(network)) {
      level_counts <- table(network$confidence)
      cat("   By confidence level:\n")
      for (lvl in names(level_counts)) {
        cat("     ", lvl, ":", level_counts[lvl], "\n")
      }
    }
  }

  # Add metadata as attributes
  attr(network, "database") <- "dorothea"
  attr(network, "organism") <- organism
  attr(network, "levels") <- levels
  attr(network, "network_type") <- "tf_target"

  return(network)
}


# ==============================================================================
# PATHWAY NETWORKS
# ==============================================================================

#' Get PROGENy pathway signatures
#'
#' Retrieves PROGENy pathway-responsive gene signatures. PROGENy provides
#' genes that are responsive to pathway activity, derived from perturbation
#' experiments.
#'
#' @param organism Character. Organism to use: "human" or "mouse". Default: "human".
#' @param top Integer. Number of top responsive genes per pathway to include.
#'   Default: 500. Use higher values for more coverage, lower for more specificity.
#' @param verbose Logical. Print information about the network. Default: TRUE.
#' @param ... Additional arguments passed to decoupleR::get_progeny().
#'
#' @return A tibble with columns:
#'   \item{source}{Pathway name}
#'   \item{target}{Responsive gene name}
#'   \item{weight}{Gene weight/importance for pathway activity}
#'
#' @details
#' PROGENy pathways include:
#' - Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K,
#'   TGFb, TNFa, Trail, VEGF, WNT
#'
#' The `top` parameter controls how many genes per pathway are included.
#' More genes = better coverage but potentially more noise.
#' Fewer genes = more specific but may miss relevant signals.
#'
#' @examples
#' \dontrun{
#' # Get PROGENy signatures (default 500 genes per pathway)
#' network <- APOLLO_get_progeny()
#'
#' # Get more specific signatures (100 genes per pathway)
#' network_specific <- APOLLO_get_progeny(top = 100)
#'
#' }
#' @export
APOLLO_get_progeny <- function(organism = "human",
                                top = 500,
                                verbose = TRUE,
                                ...) {

  # Validate organism
  organism <- tolower(organism)
  if (!organism %in% c("human", "mouse")) {
    stop("organism must be 'human' or 'mouse'")
  }

  # Validate top
  if (!is.numeric(top) || top < 1) {
    stop("top must be a positive integer")
  }
  top <- as.integer(top)

  if (verbose) {
    cat("[APOLLO] Fetching PROGENy signatures for", organism, "...\n")
    cat("   Top genes per pathway:", top, "\n")
  }

  # Get network
  network <- decoupleR::get_progeny(
    organism = organism,
    top = top,
    ...
  )

  if (verbose) {
    n_pathways <- length(unique(network$source))
    n_genes <- length(unique(network$target))
    n_interactions <- nrow(network)
    cat("   Pathways:", n_pathways, "\n")
    cat("   Genes:", n_genes, "\n")
    cat("   Pathway-gene links:", n_interactions, "\n")

    # List pathways
    pathways <- sort(unique(network$source))
    cat("   Available pathways:", paste(pathways, collapse = ", "), "\n")
  }

  # Add metadata as attributes
  attr(network, "database") <- "progeny"
  attr(network, "organism") <- organism
  attr(network, "top") <- top
  attr(network, "network_type") <- "pathway"

  return(network)
}


# ==============================================================================
# PROTEIN-PROTEIN INTERACTION NETWORKS
# ==============================================================================

#' Get protein-protein interaction network from OmniPath
#'
#' Retrieves PPI data from OmniPath (aggregating STRING, BioGRID, IntAct,
#' and other sources) for a set of genes and returns an igraph network object.
#'
#' @param genes Character vector of gene symbols to include in the network.
#' @param organism Integer. NCBI taxonomy ID. Default: 9606 (human).
#'   Common values: 9606 (human), 10090 (mouse), 10116 (rat).
#' @param resources Character vector. Specific interaction databases to query.
#'   Default: NULL (all available). Examples: "STRING", "BioGRID", "IntAct",
#'   "SIGNOR", "PhosphoSite".
#' @param min_resources Integer. Minimum number of databases supporting an
#'   interaction for it to be included. Default: 1. Higher values = more
#'   stringent filtering.
#' @param drop_isolates Logical. Remove query genes with zero interactions
#'   (degree 0) from the returned graph, instead of keeping them as
#'   disconnected nodes. Default: FALSE.
#' @param directed Logical. Whether to return directed interactions.
#'   Default: FALSE (undirected PPI network).
#' @param verbose Logical. Print network summary. Default: TRUE.
#'
#' @return An igraph graph object with:
#'   \describe{
#'     \item{Vertex attributes}{name (gene symbol), degree, betweenness}
#'     \item{Edge attributes}{source, target, n_resources, resources (database names)}
#'     \item{Graph attributes}{organism, n_query_genes, n_mapped_genes}
#'   }
#'   Also has attribute "interaction_df" with the raw interaction data.frame.
#'
#' @details
#' OmniPath aggregates protein-protein interactions from dozens of databases
#' including STRING, BioGRID, IntAct, SIGNOR, PhosphoSite, and more. This
#' provides broader coverage than any single database.
#'
#' Only interactions where BOTH endpoints are in the provided gene list are
#' returned (within-list network). To include first neighbors, add them to
#' the gene list before calling.
#'
#' Network metrics (degree, betweenness) are computed automatically and stored
#' as vertex attributes for use in downstream visualization.
#'
#' @examples
#' \dontrun{
#' # PPI network for a set of DEGs
#' graph <- APOLLO_get_ppi(c("TP53", "BRCA1", "EGFR", "MYC", "CDK2"))
#'
#' # Mouse network from specific databases
#' graph <- APOLLO_get_ppi(genes, organism = 10090, resources = c("STRING", "BioGRID"))
#'
#' # Stricter filtering: require at least 2 supporting databases
#' graph <- APOLLO_get_ppi(genes, min_resources = 2)
#'
#' }
#' @export
APOLLO_get_ppi <- function(genes,
                            organism = 9606,
                            resources = NULL,
                            min_resources = 1,
                            drop_isolates = FALSE,
                            directed = FALSE,
                            verbose = TRUE) {

  if (!requireNamespace("OmnipathR", quietly = TRUE)) {
    stop("Package 'OmnipathR' is required. Install from Bioconductor:\n",
         "  BiocManager::install('OmnipathR')", call. = FALSE)
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install with:\n",
         "  install.packages('igraph')", call. = FALSE)
  }

  genes <- unique(genes)
  n_query <- length(genes)

  if (n_query < 2) {
    stop("At least 2 genes are required to build a network.", call. = FALSE)
  }

  if (verbose) {
    cat("[APOLLO] Fetching PPI network from OmniPath...\n")
    cat("   Query genes:", n_query, "\n")
    if (!is.null(resources)) cat("   Resources:", paste(resources, collapse = ", "), "\n")
  }

  # --- Fetch interactions ---
  interactions <- tryCatch({
    args <- list(organism = organism)
    if (!is.null(resources)) args$resources <- resources
    do.call(OmnipathR::import_all_interactions, args)
  }, error = function(e) {
    stop("Failed to fetch interactions from OmniPath: ", e$message, call. = FALSE)
  })

  if (is.null(interactions) || nrow(interactions) == 0) {
    stop("No interactions returned from OmniPath.", call. = FALSE)
  }

  # --- Filter to query genes (both endpoints must be in the list) ---
  interactions <- interactions[
    interactions$source_genesymbol %in% genes &
    interactions$target_genesymbol %in% genes,
  ]

  if (nrow(interactions) == 0) {
    warning("No interactions found between the provided genes.")
    # Return empty graph with the query genes as isolated nodes
    graph <- igraph::make_empty_graph(directed = directed)
    graph <- igraph::add_vertices(graph, n_query, name = genes)
    igraph::V(graph)$degree <- 0
    igraph::V(graph)$betweenness <- 0
    igraph::graph_attr(graph, "organism") <- organism
    igraph::graph_attr(graph, "n_query_genes") <- n_query
    igraph::graph_attr(graph, "n_mapped_genes") <- 0
    attr(graph, "interaction_df") <- data.frame()
    return(graph)
  }

  # --- Filter by minimum number of supporting resources ---
  if (min_resources > 1 && "n_resources" %in% colnames(interactions)) {
    interactions <- interactions[interactions$n_resources >= min_resources, ]
    if (nrow(interactions) == 0) {
      warning("No interactions remain after min_resources filter (", min_resources, ").")
      graph <- igraph::make_empty_graph(directed = directed)
      graph <- igraph::add_vertices(graph, n_query, name = genes)
      igraph::V(graph)$degree <- 0
      igraph::V(graph)$betweenness <- 0
      return(graph)
    }
  }

  # --- Build edge data.frame ---
  edge_df <- data.frame(
    from       = interactions$source_genesymbol,
    to         = interactions$target_genesymbol,
    stringsAsFactors = FALSE
  )

  # Add available edge metadata
  if ("n_resources" %in% colnames(interactions)) {
    edge_df$n_resources <- interactions$n_resources
  }
  if ("sources" %in% colnames(interactions)) {
    edge_df$resources <- interactions$sources
  }

  # Remove self-loops
  edge_df <- edge_df[edge_df$from != edge_df$to, ]

  # Remove duplicate edges (keep the one with most resources)
  if (!directed) {
    # Normalize undirected edges (alphabetical order)
    edge_df$key <- apply(edge_df[, c("from", "to")], 1, function(x) {
      paste(sort(x), collapse = "_")
    })
    if ("n_resources" %in% colnames(edge_df)) {
      edge_df <- edge_df[order(-edge_df$n_resources), ]
    }
    edge_df <- edge_df[!duplicated(edge_df$key), ]
    edge_df$key <- NULL
  }

  # --- Build igraph ---
  graph <- igraph::graph_from_data_frame(edge_df, directed = directed, vertices = NULL)

  # Add isolated query genes as vertices
  missing_genes <- setdiff(genes, igraph::V(graph)$name)
  if (length(missing_genes) > 0) {
    graph <- igraph::add_vertices(graph, length(missing_genes), name = missing_genes)
  }

  # --- Compute network metrics ---
  igraph::V(graph)$degree <- igraph::degree(graph)
  igraph::V(graph)$betweenness <- igraph::betweenness(graph)

  # --- Graph attributes ---
  igraph::graph_attr(graph, "organism") <- organism
  igraph::graph_attr(graph, "n_query_genes") <- n_query
  igraph::graph_attr(graph, "n_mapped_genes") <- sum(igraph::V(graph)$degree > 0)

  # Store raw interaction data
  attr(graph, "interaction_df") <- edge_df

  if (verbose) {
    cat("   Nodes:", igraph::vcount(graph), "\n")
    cat("   Edges:", igraph::ecount(graph), "\n")
    n_connected <- sum(igraph::V(graph)$degree > 0)
    n_isolated <- igraph::vcount(graph) - n_connected
    cat("   Connected:", n_connected, "| Isolated:", n_isolated, "\n")
    if (igraph::ecount(graph) > 0) {
      cat("   Mean degree:", round(mean(igraph::V(graph)$degree), 1), "\n")
      top_hubs <- sort(igraph::V(graph)$degree, decreasing = TRUE)
      top_hubs <- head(top_hubs[top_hubs > 0], 5)
      cat("   Top hubs:", paste(names(top_hubs), paste0("(", top_hubs, ")"),
                               collapse = ", "), "\n")
    }
  }

  # --- Drop isolated nodes (degree 0) if requested ---
  if (drop_isolates) {
    n_isolated <- sum(igraph::V(graph)$degree == 0)
    if (n_isolated > 0) {
      graph <- igraph::delete_vertices(graph, igraph::V(graph)[igraph::V(graph)$degree == 0])
      if (verbose) cat("   Dropped", n_isolated, "isolated node(s) (degree 0)\n")
    }
  }

  return(graph)
}


#' Get protein-protein interaction network directly from STRING
#'
#' @description Queries the STRING REST API directly (rather than via
#' OmniPath) for a set of genes, returning both an igraph network object
#' (with the same shape as \code{\link{APOLLO_get_ppi}}, so it can be used
#' as a drop-in replacement) and, optionally, STRING's own rendered network
#' image. STRING's image renderer uses its own layout engine, which is
#' often better at avoiding node/label overlap than re-laying out the same
#' edges with \code{\link{AETHER_plot_ppi_network}}.
#'
#' @param genes Character vector of gene symbols to include in the network.
#' @param organism Integer. NCBI taxonomy ID. Default: 9606 (human).
#' @param score_threshold Integer (0-1000). STRING's \code{required_score}
#'   (combined confidence score cutoff). Default: 400 (STRING's own
#'   "medium confidence" default).
#' @param network_type Character. "functional" (default; includes predicted/
#'   indirect functional associations, STRING's website default) or
#'   "physical" (direct physical binding interactions only).
#' @param add_nodes Integer. Number of extra first-shell interactors STRING
#'   should add beyond the query gene list. Default: 0 (within-list network
#'   only, matching \code{APOLLO_get_ppi()}'s default behavior).
#' @param image_format Character. Which STRING-rendered image to download:
#'   "image" (default, PNG), "highres_image" (high-resolution PNG), "svg"
#'   (vector), or "none" (skip image download entirely).
#' @param image_path Character. File path to save the downloaded image.
#'   Required when \code{image_format != "none"}.
#' @param verbose Logical. Print progress and network summary. Default: TRUE.
#'
#' @return A list of class "apollo_string_ppi" with:
#'   \describe{
#'     \item{graph}{igraph object — same vertex/edge/graph attribute shape
#'       as \code{APOLLO_get_ppi()}'s return value.}
#'     \item{string_ids}{data.frame. Raw STRING ID mapping for the query genes.}
#'     \item{interaction_df}{data.frame. Edge list with STRING combined scores.}
#'     \item{image_path}{Character or NULL. Path to the downloaded image, if any.}
#'     \item{params}{List of the call parameters used.}
#'   }
#'
#' @details
#' This calls the STRING REST API (string-db.org/api) directly over HTTP
#' using POST requests (via the \code{curl} package), with identifiers sent
#' in the request body rather than a GET query string — this avoids the
#' URL-length ceiling that GET hits once gene lists run into the hundreds
#' (STRING's own documented recommendation for large identifier lists).
#' Gene symbols are first resolved to STRING IDs via the
#' \code{get_string_ids} endpoint; genes that fail to resolve are reported
#' and kept as isolated nodes in the returned graph (for parity with
#' \code{APOLLO_get_ppi()}'s \code{n_query_genes}/\code{n_mapped_genes}
#' graph attributes).
#'
#' On a non-200 STRING API response, the actual HTTP status code and
#' response body are included in the error (rather than the generic
#' connection-failure message a GET-based \code{url()} connection would
#' give), to make transient vs. persistent failures easier to tell apart.
#'
#' @examples
#' \dontrun{
#' # Basic STRING network, with STRING's own rendered image
#' res <- APOLLO_get_string_ppi(c("TP53", "BRCA1", "EGFR", "MYC", "CDK2"),
#'                               image_path = "string_network.png")
#'
#' # Data only, no image, stricter confidence threshold
#' res <- APOLLO_get_string_ppi(genes, score_threshold = 700,
#'                               image_format = "none")
#'
#' # Feed the graph into the existing ggraph-based plot for comparison
#' p <- AETHER_plot_ppi_network(res$graph, layout = "stress")
#' }
#' @export
APOLLO_get_string_ppi <- function(genes,
                                   organism = 9606,
                                   score_threshold = 400,
                                   network_type = c("functional", "physical"),
                                   add_nodes = 0,
                                   image_format = c("image", "highres_image", "svg", "none"),
                                   image_path = NULL,
                                   verbose = TRUE) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.", call. = FALSE)
  }
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop("Package 'curl' is required. Install with:\n",
         "  install.packages('curl')", call. = FALSE)
  }

  network_type <- match.arg(network_type)
  image_format <- match.arg(image_format)

  if (image_format != "none" && is.null(image_path)) {
    stop("'image_path' must be provided when image_format != 'none'.", call. = FALSE)
  }

  genes <- unique(genes)
  n_query <- length(genes)

  if (n_query < 2) {
    stop("At least 2 genes are required to build a network.", call. = FALSE)
  }

  caller_id <- "GAIA_ZERO_DAWN"
  api_base <- "https://string-db.org/api"

  if (verbose) {
    cat("[APOLLO] Fetching PPI network from STRING API...\n")
    cat("   Query genes:", n_query, "\n")
    cat("   Network type:", network_type, "| Required score:", score_threshold, "\n")
  }

  # --- Resolve gene symbols to STRING IDs ---
  id_resp <- .apollo_string_api_post(
    "tsv/get_string_ids",
    list(
      identifiers     = paste(utils::URLencode(genes, reserved = TRUE), collapse = "%0d"),
      species         = organism,
      limit           = 1,
      echo_query      = 1,
      caller_identity = caller_id
    ),
    api_base = api_base
  )
  id_map <- utils::read.delim(text = rawToChar(id_resp$content), stringsAsFactors = FALSE)

  if (is.null(id_map) || nrow(id_map) == 0 || !"stringId" %in% colnames(id_map)) {
    stop("STRING API returned no ID mappings for the provided genes.", call. = FALSE)
  }

  mapped_genes <- unique(id_map$queryItem)
  unmapped <- setdiff(genes, mapped_genes)

  if (verbose) {
    cat("   Mapped:", length(mapped_genes), "/", n_query, "genes\n")
    if (length(unmapped) > 0) {
      cat("   Unmapped (kept as isolated nodes):", paste(unmapped, collapse = ", "), "\n")
    }
  }

  # --- Fetch interaction network ---
  string_id_params <- list(
    identifiers     = paste(utils::URLencode(id_map$stringId, reserved = TRUE), collapse = "%0d"),
    species         = organism,
    required_score  = score_threshold,
    network_type    = network_type,
    add_nodes       = add_nodes,
    caller_identity = caller_id
  )

  net_resp <- .apollo_string_api_post("tsv/network", string_id_params, api_base = api_base)
  net_df <- if (length(net_resp$content) == 0) {
    data.frame()
  } else {
    utils::read.delim(text = rawToChar(net_resp$content), stringsAsFactors = FALSE)
  }

  # --- Build edge data.frame ---
  if (is.null(net_df) || nrow(net_df) == 0) {
    edge_df <- data.frame(from = character(0), to = character(0),
                           score = numeric(0), stringsAsFactors = FALSE)
  } else {
    edge_df <- data.frame(
      from  = net_df$preferredName_A,
      to    = net_df$preferredName_B,
      score = net_df$score,
      stringsAsFactors = FALSE
    )
    edge_df <- edge_df[edge_df$from != edge_df$to, ]

    # Normalize undirected edges (alphabetical order), keep highest score
    edge_df$key <- apply(edge_df[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "_"))
    edge_df <- edge_df[order(-edge_df$score), ]
    edge_df <- edge_df[!duplicated(edge_df$key), ]
    edge_df$key <- NULL
  }

  # --- Build igraph ---
  if (nrow(edge_df) > 0) {
    graph <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = NULL)
  } else {
    graph <- igraph::make_empty_graph(directed = FALSE)
  }

  missing_genes <- setdiff(genes, igraph::V(graph)$name)
  if (length(missing_genes) > 0) {
    graph <- igraph::add_vertices(graph, length(missing_genes), name = missing_genes)
  }

  igraph::V(graph)$degree <- igraph::degree(graph)
  igraph::V(graph)$betweenness <- igraph::betweenness(graph)

  igraph::graph_attr(graph, "organism") <- organism
  igraph::graph_attr(graph, "n_query_genes") <- n_query
  igraph::graph_attr(graph, "n_mapped_genes") <- sum(igraph::V(graph)$degree > 0)
  igraph::graph_attr(graph, "source") <- "STRING_API"
  igraph::graph_attr(graph, "network_type") <- network_type
  igraph::graph_attr(graph, "score_threshold") <- score_threshold

  attr(graph, "interaction_df") <- edge_df

  if (verbose) {
    cat("   Nodes:", igraph::vcount(graph), "\n")
    cat("   Edges:", igraph::ecount(graph), "\n")
    n_connected <- sum(igraph::V(graph)$degree > 0)
    n_isolated <- igraph::vcount(graph) - n_connected
    cat("   Connected:", n_connected, "| Isolated:", n_isolated, "\n")
    if (igraph::ecount(graph) > 0) {
      cat("   Mean degree:", round(mean(igraph::V(graph)$degree), 1), "\n")
      top_hubs <- sort(igraph::V(graph)$degree, decreasing = TRUE)
      top_hubs <- head(top_hubs[top_hubs > 0], 5)
      cat("   Top hubs:", paste(names(top_hubs), paste0("(", top_hubs, ")"),
                               collapse = ", "), "\n")
    }
  }

  # --- Download STRING's rendered network image ---
  if (image_format != "none") {
    tryCatch({
      img_resp <- .apollo_string_api_post(paste0(image_format, "/network"),
                                           string_id_params, api_base = api_base)
      writeBin(img_resp$content, image_path)
    }, error = function(e) {
      warning("Failed to download STRING network image: ", e$message, call. = FALSE)
      image_path <<- NULL
    })
    if (verbose && !is.null(image_path)) cat("[APOLLO] STRING network image saved:", image_path, "\n")
  } else {
    image_path <- NULL
  }

  result <- list(
    graph = graph,
    string_ids = id_map,
    interaction_df = edge_df,
    image_path = image_path,
    params = list(
      organism = organism,
      score_threshold = score_threshold,
      network_type = network_type,
      add_nodes = add_nodes,
      image_format = image_format
    )
  )
  class(result) <- "apollo_string_ppi"
  return(result)
}


#' POST a form-encoded request to the STRING API
#'
#' @description Sends \code{params} as an \code{application/x-www-form-urlencoded}
#' POST body (the same key=value pairs STRING's GET endpoints accept as query
#' parameters, just relocated to the request body). Used by
#' \code{\link{APOLLO_get_string_ppi}} instead of GET so identifier lists in
#' the hundreds/thousands don't hit the URL-length ceiling GET requests run
#' into. Values already containing percent-encoded characters (e.g. the
#' \code{identifiers} field, joined with \code{\%0d}) are passed through as-is.
#'
#' @param path Character. API path after \code{api_base/}, e.g.
#'   \code{"tsv/network"} or \code{"image/network"}.
#' @param params Named list of POST body fields (character/numeric scalars).
#' @param api_base Character. STRING API base URL.
#'
#' @return The raw response list from \code{curl::curl_fetch_memory()}
#'   (\code{$status_code}, \code{$content}, ...).
#' @keywords internal
.apollo_string_api_post <- function(path, params, api_base = "https://string-db.org/api") {
  body <- paste(paste0(names(params), "=", params), collapse = "&")
  handle <- curl::new_handle(postfields = body)

  resp <- tryCatch(
    curl::curl_fetch_memory(paste0(api_base, "/", path), handle = handle),
    error = function(e) stop("Failed to reach STRING API (", path, "): ", e$message, call. = FALSE)
  )

  if (resp$status_code != 200) {
    stop("STRING API request to '", path, "' failed with HTTP ", resp$status_code,
         ": ", rawToChar(resp$content), call. = FALSE)
  }

  resp
}


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' List available decoupleR databases
#'
#' Prints information about available prior knowledge databases
#' for activity inference.
#'
#' @return Invisibly returns a data.frame of database information.
#'
#' @examples
#' \dontrun{
#' APOLLO_list_activity_databases()
#'
#' }
#' @export
APOLLO_list_activity_databases <- function() {

  cat("[APOLLO] Available Prior Knowledge Databases n\n")

  cat(" TF-TARGET NETWORKS (for TF activity inference):\n")
  cat(" ------------------------------------------------\n")
  cat(" CollecTRI (default):\n")
  cat("   - Comprehensive TF-target collection\n")
  cat("   - Good coverage and accuracy balance\n")
  cat("   - Function: APOLLO_get_collectri()\n\n")

  cat(" DoRothEA:\n")
  cat("   - Curated TF regulons with confidence levels (A-E)\n")
  cat("   - Higher confidence = fewer but more reliable interactions\n")
  cat("   - Function: APOLLO_get_dorothea()\n\n")

  cat(" PATHWAY NETWORKS (for pathway activity inference):\n")
  cat(" --------------------------------------------------\n")
  cat(" PROGENy:\n")
  cat("   - Pathway-responsive gene signatures\n")
  cat("   - 14 cancer-relevant pathways\n")
  cat("   - Derived from perturbation experiments\n")
  cat("   - Function: APOLLO_get_progeny()\n\n")

  cat(" Note: Custom networks can also be used with ARTEMIS_run_decoupler()\n")
  cat(" Format: data.frame with 'source', 'target', and 'mor' or 'weight' columns\n")

  # Return summary table invisibly
  db_info <- data.frame(
    database = c("collectri", "dorothea", "progeny"),
    type = c("tf_target", "tf_target", "pathway"),
    function_name = c("APOLLO_get_collectri", "APOLLO_get_dorothea", "APOLLO_get_progeny"),
    organisms = c("human, mouse", "human, mouse", "human, mouse"),
    stringsAsFactors = FALSE
  )

  invisible(db_info)
}
