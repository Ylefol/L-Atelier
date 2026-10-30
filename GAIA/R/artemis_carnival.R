# GAIA/Artemis/carnival.R
# Causal network contextualization via CARNIVAL
#
# Takes per-node activity scores (TF activity from ARTEMIS_run_decoupler(),
# or a MOFA factor's per-gene signature via HEPHAESTUS_mofa_factor_signature()
# + ARTEMIS_run_decoupler()) and a signed/directed prior knowledge network
# (APOLLO_get_ppi(signed_only = TRUE)), and searches for the causal subnetwork
# that best explains the measured activities -- Standard CARNIVAL when known
# perturbation targets are supplied, Inverse CARNIVAL otherwise.
#
# CARNIVAL was not installed when this file was written, so its current
# (runVanillaCarnival/runInverseCarnival) API and I/O shapes were confirmed
# against the Bioconductor release vignette source rather than by running the
# package -- see the function-level comments below for exactly what was and
# wasn't verified. Confirm against the installed version on first real run.


# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Resolve the pkn argument to a CARNIVAL-shaped data.frame
#' @keywords internal
.artemis_carnival_resolve_pkn <- function(pkn) {
  if (inherits(pkn, "igraph")) {
    pkn_df <- attr(pkn, "carnival_pkn")
    if (is.null(pkn_df)) {
      stop("'pkn' is an igraph object but has no 'carnival_pkn' attribute -- ",
           "build it with APOLLO_get_ppi(genes, signed_only = TRUE).", call. = FALSE)
    }
    return(pkn_df)
  }

  if (is.data.frame(pkn)) {
    required_cols <- c("source", "interaction", "target")
    if (!all(required_cols %in% colnames(pkn))) {
      stop("'pkn' data.frame must have columns: ",
           paste(required_cols, collapse = ", "),
           " ('interaction' holding the sign, +1/-1) -- the shape ",
           "APOLLO_get_ppi(signed_only = TRUE) attaches as attr(graph, 'carnival_pkn').",
           call. = FALSE)
    }
    return(pkn)
  }

  stop("'pkn' must be an igraph object from APOLLO_get_ppi(signed_only = TRUE), ",
       "or a data.frame with columns source/interaction/target.", call. = FALSE)
}


#' Discretize a named activity vector to {-1, 0, 1} and drop sub-threshold entries
#' @keywords internal
.artemis_carnival_discretize <- function(x, threshold, verbose) {
  n_before <- length(x)
  keep <- abs(x) >= threshold
  x <- x[keep]
  x <- vapply(x, function(v) if (v > 0) 1 else -1, numeric(1))
  if (verbose) {
    cat("    Discretized (threshold =", threshold, "):", length(x), "/", n_before,
        "nodes kept\n")
  }
  x
}


# ==============================================================================
# CORE FUNCTION
# ==============================================================================

#' Contextualize a causal signaling network with CARNIVAL
#'
#' Given per-node activity scores (e.g. TF activity from
#' \code{ARTEMIS_run_decoupler()}) and a signed/directed prior knowledge
#' network, runs CARNIVAL to find the subnetwork of causal paths that best
#' explains the measured activities. Runs Standard CARNIVAL when
#' \code{perturbations} is supplied (known input nodes), Inverse CARNIVAL
#' otherwise (broader, weaker search over possible upstream explanations).
#'
#' @param measurements Named numeric vector of per-node activity/signature
#'   scores (e.g. one column of a \code{decoupler_result$activities} matrix).
#'   Names must match node names in \code{pkn}.
#' @param pkn Prior knowledge network: either an igraph object from
#'   \code{APOLLO_get_ppi(genes, signed_only = TRUE)}, or a data.frame with
#'   columns \code{source}, \code{interaction} (+1 activation / -1 inhibition),
#'   \code{target}.
#' @param perturbations Named numeric vector of known perturbation/input
#'   nodes (e.g. a drug target forced to -1). Default \code{NULL} runs
#'   Inverse CARNIVAL (no known inputs); supplying this runs Standard
#'   CARNIVAL.
#' @param discretize Logical. CARNIVAL's ILP formulation expects discretized
#'   node states, not continuous scores -- see \code{CARNIVAL}'s own
#'   documentation ("TF activities ... subsequently discretized"). Default
#'   \code{TRUE} applies \code{discretize_threshold} and \code{sign()}.
#'   Set \code{FALSE} only if \code{measurements}/\code{perturbations} are
#'   already in \{-1, 0, 1\}.
#' @param discretize_threshold Numeric. Minimum \code{abs(score)} to keep a
#'   node when \code{discretize = TRUE}; nodes below this are dropped rather
#'   than fixed at 0. Default \code{1}, chosen as a generic
#'   "one-score-unit" cutoff with no reference-backed justification --
#'   this is NOT an established threshold for TF-activity scores from any
#'   particular decoupleR method. Inspect the distribution of
#'   \code{measurements} and adjust for your data before trusting results
#'   downstream of this step.
#' @param solver Character. CARNIVAL ILP solver. Only \code{"lpsolve"} (the
#'   default) is confirmed against the current CARNIVAL vignette -- it's
#'   pure-R (CRAN \code{lpSolve}), needs no external binary, but only
#'   practical for small networks. \code{"cbc"}/\code{"cplex"}/\code{"gurobi"}
#'   are dispatched by naming convention
#'   (\code{defaultCbcCarnivalOptions()} etc.) and need that solver's binary
#'   installed separately -- not individually confirmed here; if dispatch
#'   fails, build \code{carnivalOptions} manually per
#'   \code{?CARNIVAL::runVanillaCarnival}.
#' @param output_dir Character. Directory CARNIVAL writes its working files
#'   to (\code{carnivalOptions$outputFolder}). Default \code{NULL} uses a
#'   fresh \code{tempfile()} directory.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{CARNIVAL::runVanillaCarnival()}
#'   / \code{CARNIVAL::runInverseCarnival()}.
#'
#' @return An \code{artemis_carnival} S3 object:
#'   \item{edges}{data.frame: source, sign, target, weight}
#'   \item{nodes}{data.frame: node, zero_act, up_act, down_act, avg_act}
#'   \item{graph}{igraph object built from \code{edges}, with \code{sign} and
#'     \code{carnival_weight} edge attributes (deliberately NOT named
#'     "weight" -- see source comment) and \code{avg_activity} vertex
#'     attribute -- pass directly to \code{AETHER_plot_ppi_network()}}
#'   \item{mode}{\code{"standard"} or \code{"inverse"}}
#'   \item{n_measurements_used}{how many measurement nodes overlapped \code{pkn}
#'     and survived discretization}
#'   \item{params}{the call's key arguments, for reference}
#'   \item{raw}{the unmodified list returned by CARNIVAL, for anything not
#'     surfaced above}
#'
#' @details
#' Node names in \code{measurements}/\code{perturbations} that aren't present
#' in \code{pkn} are dropped (reported when \code{verbose = TRUE}) -- CARNIVAL
#' requires every measured/perturbed node to exist in the network.
#'
#' @examples
#' \dontrun{
#' # TF activity from a DE contrast, standalone (no MOFA/COSMOS involved)
#' tf <- ARTEMIS_infer_tf_activity(contrast_mat, method = "ulm")
#' graph <- APOLLO_get_ppi(rownames(tf$activities), signed_only = TRUE)
#' result <- ARTEMIS_run_carnival(
#'   measurements = tf$activities[, "contrast1"],
#'   pkn = graph
#' )
#' AETHER_plot_ppi_network(result$graph, color_by = "avg_activity" )
#'
#' # TF activity from a MOFA factor's loadings
#' sig <- HEPHAESTUS_mofa_factor_signature(mofa_result, factor = "Factor1", view = "rna")
#' tf_factor <- ARTEMIS_run_decoupler(sig, APOLLO_get_collectri(), method = "ulm")
#' result <- ARTEMIS_run_carnival(
#'   measurements = tf_factor$activities[, 1],
#'   pkn = graph
#' )
#'
#' }
#' @export
ARTEMIS_run_carnival <- function(measurements,
                                  pkn,
                                  perturbations = NULL,
                                  discretize = TRUE,
                                  discretize_threshold = 1,
                                  solver = "lpsolve",
                                  output_dir = NULL,
                                  verbose = TRUE,
                                  ...) {

  if (!requireNamespace("CARNIVAL", quietly = TRUE)) {
    stop("Package 'CARNIVAL' is required. Install from Bioconductor:\n",
         "  BiocManager::install('CARNIVAL')\n",
         "and, for the default solver, install CRAN 'lpSolve'.", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.", call. = FALSE)
  }

  if (is.null(names(measurements)) || any(names(measurements) == "")) {
    stop("'measurements' must be a named numeric vector (node name -> activity score).",
         call. = FALSE)
  }
  if (!is.null(perturbations) &&
      (is.null(names(perturbations)) || any(names(perturbations) == ""))) {
    stop("'perturbations' must be a named numeric vector (node name -> state).", call. = FALSE)
  }

  pkn_df <- .artemis_carnival_resolve_pkn(pkn)
  mode <- if (is.null(perturbations)) "inverse" else "standard"

  if (verbose) {
    cat("[ARTEMIS] CARNIVAL causal network contextualization\n")
    cat("    Mode:", mode, if (mode == "inverse") "(no known perturbation targets)" else "", "\n")
    cat("    PKN:", nrow(pkn_df), "edges,",
        length(unique(c(pkn_df$source, pkn_df$target))), "nodes\n")
    cat("    Measurements supplied:", length(measurements), "\n")
  }

  if (discretize) {
    measurements <- .artemis_carnival_discretize(measurements, discretize_threshold, verbose)
    if (!is.null(perturbations)) {
      perturbations <- .artemis_carnival_discretize(perturbations, discretize_threshold, verbose = FALSE)
    }
  }

  pkn_nodes <- unique(c(pkn_df$source, pkn_df$target))

  n_meas_before <- length(measurements)
  measurements <- measurements[names(measurements) %in% pkn_nodes]
  if (verbose) {
    cat("    Measurements in PKN:", length(measurements), "/", n_meas_before, "\n")
  }
  if (length(measurements) == 0) {
    stop("No measurement node names overlap the PKN -- check identifier types match ",
         "(e.g. both gene symbols) between 'measurements' and 'pkn'.", call. = FALSE)
  }

  if (!is.null(perturbations)) {
    n_pert_before <- length(perturbations)
    perturbations <- perturbations[names(perturbations) %in% pkn_nodes]
    if (verbose) cat("    Perturbations in PKN:", length(perturbations), "/", n_pert_before, "\n")
    if (length(perturbations) == 0) {
      stop("No perturbation node names overlap the PKN.", call. = FALSE)
    }
  }

  # --- carnivalOptions ---
  solver <- tolower(solver)
  opts_fn_name <- switch(solver,
    lpsolve = "defaultLpSolveCarnivalOptions",
    cbc     = "defaultCbcCarnivalOptions",
    cplex   = "defaultCplexCarnivalOptions",
    gurobi  = "defaultGurobiCarnivalOptions",
    stop("Unknown solver '", solver, "'. Options: lpsolve, cbc, cplex, gurobi.", call. = FALSE)
  )
  if (!exists(opts_fn_name, where = asNamespace("CARNIVAL"), inherits = FALSE)) {
    stop("CARNIVAL's installed version has no '", opts_fn_name, "' function -- only ",
         "'lpsolve' was confirmed against the current Bioconductor vignette. Build ",
         "carnivalOptions manually and call CARNIVAL::runVanillaCarnival()/",
         "runInverseCarnival() directly, or check ?CARNIVAL for the current solver ",
         "option helpers.", call. = FALSE)
  }
  carnival_options <- do.call(getExportedValue("CARNIVAL", opts_fn_name), list())

  if (is.null(output_dir)) {
    output_dir <- tempfile("carnival_")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  carnival_options$outputFolder <- output_dir
  if (verbose) cat("    Solver:", solver, "| outputFolder:", output_dir, "\n")

  if (verbose) cat("[ARTEMIS] Running CARNIVAL (this can take a while)...\n")

  raw <- if (mode == "standard") {
    CARNIVAL::runVanillaCarnival(
      perturbations = perturbations,
      measurements = measurements,
      priorKnowledgeNetwork = pkn_df,
      carnivalOptions = carnival_options,
      ...
    )
  } else {
    CARNIVAL::runInverseCarnival(
      measurements = measurements,
      priorKnowledgeNetwork = pkn_df,
      carnivalOptions = carnival_options,
      ...
    )
  }

  # --- Parse result ---
  # weightedSIF (Node1/Sign/Node2/Weight) and nodesAttributes (Node/ZeroAct/
  # UpAct/DownAct/AvgAct) are CARNIVAL's long-standing output shape across
  # past versions, but NOT individually re-confirmed against the current
  # runVanillaCarnival()/runInverseCarnival() release -- this check fails
  # loudly rather than silently mis-extracting if that's changed.
  expected <- c("weightedSIF", "nodesAttributes")
  if (!all(expected %in% names(raw))) {
    stop("Unexpected CARNIVAL result structure -- expected elements '",
         paste(expected, collapse = "', '"), "' but the result has: ",
         paste(names(raw), collapse = ", "),
         ". Inspect the installed CARNIVAL version's return shape directly ",
         "(?CARNIVAL::runVanillaCarnival) and update this parsing step.", call. = FALSE)
  }

  sif <- raw$weightedSIF
  expected_sif_cols <- c("Node1", "Sign", "Node2", "Weight")
  if (!all(expected_sif_cols %in% colnames(sif))) {
    stop("Unexpected weightedSIF columns -- expected ",
         paste(expected_sif_cols, collapse = ", "), " but got: ",
         paste(colnames(sif), collapse = ", "), call. = FALSE)
  }
  edges <- data.frame(
    source = sif$Node1,
    sign   = sif$Sign,
    target = sif$Node2,
    weight = sif$Weight,
    stringsAsFactors = FALSE
  )

  node_attr <- raw$nodesAttributes
  node_name_col <- intersect(c("Node", "node"), colnames(node_attr))[1]
  if (is.na(node_name_col)) {
    stop("Unexpected nodesAttributes columns -- expected a 'Node' column but got: ",
         paste(colnames(node_attr), collapse = ", "), call. = FALSE)
  }
  nodes <- data.frame(
    node     = node_attr[[node_name_col]],
    zero_act = node_attr$ZeroAct,
    up_act   = node_attr$UpAct,
    down_act = node_attr$DownAct,
    avg_act  = node_attr$AvgAct,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("[ARTEMIS] Done. Network:", nrow(edges), "edges,", nrow(nodes), "nodes\n")
  }

  # --- Build plottable igraph (drop-in for AETHER_plot_ppi_network()) ---
  # NOT named "weight" -- igraph/ggraph's default layouts (e.g. "fr") auto-
  # detect an edge attribute literally called "weight" and use it to drive
  # the layout, which requires strictly positive values. CARNIVAL's weight
  # is an edge confidence/frequency score (0-100, 0 is valid) and isn't
  # meant to be a layout weight -- confirmed this breaks AETHER_plot_ppi_network()'s
  # default "fr" layout with "Weights must be positive" when any edge is 0.
  graph <- igraph::graph_from_data_frame(
    data.frame(from = edges$source, to = edges$target,
               sign = edges$sign, carnival_weight = edges$weight, stringsAsFactors = FALSE),
    directed = TRUE
  )
  igraph::V(graph)$avg_activity <- nodes$avg_act[match(igraph::V(graph)$name, nodes$node)]

  result <- list(
    edges = edges,
    nodes = nodes,
    graph = graph,
    mode = mode,
    n_measurements_used = length(measurements),
    params = list(
      discretize = discretize,
      discretize_threshold = discretize_threshold,
      solver = solver,
      output_dir = output_dir
    ),
    raw = raw
  )
  class(result) <- c("artemis_carnival", "list")
  result
}


# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' @method print artemis_carnival
#' @export
print.artemis_carnival <- function(x, ...) {
  cat("CARNIVAL Causal Network Result\n")
  cat("-------------------------------\n")
  cat("Mode:", x$mode, "\n")
  cat("Measurement nodes used:", x$n_measurements_used, "\n")
  cat("Network: ", nrow(x$edges), "edges,", nrow(x$nodes), "nodes\n")
  cat("Solver:", x$params$solver, "\n")
  cat("\nUse $edges/$nodes for tables, $graph for AETHER_plot_ppi_network().\n")
}
