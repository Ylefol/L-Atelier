###############################################################################
############### MOFA2 - Multi-Omics Factor Analysis Integration ###############
###############################################################################
# Unlike ConvCCA/RelPMDCCA in this module, MOFA2 natively tolerates samples
# missing entirely from one or more views (partial sample overlap across
# omics), rather than requiring every view to share the exact same sample
# set. It also manages its own basilisk-provisioned Python environment for
# its mofapy2 backend internally (MOFA2::run_mofa(use_basilisk = TRUE)) --
# see GAIA/R/basilisk.R for why no .gaia_mofa_env is defined there.


# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Validate list-of-matrices input for HEPHAESTUS_run_mofa()
#' @keywords internal
.hephaestus_mofa_validate <- function(X) {
    if (!is.list(X) || is.null(names(X)) || any(names(X) == "")) {
        stop("X must be a named list of matrices, one per view -- names become ",
             "MOFA2 view names.", call. = FALSE)
    }

    for (view in names(X)) {
        mat <- X[[view]]
        if (!is.matrix(mat) && !is.data.frame(mat)) {
            stop("View '", view, "' must be a matrix or data.frame.", call. = FALSE)
        }
        if (is.null(rownames(mat)) || is.null(colnames(mat))) {
            stop("View '", view, "' must have both rownames (samples) and ",
                 "colnames (features) set.", call. = FALSE)
        }
    }

    # Unlike POSEIDON_validate_multi_dataset(), views are NOT required to share
    # the same samples/order -- that's the whole point of using MOFA2 here.
    # We only check that at least one pair of views has some sample overlap,
    # since zero overlap anywhere means there is nothing to integrate.
    if (length(X) > 1) {
        pairs <- utils::combn(names(X), 2, simplify = FALSE)
        any_overlap <- any(vapply(pairs, function(p) {
            length(intersect(rownames(X[[p[1]]]), rownames(X[[p[2]]]))) > 0
        }, logical(1)))
        if (!any_overlap) {
            stop("No sample-ID overlap between any pair of views -- check that ",
                 "rownames() match (at least partially) across X.", call. = FALSE)
        }
    }

    invisible(TRUE)
}

#' Convert GAIA's samples x features list-of-matrices into MOFA2's long format
#'
#' @description MOFA2's list-of-matrices input path
#' (\code{create_mofa_from_matrix}) requires identical, identically-ordered
#' colnames (samples) across every view and does not pad partial overlap.
#' The long-format path (\code{sample, feature, view, value[, group]}) has no
#' such requirement -- rows for missing sample/view combinations are simply
#' omitted, which is exactly what's needed for GAIA's list-of-matrices inputs
#' where views may only partially overlap in samples.
#'
#' @keywords internal
.hephaestus_mofa_to_long <- function(X, groups = NULL) {
    long_list <- lapply(names(X), function(view) {
        mat <- as.matrix(X[[view]])
        df <- data.frame(
            sample  = rownames(mat)[row(mat)],
            feature = colnames(mat)[col(mat)],
            view    = view,
            value   = as.vector(mat),
            stringsAsFactors = FALSE
        )
        df[!is.na(df$value), ]
    })
    long_df <- do.call(rbind, long_list)

    if (!is.null(groups)) {
        long_df$group <- groups[long_df$sample]
        if (any(is.na(long_df$group))) {
            stop("groups is missing an entry for at least one sample present in ",
                 "X -- every sample across all views must have a group label.",
                 call. = FALSE)
        }
    }

    long_df
}

#' Build a hephaestus_mofa result object from a trained/loaded MOFA model
#'
#' @details MOFA2::get_factors()/get_weights()/calculate_variance_explained()
#' shapes used below (list-per-group of sample x factor / feature x factor /
#' view x factor matrices) match MOFA2's own documentation and vignettes as
#' of this writing -- not independently verified against every MOFA2 release,
#' so confirm dimnames/orientation on first real run if a future MOFA2
#' version changes them.
#'
#' @keywords internal
.hephaestus_mofa_build_result <- function(trained, params = NULL) {
    factors_list <- MOFA2::get_factors(trained)
    factors_df <- do.call(rbind, lapply(names(factors_list), function(grp) {
        mat <- factors_list[[grp]]
        data.frame(
            sample = rep(rownames(mat), times = ncol(mat)),
            group  = grp,
            factor = rep(colnames(mat), each = nrow(mat)),
            value  = as.vector(mat),
            stringsAsFactors = FALSE
        )
    }))

    weights_list <- MOFA2::get_weights(trained)
    weights_df <- do.call(rbind, lapply(names(weights_list), function(view) {
        mat <- weights_list[[view]]
        data.frame(
            view    = view,
            feature = rep(rownames(mat), times = ncol(mat)),
            factor  = rep(colnames(mat), each = nrow(mat)),
            value   = as.vector(mat),
            stringsAsFactors = FALSE
        )
    }))

    # calculate_variance_explained()'s r2_per_factor matrices are factors x
    # views (rownames = factors, colnames = views) -- the opposite orientation
    # from get_weights()/get_factors() (features/samples x factors) used
    # above, so view/factor are NOT extracted the same way here.
    ve_list <- MOFA2::calculate_variance_explained(trained)$r2_per_factor
    variance_df <- do.call(rbind, lapply(names(ve_list), function(grp) {
        mat <- ve_list[[grp]]
        data.frame(
            group               = grp,
            factor              = rep(rownames(mat), times = ncol(mat)),
            view                = rep(colnames(mat), each = nrow(mat)),
            variance_explained  = as.vector(mat),
            stringsAsFactors = FALSE
        )
    }))

    result <- list(
        model               = trained,
        factors             = factors_df,
        weights             = weights_df,
        variance_explained  = variance_df,
        params              = params
    )
    class(result) <- c("hephaestus_mofa", "list")
    result
}


# ==============================================================================
# TRAIN A MOFA2 MODEL
# ==============================================================================

#' Multi-Omics Factor Analysis via MOFA2
#'
#' @description Trains a MOFA2 model on multiple omics views, following the
#' same list-of-matrices convention as \code{HEPHAESTUS_multi_convCCA()} --
#' but unlike the sparse-CCA methods in this module, MOFA2 natively tolerates
#' samples missing entirely from one or more views. Internally the views are
#' converted to MOFA2's long-format input (\code{sample, feature, view,
#' value}) rather than its stricter matrix-input path (which requires
#' identical, identically-ordered sample sets across every view), so callers
#' can pass views with only partial sample overlap without manually
#' aligning/padding them first.
#'
#' @param X Named list of matrices, one per view, samples in rows and
#'   features in columns (GAIA's standard convention, same as
#'   \code{HEPHAESTUS_multi_convCCA()}). List names become MOFA2 view names.
#'   Views do not need the same samples or sample order -- rownames are the
#'   join key across views.
#' @param groups Optional named character vector or factor, sample ID ->
#'   group label, for MOFA2's multi-group framework. Default \code{NULL}
#'   (single group). MOFA2's own documentation describes multi-group as an
#'   advanced option not recommended as a first approach -- treat this as
#'   experimental.
#' @param num_factors Integer. Number of latent factors to fit. Default
#'   \code{10} (MOFA2's own package default). For views with few samples,
#'   check \code{variance_explained} in the result rather than assuming 10 is
#'   appropriate -- this default is not tuned to any particular dataset.
#' @param likelihoods Optional named character vector, view name ->
#'   \code{"gaussian"}, \code{"poisson"}, or \code{"bernoulli"}. Default
#'   \code{NULL} uses \code{"gaussian"} for every view. MOFA2's own guidance
#'   is to prefer transforming count/binary data to continuous values rather
#'   than relying on non-gaussian likelihoods, which give non-optimal results.
#' @param spikeslab_weights Logical. Use a spike-and-slab sparsity prior on
#'   factor weights. Default \code{TRUE} (MOFA2 default).
#' @param spikeslab_factors Logical. Use a spike-and-slab sparsity prior on
#'   factors. Default \code{FALSE} (MOFA2 default).
#' @param scale_views Logical. Scale each view to unit variance before
#'   training. Default \code{FALSE}.
#' @param convergence_mode Character. \code{"fast"}, \code{"medium"}, or
#'   \code{"slow"} -- controls the ELBO-change convergence threshold. Default
#'   \code{"fast"}.
#' @param maxiter Integer. Maximum training iterations. Default \code{1000}.
#' @param drop_factor_threshold Numeric. Fraction of variance explained below
#'   which a factor is dropped during training; \code{-1} disables dropping.
#'   Default \code{-1} (MOFA2 default).
#' @param seed Integer. Random seed for reproducibility. Default \code{42}
#'   (MOFA2 default).
#' @param use_basilisk Logical. Let MOFA2 provision its own basilisk-managed
#'   Python environment for its \code{mofapy2} backend. Default \code{TRUE}
#'   -- unlike GAIA's other basilisk-wrapped tools (CPAT, SignalP, Pfam,
#'   HOMER), no \code{.gaia_mofa_env} exists in \code{GAIA/R/basilisk.R};
#'   MOFA2 manages this itself. Set \code{FALSE} only if you have manually
#'   configured a Python/\code{mofapy2} installation via
#'   \code{reticulate::use_python()}.
#' @param outfile Character or \code{NULL}. Path to save the trained model
#'   (\code{.hdf5}). Default \code{NULL} uses a temp file -- pass an explicit
#'   path if you intend to reload the model later with
#'   \code{HEPHAESTUS_load_mofa_model()}, since temp files are not guaranteed
#'   to persist across R sessions.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A classed \code{hephaestus_mofa} list:
#'   \item{model}{The raw trained \code{MOFA} object, for advanced use.}
#'   \item{factors}{Tidy data.frame of factor scores (sample, group, factor,
#'     value).}
#'   \item{weights}{Tidy data.frame of feature weights (view, feature,
#'     factor, value).}
#'   \item{variance_explained}{Tidy data.frame of variance explained
#'     (group, view, factor, variance_explained).}
#'   \item{params}{List of the training/model/data options actually used.}
#'
#' @details
#' Field names for \code{model_options}/\code{data_options}/
#' \code{training_options} are set from MOFA2's own
#' \code{get_default_model_options()}/\code{get_default_data_options()}/
#' \code{get_default_training_options()} and were not independently pinned
#' against a citable spec beyond MOFA2's own documentation -- if a field name
#' changes in a future MOFA2 release, the option-setting step below will fail
#' loudly (unknown list element) rather than silently doing the wrong thing.
#'
#' This function is intentionally not wired into \code{MINERVA_scca_fit()}/
#' \code{MINERVA_cv_scca()}: those normalize every method to a
#' \code{list(W = list(...))} canonical-weight-vector contract and search
#' sparsity parameters via k-fold CV, neither of which applies to a shared
#' latent-factor model selected by number-of-factors/ELBO.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' # view_a and view_b share all 50 samples; view_c only covers a 15-sample
#' # subset -- MOFA2 handles this partial overlap natively.
#' view_a <- matrix(rnorm(50 * 100), 50, 100,
#'                   dimnames = list(paste0("S", 1:50), paste0("gA", 1:100)))
#' view_b <- matrix(rnorm(50 * 80), 50, 80,
#'                   dimnames = list(paste0("S", 1:50), paste0("gB", 1:80)))
#' view_c <- matrix(rnorm(15 * 60), 15, 60,
#'                   dimnames = list(paste0("S", 1:15), paste0("gC", 1:60)))
#'
#' result <- HEPHAESTUS_run_mofa(
#'   X = list(rna = view_a, atac = view_b, cutntag = view_c),
#'   num_factors = 5
#' )
#' print(result)
#' head(result$factors)
#'
#' }
HEPHAESTUS_run_mofa <- function(X,
                                 groups = NULL,
                                 num_factors = 10,
                                 likelihoods = NULL,
                                 spikeslab_weights = TRUE,
                                 spikeslab_factors = FALSE,
                                 scale_views = FALSE,
                                 convergence_mode = "fast",
                                 maxiter = 1000,
                                 drop_factor_threshold = -1,
                                 seed = 42,
                                 use_basilisk = TRUE,
                                 outfile = NULL,
                                 verbose = TRUE) {

    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        stop("Package 'MOFA2' is required for HEPHAESTUS_run_mofa(). Install with: ",
             "BiocManager::install('MOFA2')", call. = FALSE)
    }

    .hephaestus_mofa_validate(X)

    if (!is.null(likelihoods) &&
        (is.null(names(likelihoods)) || !all(names(X) %in% names(likelihoods)))) {
        stop("likelihoods must be a named character vector covering every view in X.",
             call. = FALSE)
    }

    if (verbose) {
        cat("[HEPHAESTUS] MOFA2 Multi-Omics Factor Analysis\n")
        cat("    Views: ", paste(names(X), collapse = ", "), "\n", sep = "")
        for (view in names(X)) {
            cat("        ", view, ": ", nrow(X[[view]]), " samples x ",
                ncol(X[[view]]), " features\n", sep = "")
        }
        if (!is.null(groups)) {
            cat("    Groups: ", length(unique(groups)), " (multi-group, experimental)\n", sep = "")
        }
        cat("    Factors requested: ", num_factors, "\n\n", sep = "")
    }

    long_df  <- .hephaestus_mofa_to_long(X, groups = groups)
    mofa_obj <- MOFA2::create_mofa(data = long_df)

    data_opts  <- MOFA2::get_default_data_options(mofa_obj)
    model_opts <- MOFA2::get_default_model_options(mofa_obj)
    train_opts <- MOFA2::get_default_training_options(mofa_obj)

    data_opts$scale_views <- scale_views

    model_opts$num_factors       <- num_factors
    model_opts$spikeslab_weights <- spikeslab_weights
    model_opts$spikeslab_factors <- spikeslab_factors
    if (!is.null(likelihoods)) {
        model_opts$likelihoods[names(likelihoods)] <- likelihoods
    }

    train_opts$convergence_mode      <- convergence_mode
    train_opts$maxiter               <- maxiter
    train_opts$drop_factor_threshold <- drop_factor_threshold
    train_opts$seed                  <- seed
    train_opts$verbose               <- verbose

    mofa_obj <- MOFA2::prepare_mofa(mofa_obj, data_options = data_opts,
                                     model_options = model_opts,
                                     training_options = train_opts)

    if (is.null(outfile)) outfile <- tempfile(fileext = ".hdf5")

    if (verbose) cat("[HEPHAESTUS] Training MOFA2 model (this may take a while)...\n")

    start_time <- proc.time()
    trained <- MOFA2::run_mofa(mofa_obj, outfile = outfile, use_basilisk = use_basilisk)
    elapsed <- (proc.time() - start_time)["elapsed"]

    if (verbose) {
        cat("    Training completed in ", round(elapsed, 1), " seconds. Model saved to: ",
            outfile, "\n\n", sep = "")
    }

    result <- .hephaestus_mofa_build_result(
        trained,
        params = list(num_factors = num_factors, likelihoods = model_opts$likelihoods,
                      spikeslab_weights = spikeslab_weights, spikeslab_factors = spikeslab_factors,
                      scale_views = scale_views, convergence_mode = convergence_mode,
                      maxiter = maxiter, drop_factor_threshold = drop_factor_threshold,
                      seed = seed, outfile = outfile)
    )

    if (verbose) print(result)

    result
}


#' Print method for hephaestus_mofa objects
#' @param x A hephaestus_mofa object
#' @param ... Additional arguments (unused)
#' @method print hephaestus_mofa
#' @export
print.hephaestus_mofa <- function(x, ...) {
    cat("MOFA2 Multi-Omics Factor Analysis\n")
    cat("------------------------------\n")
    views <- unique(x$weights$view)
    cat("Views:", paste(views, collapse = ", "), "\n")
    cat("Samples:", length(unique(x$factors$sample)), "\n")
    cat("Factors:", length(unique(x$factors$factor)), "\n")
    cat("\nTotal variance explained per view:\n")
    total_ve <- stats::aggregate(variance_explained ~ view, data = x$variance_explained, FUN = sum)
    total_ve <- total_ve[order(-total_ve$variance_explained), ]
    print(total_ve, row.names = FALSE)
}


# ==============================================================================
# LOAD A TRAINED MOFA2 MODEL
# ==============================================================================

#' Load a Previously Trained MOFA2 Model
#'
#' @description Reloads a MOFA2 model saved by \code{HEPHAESTUS_run_mofa()}
#' (or by \code{MOFA2::run_mofa()} directly) from its \code{.hdf5} file,
#' without retraining. Training is slow and stochastic, so keeping a loader
#' separate from the run step matters in practice -- mirrors the run/parse
#' split already used for HOMER (\code{.homer_run()} vs
#' \code{APOLLO_load_homer_results()}).
#'
#' @param file Character. Path to a trained MOFA2 \code{.hdf5} model file.
#'
#' @return A classed \code{hephaestus_mofa} list, identical in shape to the
#'   return value of \code{HEPHAESTUS_run_mofa()} -- \code{params} will be
#'   \code{NULL} since training options aren't recoverable from the saved
#'   model file alone.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result   <- HEPHAESTUS_run_mofa(X = list(...), outfile = "my_model.hdf5")
#' reloaded <- HEPHAESTUS_load_mofa_model("my_model.hdf5")
#'
#' }
HEPHAESTUS_load_mofa_model <- function(file) {
    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        stop("Package 'MOFA2' is required for HEPHAESTUS_load_mofa_model(). Install with: ",
             "BiocManager::install('MOFA2')", call. = FALSE)
    }

    if (!file.exists(file)) stop("file not found: ", file, call. = FALSE)

    trained <- MOFA2::load_model(file)
    .hephaestus_mofa_build_result(trained, params = NULL)
}


# ==============================================================================
# PROJECT NEW SAMPLES ONTO A TRAINED MOFA2 MODEL (OUT-OF-SAMPLE)
# ==============================================================================

#' Project New Samples onto a Trained MOFA2 Model
#'
#' @description MOFA2 has no out-of-sample projection of its own:
#' \code{predict()}/\code{impute()} only reconstruct the model's own training
#' data via its already-fitted factor scores, and
#' \code{interpolate_factors()} is MEFISTO-only (interpolates along a trained
#' covariate for existing training samples, not new ones). This fills that
#' gap: holding the trained per-view feature weights fixed, it solves for
#' each new sample's factor vector by ordinary least squares against
#' \code{Y_view = Z \%*\% t(W_view)}, after centering the new sample's data
#' the same way MOFA centered the training data internally (per-feature
#' training mean; confirmed empirically -- reconstructing training data as
#' \code{Z \%*\% t(W)} without adding the feature mean back leaves residuals
#' on the order of the feature means themselves, and centered residuals are
#' small). A view a sample has no data for simply drops out of that
#' sample's fit -- no imputation is attempted.
#'
#' @param mofa_result A \code{hephaestus_mofa} object from
#'   \code{HEPHAESTUS_run_mofa()} or \code{HEPHAESTUS_load_mofa_model()}.
#' @param train_data The named list of view matrices (samples x features)
#'   originally passed as \code{X} to \code{HEPHAESTUS_run_mofa()} to train
#'   this model. Used only to recover each view's per-feature training mean
#'   for centering -- MOFA2 does not expose this on the fitted model object.
#' @param new_data Named list of view matrices (samples x features) for the
#'   new samples to project. Must use the SAME normalization/transform as
#'   \code{train_data} (e.g. the same VST dispersion fit, the same batch
#'   correction) -- this function does not check for that, only for shared
#'   feature names. Feature sets are intersected against \code{train_data}
#'   and the model's weights per view; extra features are silently dropped.
#'   Views may cover different, partially-overlapping sets of samples, same
#'   as MOFA2 tolerates for training.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A matrix (new samples x factors) of projected factor scores, in
#'   the same units as \code{mofa_result$factors}.
#'
#' @details
#' Ordinary least squares, solved independently per sample via
#' \code{qr.solve()} on that sample's stacked available-view design. No
#' regularization is applied -- if a sample's total usable feature count
#' (summed across its available views) is small relative to the number of
#' factors, the fit is underdetermined/unstable; this is a known limitation,
#' not handled automatically. Samples with zero usable features in every
#' view are returned as all-NA rows, with a warning.
#'
#' @export
HEPHAESTUS_project_mofa_samples <- function(mofa_result, train_data, new_data,
                                             verbose = TRUE) {

    if (!inherits(mofa_result, "hephaestus_mofa")) {
        stop("mofa_result must be a hephaestus_mofa object from HEPHAESTUS_run_mofa() ",
             "or HEPHAESTUS_load_mofa_model().", call. = FALSE)
    }
    if (is.null(mofa_result$model)) {
        stop("mofa_result$model is NULL -- the raw MOFA2 model object is required.",
             call. = FALSE)
    }
    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        stop("Package 'MOFA2' is required.", call. = FALSE)
    }

    common_views <- intersect(names(train_data), names(new_data))
    if (length(common_views) == 0) {
        stop("No view names in common between train_data and new_data.", call. = FALSE)
    }

    weights_list <- MOFA2::get_weights(mofa_result$model)
    factor_names <- colnames(weights_list[[1]])
    n_factors    <- length(factor_names)

    if (verbose) {
        cat("[HEPHAESTUS] MOFA2 out-of-sample projection\n")
        cat("    Views used: ", paste(common_views, collapse = ", "), "\n", sep = "")
        cat("    Factors: ", n_factors, "\n", sep = "")
    }

    # Per-view: restrict weights + new/train data to the shared feature set,
    # and compute training per-feature means for centering.
    view_prep <- lapply(common_views, function(view) {
        w <- weights_list[[view]]
        shared_feat <- Reduce(intersect, list(rownames(w),
                                               colnames(train_data[[view]]),
                                               colnames(new_data[[view]])))
        if (length(shared_feat) == 0) {
            if (verbose) cat("    ", view, ": no shared features, dropped\n", sep = "")
            return(NULL)
        }
        list(
            W       = w[shared_feat, , drop = FALSE],
            means   = colMeans(train_data[[view]][, shared_feat, drop = FALSE]),
            new_mat = new_data[[view]][, shared_feat, drop = FALSE]
        )
    })
    names(view_prep) <- common_views
    view_prep <- view_prep[!vapply(view_prep, is.null, logical(1))]

    if (length(view_prep) == 0) {
        stop("No usable (shared-feature) views remain after intersecting train_data, ",
             "new_data, and the model's weights.", call. = FALSE)
    }

    all_samples <- unique(unlist(lapply(view_prep, function(v) rownames(v$new_mat))))

    if (verbose) {
        for (view in names(view_prep)) {
            cat("    ", view, ": ", nrow(view_prep[[view]]$W), " shared features, ",
                nrow(view_prep[[view]]$new_mat), " samples\n", sep = "")
        }
        cat("    New samples (any view): ", length(all_samples), "\n\n", sep = "")
    }

    Z <- matrix(NA_real_, nrow = length(all_samples), ncol = n_factors,
                dimnames = list(all_samples, factor_names))

    n_skipped <- 0
    for (s in all_samples) {
        y_blocks <- list()
        W_blocks <- list()
        for (view in names(view_prep)) {
            vp <- view_prep[[view]]
            if (!s %in% rownames(vp$new_mat)) next
            y_blocks[[view]] <- as.numeric(vp$new_mat[s, ]) - vp$means
            W_blocks[[view]] <- vp$W
        }
        if (length(y_blocks) == 0) {
            n_skipped <- n_skipped + 1
            next
        }
        y <- unlist(y_blocks, use.names = FALSE)
        W <- do.call(rbind, W_blocks)
        Z[s, ] <- qr.solve(W, y)
    }

    if (n_skipped > 0) {
        warning(n_skipped, " sample(s) had no usable data in any shared-feature view ",
                "and were left as NA.", call. = FALSE)
    }

    if (verbose) {
        cat("[HEPHAESTUS] Projection complete. ",
            length(all_samples) - n_skipped, "/", length(all_samples),
            " samples projected.\n", sep = "")
    }

    Z
}
