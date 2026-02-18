###############################################################################
########### Data Simulation Functions ###########
###############################################################################

#' Generate Multi-Dataset Simulations with Shared Latent Structure
#'
#' @description Generates multiple synthetic datasets that share a common latent
#' variable structure. Used for testing and validating multi-omics integration
#' methods like ConvCCA and RelPMDCCA.
#'
#' @param n_datasets Integer. Number of datasets to generate (default = 3).
#' @param n_samples Integer. Number of samples in each dataset (default = 100).
#' @param n_features Integer vector. Number of features for each dataset.
#'   Must have length equal to n_datasets (default = c(500, 1000, 750)).
#' @param signal_samples Integer. Number of samples containing true signal
#'   in the shared latent variable (default = 25). Remaining samples are noise.
#' @param signal_features Integer vector. Number of features containing true
#'   signal in each dataset. Must have length equal to n_datasets
#'   (default = c(50, 50, 40)).
#' @param seed Integer. Random seed for reproducibility (default = 1234).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{datasets}{List of length n_datasets containing simulated data matrices.
#'     Each matrix has dimensions n_samples x n_features.}
#'   \item{true_canonical}{List of true canonical vectors (one per dataset).
#'     Signal features have weight 1, noise features have weight 0.}
#'   \item{shared_latent}{Matrix (n_samples x 1) containing the shared latent
#'     variable. First signal_samples rows are 1, remaining rows are 0.}
#'   \item{signal_samples}{Integer indicating number of signal samples.}
#'   \item{signal_features}{Integer vector indicating number of signal features
#'     per dataset.}
#' }
#'
#' @details
#' The simulation model follows the sparse CCA framework:
#'
#' For each dataset i:
#' \deqn{X_i = u \cdot v_i^T + E_i}
#'
#' Where:
#' \itemize{
#'   \item u is the shared latent variable (n_samples x 1)
#'   \item v_i is the canonical vector for dataset i (n_features x 1)
#'   \item E_i is Gaussian noise ~ N(0, 1)
#' }
#'
#' The true signal is sparse:
#' \itemize{
#'   \item Only the first signal_samples samples have non-zero latent values
#'   \item Only the first `signal_features[i]` features in each dataset
#'     participate in the shared structure
#' }
#'
#' This simulation allows controlled testing of:
#' \itemize{
#'   \item Feature selection accuracy (can methods recover signal features?)
#'   \item Sample selection accuracy (can methods identify signal samples?)
#'   \item Robustness to noise
#'   \item Scalability with varying dimensions
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate 3 datasets with default parameters
#' sim <- DEMETER_generate_multi_datasets()
#'
#' # Generate 2 high-dimensional datasets
#' sim <- DEMETER_generate_multi_datasets(
#'   n_datasets = 2,
#'   n_samples = 50,
#'   n_features = c(2000, 3000),
#'   signal_samples = 10,
#'   signal_features = c(100, 150),
#'   seed = 42
#' )
#'
#' # Access simulated data
#' X_1 <- sim$datasets[[1]]
#' X_2 <- sim$datasets[[2]]
#' true_weights <- sim$true_canonical[[1]]
#'
#' }
DEMETER_generate_multi_datasets <- function(n_datasets = 3,n_samples = 100,n_features = c(500, 1000, 750),signal_samples = 25,signal_features = c(50, 50, 40),seed = 1234) {

  set.seed(seed)

  # Input validation
  if (length(n_features) != n_datasets) {
    stop("n_features must have length equal to n_datasets")
  }
  if (length(signal_features) != n_datasets) {
    stop("signal_features must have length equal to n_datasets")
  }

  # Shared latent variable: first 'signal_samples' have signal, rest are noise
  u <- matrix(c(rep(1, signal_samples), rep(0, n_samples - signal_samples)), ncol = 1)

  # Generate datasets
  datasets <- list()
  canonical_vectors <- list()

  for (i in 1:n_datasets) {
    # Canonical vector: first 'signal_features[i]' are signal
    v <- matrix(c(rep(1, signal_features[i]),
                  rep(0, n_features[i] - signal_features[i])),
                ncol = 1)

    # Dataset = shared signal + noise
    X <- u %*% t(v) + matrix(rnorm(n_samples * n_features[i]),
                             ncol = n_features[i])

    datasets[[i]] <- X
    canonical_vectors[[i]] <- v
  }

  simulation_list<-list(
    datasets = datasets,
    true_canonical = canonical_vectors,
    shared_latent = u,
    signal_samples = signal_samples,
    signal_features = signal_features
  )

  cat("Dataset dimensions:\n")
  for (i in 1:length(simulation_list$datasets)) {
    cat(sprintf("  Dataset %d: %d samples x %d features (signal in first %d features)\n",
                i, nrow(simulation_list$datasets[[i]]), ncol(simulation_list$datasets[[i]]),
                simulation_list$signal_features[i]))
  }

  return(simulation_list)
}