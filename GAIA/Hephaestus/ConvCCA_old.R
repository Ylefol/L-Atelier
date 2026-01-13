###############################################################################
########### ConvCCA - Convex Sparse Canonical Correlation Analysis ###########
###############################################################################

# Source dependencies from other GAIA modules
source("../Demeter/matrix_utils.R")
source("../Minerva/penalty_functions.R")

#' Update Canonical Vector for Two-Dataset ConvCCA
#'
#' @description Updates one canonical vector while keeping the other fixed,
#' using soft-thresholding with LASSO or SCAD penalty
#'
#' @param K Correlation matrix (transformed covariance structure)
#' @param v Current canonical vector for the other dataset
#' @param tau Sparsity tuning parameter (higher = more sparse)
#' @param penalty Penalty function to use: "LASSO" (default) or "SCAD"
#'
#' @return Updated canonical vector (normalized)
#'
#' @keywords internal
updateW.convCCA <- function(K, v, tau, penalty = "LASSO"){
    # Step 1 -- Initialize
    K <- as.matrix(K)
    v <- as.numeric(v)
    vx <- K %*% v

    # Step 2 -- Normalize
    norm.vx <- as.numeric(sqrt(t(vx) %*% vx))
    if(norm.vx == 0) norm.vx <- 1
    vx <- vx / norm.vx

    # Step 3 -- Apply penalty function
    if (penalty == "SCAD"){
        u.new <- soft_threshold_scad(vx, tau)
    } else if (penalty == "LASSO"){
        u.new <- soft_threshold_lasso(vx, tau)
    } else {
        stop("Please provide the penalty function to be either LASSO or SCAD")
    }

    # Step 4 -- Normalize again
    norm.u.new <- as.numeric(sqrt(t(u.new) %*% u.new))
    if(norm.u.new == 0) norm.u.new <- 1
    u.new <- u.new / norm.u.new
    return(u.new)
}

#' Update Canonical Vector for Multi-Dataset ConvCCA
#'
#' @description Updates one canonical vector in a multi-dataset setting by
#' maximizing correlation with all other datasets
#'
#' @param K List of correlation matrices between datasets
#' @param w List of current canonical vectors for all datasets
#' @param whichDSet Index of the dataset being updated
#' @param tau Sparsity tuning parameter for this dataset
#' @param penalty Penalty function: "LASSO" (default) or "SCAD"
#'
#' @return Updated canonical vector for dataset whichDSet
#'
#' @keywords internal
updateWmulti.convCCA <- function(K, w, whichDSet, tau, penalty = "LASSO"){
    if (class(K) != "list"){
        stop("Please provide a valid K Correlation matrix. It should be a list")
    }

    # Step 1 -- Initialize
    vx <- 0
    # Take the sum of canonical pairs to use in updating w
    for (i in 1:length(K)){
        if (i != whichDSet){
            K[[i]] <- as.matrix(K[[i]])
            w[[i]] <- as.numeric(w[[i]])
            vx <- vx + K[[i]] %*% w[[i]]
        }
    }

    # Step 2 -- Normalize
    norm.vx <- as.numeric(sqrt(t(vx) %*% vx))
    if (is.na(norm.vx)) norm.vx <- 1
    if(norm.vx == 0) norm.vx <- 1
    vx <- vx / norm.vx

    # Step 3 - Apply penalty functions
    if (penalty == "SCAD"){
        u.new <- soft_threshold_scad(vx, tau)
    } else if (penalty == "LASSO"){
        u.new <- soft_threshold_lasso(vx, tau)
    } else {
        stop("Please provide the penalty function to be either LASSO or SCAD")
    }

    # Step 4 -- Normalize again
    norm.u.new <- as.numeric(sqrt(t(u.new) %*% u.new))
    if (is.na(norm.u.new)) norm.u.new <- 1
    if(norm.u.new == 0) norm.u.new <- 1
    u.new <- u.new / norm.u.new
    return(u.new)
}

#' Core ConvCCA Algorithm for Two Datasets
#'
#' @keywords internal
convCCA.algo <- function(X_1, X_2, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, nIter = NULL, penalty = "LASSO"){
    # Initializations
    X <- scale(X_1)
    Y <- scale(X_2)
    SigmaX <- cov(X)
    SigmaY <- cov(Y)
    SigmaXY <- cov(X, Y)

    if (ncol(X) > nrow(X)){
        SigmaXforK <- diag(ncol(X))
    } else {
        SigmaXforK <- SigmaX %^% (-0.5)
    }
    if (ncol(Y) > nrow(Y)){
        SigmaYforK <- diag(ncol(Y))
    } else {
        SigmaYforK <- SigmaY %^% (-0.5)
    }

    # Compute K
    K <- SigmaXforK %*% SigmaXY %*% SigmaYforK
    if (!any(is.na(K))){
        ee <- eigen(t(K) %*% K)
        d <- ee$values[1]
    }

    # Initialize
    if (is.null(initW_1)){
        newW_1 = rep(0.1, ncol(X))
    } else {
        newW_1 <- initW_1
    }
    if (is.null(initW_2)){
        newW_2 = rep(0.1, ncol(Y))
    } else {
        newW_2 <- initW_2
    }

    # Run the algorithm:
    if (is.null(nIter)){
        diff <- 1
        diffTemp <- 1
        while (diff > 1e-05 && !is.na(diffTemp)) {
            oldW_1 <- newW_1
            oldW_2 <- newW_2
            ## Update u:
            newW_1 <- updateW.convCCA(K = K, v = oldW_2, tau = tauW_1, penalty = penalty)
            ## Update v:
            newW_2 <- updateW.convCCA(K = t(K), v = newW_1, tau = tauW_2, penalty = penalty)
            diff0.v <- sqrt(sum((oldW_2 - newW_2) ^ 2))
            diff0.u <- sqrt(sum((oldW_1 - newW_1) ^ 2))
            diffTemp <- diff0.u + diff0.v
            if (!is.na(diffTemp)){
                diffU = diffTemp
            }
        }
    } else if (!is.numeric(nIter)){
        stop("Please provide a numerical value for the number of iterations(nIter)")
    } else {
        for (i in 1:nIter){
          oldW_1 <- newW_1
          oldW_2 <- newW_2
          newW_1 <- updateW.convCCA(K = K, v = oldW_2, tau = tauW_1, penalty = penalty)
          newW_2 <- updateW.convCCA(K = t(K), v = newW_1, tau = tauW_2, penalty = penalty)
        }
    }
    return(list("w_1" = newW_1, "w_2" = newW_2, "d" = d))
}

#' Core Multi-Dataset ConvCCA Algorithm
#'
#' @keywords internal
multi.convCCA.algo <- function(X, tau, nIter = 1000, penalty = "LASSO", initW = NULL){
    # Initializations
    length.list <- length(X)
    n.samples <- nrow(X[[1]])
    for (i in 1:length.list){
        X[[i]] <- scale(X[[i]])
        if (any(is.na(X[[i]]))){
            stop("Please remove/replace any missing data before implementing CCA.")
        }
        if (nrow(X[[i]]) != n.samples){
            stop("Please make sure that all data-sets in X have the same number of samples.")
        }
    }
    wNew <- vector("list", length = length.list)
    wOld <- vector("list", length = length.list)
    SigmaX <- vector("list", length = length.list)
    SigmaXforK <- vector("list", length = length.list)
    if (is.null(initW)){
        for (i in 1:length.list){
            wNew[[i]] <- rep(0.1, ncol(X[[i]]))
            wOld[[i]] <- wNew[[i]]
            SigmaX[[i]] <- cov(X[[i]])
            if (ncol(X[[i]]) > nrow(X[[i]])){
                SigmaXforK[[i]] <- diag(ncol(X[[i]]))
            } else {
                SigmaXforK[[i]] <- SigmaX[[i]] %^% (-0.5)
            }
        }
    } else {
        for (i in 1:length.list){
          wNew[[i]] <- initW[[i]]
          wOld[[i]] <- wNew[[i]]
          SigmaX[[i]] <- cov(X[[i]])
          if (ncol(X[[i]]) > nrow(X[[i]])){
              SigmaXforK[[i]] <- diag(ncol(X[[i]]))
          } else {
              SigmaXforK[[i]] <- SigmaX[[i]] %^% (-0.5)
          }
        }
    }
    diff <- 1
    n.t <- 0
    # Run the algorithm
    for (k in 1:nIter){
        for (i in 1:length.list){
            # Compute K matrices
            K <- vector("list", length = length.list-1)
            for (j in 1:length.list){
                if (j != i){
                    K[[j]] = SigmaXforK[[i]] %*% cov(X[[i]], X[[j]]) %*% SigmaXforK[[j]]
                }
            }
            wOld[[i]] <- wNew[[i]]
            ## Update u:
            wNew[[i]] <- updateWmulti.convCCA(K = K, w = wOld, whichDSet = i, tau = tau[[i]], penalty = penalty)
        }
    }
    return(list("w" = wNew))
}

#' Convex Sparse Canonical Correlation Analysis for Two Datasets
#'
#' @description Performs sparse CCA on two datasets using convex optimization
#' with LASSO or SCAD penalization. Based on Parkhomenko et al. (2009).
#'
#' @param X_1 First dataset matrix (n x p1), samples in rows, features in columns
#' @param X_2 Second dataset matrix (n x p2), same samples as X_1
#' @param tauW_1 Sparsity parameter for X_1 (higher = more sparse, typically 0.1-0.9)
#' @param tauW_2 Sparsity parameter for X_2
#' @param initW_1 Optional initial canonical vector for X_1 (default: rep(0.1, p1))
#' @param initW_2 Optional initial canonical vector for X_2 (default: rep(0.1, p2))
#' @param nIter Maximum iterations (NULL = auto-converge until diff < 1e-05, recommended)
#' @param penalty Penalty function: "LASSO" (default, L1) or "SCAD" (smoothly clipped absolute deviation)
#' @param R Number of canonical pairs to compute (default = 1)
#'
#' @return List containing:
#' \item{W_1}{Canonical vector for X_1 (length p1)}
#' \item{W_2}{Canonical vector for X_2 (length p2)}
#'
#' @export
#'
#' @examples
#' # Simulate two correlated datasets
#' set.seed(123)
#' n <- 50; p1 <- 100; p2 <- 80
#' u <- matrix(rnorm(n), ncol = 1)
#' v1 <- c(rep(1, 20), rep(0, 80))
#' v2 <- c(rep(1, 15), rep(0, 65))
#' X_1 <- u %*% t(v1) + matrix(rnorm(n * p1), n, p1)
#' X_2 <- u %*% t(v2) + matrix(rnorm(n * p2), n, p2)
#'
#' # Run ConvCCA
#' result <- convCCA(X_1, X_2, tauW_1 = 0.3, tauW_2 = 0.3, nIter = NULL)
#'
#' # Check sparsity
#' sum(abs(result$W_1) > 1e-6)  # Number of selected features in X_1
#'
convCCA <-  function(X_1, X_2, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, nIter = NULL, penalty = "LASSO", R = 1){
    if (R == 1){
        values = convCCA.algo(X_1 = X_1, X_2 = X_2, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, nIter = nIter, penalty = penalty)
        finalW_1 <- values$w_1
        finalW_2 <- values$w_2
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = convCCA.algo(X_1 = X_1, X_2 = X_2, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, nIter = nIter, penalty = penalty)
        tempW_1 <- values$w_1
        tempW_2 <- values$w_2
        for (i in 1:R-1){
            addX_1_1 <- t(tempW_1) %*% t(X_1) %*% X_1
            addX_1_2 <- t(tempW_2) %*% t(X_2) %*% X_1
            addX_2_1 <- t(tempW_2) %*% t(X_2) %*% X_2
            addX_2_2 <- t(tempW_1) %*% t(X_1) %*% X_2
            Xtilde_1 <- rbind(X_1, addX_1_1, addX_1_2)
            Xtilde_2 <- rbind(X_2, addX_2_1, addX_2_2)
            values = convCCA.algo(X_1 = Xtilde_1, X_2 = Xtilde_2, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, nIter = nIter, penalty = penalty)
            tempW_1 <- values$w_1
            tempW_2 <- values$w_2
        }
        finalW_1 <- tempW_1
        finalW_2 <- tempW_2
    }
    return(list("W_1" = finalW_1, "W_2" = finalW_2))
}

#' Multi-Dataset Convex Sparse Canonical Correlation Analysis
#'
#' @description Performs sparse CCA on multiple (>2) datasets simultaneously,
#' finding canonical vectors that maximize mutual correlation across all datasets
#'
#' @param X List of M datasets, where X[[i]] is an (n x p_i) matrix.
#'   All datasets must have the same number of samples (rows)
#' @param tau List of M sparsity parameters, one per dataset (e.g., list(0.3, 0.5, 0.2))
#' @param nIter Maximum number of iterations (default = 1000). Unlike convCCA,
#'   multi.convCCA requires a fixed nIter (no auto-convergence)
#' @param penalty Penalty function: "LASSO" (default) or "SCAD"
#' @param initW Optional list of M initial canonical vectors
#' @param R Number of canonical tuples to compute (default = 1)
#'
#' @return List containing:
#' \item{W}{List of M canonical vectors, one per dataset}
#'
#' @export
#'
#' @examples
#' # Simulate 3 correlated datasets
#' set.seed(123)
#' n <- 50
#' u <- matrix(rnorm(n), ncol = 1)
#' X1 <- u %*% t(c(rep(1, 20), rep(0, 80))) + matrix(rnorm(n * 100), n, 100)
#' X2 <- u %*% t(c(rep(1, 15), rep(0, 85))) + matrix(rnorm(n * 100), n, 100)
#' X3 <- u %*% t(c(rep(1, 10), rep(0, 40))) + matrix(rnorm(n * 50), n, 50)
#'
#' # Run multi.convCCA
#' result <- multi.convCCA(
#'   X = list(X1, X2, X3),
#'   tau = list(0.3, 0.3, 0.3),
#'   nIter = 100
#' )
#'
#' # Check sparsity for each dataset
#' sapply(result$W, function(w) sum(abs(w) > 1e-6))
#'
multi.convCCA <- function(X, tau, nIter = 1000, penalty = "LASSO", initW = NULL, R = 1){
    if (R == 1){
        values = multi.convCCA.algo(X = X, tau = tau, nIter = nIter, penalty = penalty, initW = initW)
        finalW <- values$w
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = multi.convCCA.algo(X = X, tau = tau, nIter = nIter, penalty = penalty, initW = initW)
        tempW <- values$w
        for (k in 1:(R-1)){
            Xtilde <- vector("list", length = length(X))
            addX <- vector("list", length = length(X))
            for (i in 1:length(X)){
                for (j in 1:length(X)){
                    addX[[j]] <- t(tempW[[j]]) %*% t(X[[j]]) %*% X[[i]]
                }
                Xtilde[[i]] <- X[[i]]
                for (j in 1:length(X)){
                    Xtilde[[i]] <- rbind(Xtilde[[i]], addX[[j]])
                }
            }
            values = multi.convCCA.algo(X = Xtilde, tau = tau, nIter = nIter, penalty = penalty, initW = initW)
            for (j in 1:length(tempW)){
                tempW[[j]] <- cbind(tempW[[j]], values$w[[j]])
            }
        }
        finalW <- tempW
    }
    return(list("W" = finalW))
}
