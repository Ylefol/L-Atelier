###############################################################################
########### RelPMDCCA - Relaxed Penalized Matrix Decomposition CCA ###########
###############################################################################

# Import dependencies from other GAIA modules using source
source('~/A_Projects/ZERO_DAWN/GAIA/Minerva/penalty_functions.R')
source('~/A_Projects/ZERO_DAWN/GAIA/Demeter/matrix_utils.R')

#' Update Consensus Variable (z-update)
#'
#' @description Projects onto the unit ball for consensus constraint in ADMM.
#' Ensures ||z|| <= 1 by normalization if needed.
#'
#' @param x Current value before projection
#' @param old Previous value (fallback if NA encountered)
#'
#' @return Projected consensus variable z
#'
#' @export
HEPHAESTUS_updateZ <- function(x, old){
  check <- sqrt(sum(x^2))
  condition <- (check <= 1)
  if (!is.na(condition)){
      if (condition){
        update_z <- x
        } else {
          update_z <- x/check
      }
  } else {
      update_z = old
  }
  return(update_z)
}

#' Core RelPMDCCA Algorithm for Two Datasets
#'
#' @keywords internal
#' @export
HEPHAESTUS_relPMDCCA_algo <- function(X_1, X_2, lambda, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, a = 3.7, sd.mu = sqrt(2), nIter = 100, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL){
    # We use X notation for X_1 and Y for X_2
    X <- scale(X_1)
    Y <- scale(X_2)
    # mu must satisfy: 0 < mu <= lambda/norm(X)^2
    ## Each data-set has its own mu (and lambda) values - we fix lambda and have separate mu
    mu.x <- rnorm(1, mean = 0, sd = sd.mu)
    if ((mu.x < 0) | (mu.x > lambda/(Matrix::norm(X, "f")^2))){
        mu.x <- runif(1, min = 0, max = lambda/(Matrix::norm(X, "f")^2))
    }
    mu.y <- rnorm(1, mean = 0, sd = sd.mu)
    if ((mu.y < 0) | (mu.y > lambda/(Matrix::norm(Y, "f")^2))){
        mu.y <- runif(1, min = 0, max = lambda/(Matrix::norm(Y, "f")^2))
    }
    # Initializations
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
    newZ <- rep(0.4, nrow(X))
    newXi <- rep(1, nrow(X))
    # Initialize convergence constraint
    diff <- 1
    diffTemp <- 1
    diffW_1 <- 1
    diffW_2 <- 1
    n.t.w_1 <- 0
    n.t <- 0
    k <- 0
    k.w_1 <- 0
    k.w_2 <- 0
    # Run the algorithm - Iterative updated
    while ((diff > 1e-05) && (!is.na(diffTemp)) && (n.t < nIter)){
        k <- k + 1
        n.t.w_1 <- 0
        diffW_1 <- 1
        diffTemp <- 1
        oldW_1 <- newW_1
        oldW_2 <- newW_2
        oldZ <- newZ
        oldXi <- newXi
        # Update First dataset
        while ((diffW_1 > 1e-05) && (!is.na(diffTemp)) && (n.t.w_1 < nIter)){
          k.w_1 <- k.w_1 + 1
          currentW_1 <- newW_1
          currentW_2 <- newW_2
          currentZ <- newZ
          currentXi <- newXi
          # Compute omega:
          omega <- currentW_1 - ( (mu.x/lambda)*(t(X) %*% (X %*% currentW_1 - currentZ + currentXi)) )
          # Compute c:
          c = t(X) %*% Y %*% currentW_2
          for (j in 1:length(newW_1)){
              # Use penalty functions from Minerva
              if (penalty == "LASSO") {
                  newW_1[j] <- MINERVA_proximal_lasso(omega[j], mu.x, c[j], tauW_1, currentW_1[j])
              } else if (penalty == "SCAD") {
                  newW_1[j] <- MINERVA_proximal_scad(omega[j], mu.x, c[j], tauW_1, a, currentW_1[j])
              } else if (penalty == "ELASTIC-NET") {
                  newW_1[j] <- MINERVA_proximal_elastic_net(omega[j], mu.x, c[j], tauW_1, tau_EN, currentW_1[j])
              } else {
                  stop("Please provide the penalty function to be either LASSO, SCAD, or ELASTIC-NET")
              }
          }
          diffTemp <- sqrt(sum((currentW_1 - newW_1)^2))
          if (!is.na(diffTemp)){
            diffW_1 <- diffTemp
          }
          check_Z <- X %*% newW_1 + currentXi
          newZ <- HEPHAESTUS_updateZ(check_Z, currentZ)
          newXi <- currentXi + X %*% newW_1 - newZ
          n.t.w_1 <- n.t.w_1 + 1
        }
        diffTemp <- 1
        diffW_2 <- 1
        n.t.w_2 <- 0
        # Update Second dataset
        while ((diffW_2 > 1e-05) && (!is.na(diffTemp)) && (n.t.w_2 < nIter)){
          k.w_2 <- k.w_2 + 1
          currentW_1 <- newW_1
          currentW_2 <- newW_2
          currentZ <- newZ
          currentXi <- newXi
          # Compute omega:
          omega <- oldW_2 - ( (mu.y/lambda)*(t(Y) %*% (Y %*% oldW_2 - currentZ + currentXi)) )
          # Compute c:
          c = t(currentW_1) %*% t(X) %*% Y
          for (j in 1:length(newW_2)){
              # Use penalty functions from Minerva
              if (penalty == "LASSO") {
                  newW_2[j] <- MINERVA_proximal_lasso(omega[j], mu.y, c[j], tauW_2, currentW_2[j])
              } else if (penalty == "SCAD") {
                  newW_2[j] <- MINERVA_proximal_scad(omega[j], mu.y, c[j], tauW_2, a, currentW_2[j])
              } else if (penalty == "ELASTIC-NET") {
                  newW_2[j] <- MINERVA_proximal_elastic_net(omega[j], mu.y, c[j], tauW_2, tau_EN, currentW_2[j])
              } else {
                  stop("Please provide the penalty function to be either LASSO, SCAD, or ELASTIC-NET")
              }
          }
          diffTemp <- sqrt(sum((currentW_2 - newW_2)^2))
          if (!is.na(diffTemp)){
            diffW_2 <- diffTemp
          }
          check_Z <- Y %*% newW_2 + currentXi
          newZ <- HEPHAESTUS_updateZ(check_Z, currentZ)
          newXi <- currentXi + Y %*% newW_2 - newZ
          n.t.w_2 <- n.t.w_2 + 1
        }
        diff0.w_2 <- sqrt(sum((oldW_2 - newW_2)^2))
        diff0.w_1 <- sqrt(sum((oldW_1 - newW_1)^2))
        diff0.z <- sqrt(sum((oldZ - newZ)^2))
        diff0.xi <- sqrt(sum((oldXi - newXi)^2))
        diffTemp <- diff0.w_1 + diff0.w_2 + diff0.z + diff0.xi
        if (!is.na(diffTemp)){
          diff <- diffTemp
        }
        n.t <- n.t + 1
    }
    return(list("w_1" = newW_1, "w_2" = newW_2, "z" = newZ, "xi" = newXi))
}

#' Core Multi-Dataset RelPMDCCA Algorithm
#'
#' @keywords internal
#' @export
HEPHAESTUS_multi_relPMDCCA_algo <- function(X, lambda, tau, a = 3.7, sd.mu = sqrt(2), nIter = 1000, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL){
    # Initializations
    length.list <- length(X)
    n.samples = nrow(X[[1]])
    for (i in 1:length.list){
        X[[i]] <- scale(X[[i]])
        if (any(is.na(X[[i]]))){
            stop("Please remove/replace any missing data before implementing sCCA.")
        }
        if (nrow(X[[i]]) != n.samples){
            stop("Please make sure that all datasets in list X have the same number of samples")
        }
    }
    # mu must satisfy: 0 < mu <= lambda/norm(X)^2
    ## Each data-set has its own mu (and lambda) values - we fix lambda and have separate mu
    mu.w <- vector("list", length = length.list)
    newW <- vector("list", length = length.list)
    for (i in 1:length.list){
      mu.w[[i]] <- rnorm(1, mean = 0, sd = sd.mu)
      if ((mu.w[[i]] < 0) | (mu.w[[i]] > lambda/(Matrix::norm(X[[i]], "f")^2))){
        mu.w[[i]] <- runif(1, min = 0, max = lambda/(Matrix::norm(X[[i]], "f")^2))
      }
      newW[[i]] <- rep(0.1, ncol(X[[i]]))
    }
    newZ <- rep(0.4, nrow(X[[1]]))
    newXi <- rep(1, nrow(X[[1]]))

    # Run the algorithm
    for (i in 1:nIter){
        for (j in 1:length.list){
            oldW <- newW
            oldZ <- newZ
            oldXi <- newXi
            # Compute omega
            omega <- oldW[[j]] - ( (mu.w[[j]]/lambda)*(t(X[[j]]) %*% (X[[j]] %*% oldW[[j]] - oldZ + oldXi)) )
            # Compute c:
            Yv <- 0
            for (k in 1:length.list){
                if (k != j){
                    Yv <- Yv + X[[k]] %*% oldW[[k]]
                }
            }
            cc <- t(X[[j]]) %*% Yv
            for (k in 1:length(newW[[j]])){
                # Use penalty functions from Minerva
                if (penalty == "LASSO") {
                    newW[[j]][k] <- MINERVA_proximal_lasso(omega[k], mu.w[[j]], cc[k], tau[[j]], oldW[[j]][k])
                } else if (penalty == "SCAD") {
                    newW[[j]][k] <- MINERVA_proximal_scad(omega[k], mu.w[[j]], cc[k], tau[[j]], a, oldW[[j]][k])
                } else if (penalty == "ELASTIC-NET") {
                    newW[[j]][k] <- MINERVA_proximal_elastic_net(omega[k], mu.w[[j]], cc[k], tau[[j]], tau_EN, oldW[[j]][k])
                } else {
                    stop("Please provide the penalty function to be either LASSO, SCAD, or ELASTIC-NET")
                }
            }
            check_Z <- X[[j]] %*% newW[[j]] + oldXi
            newZ <- HEPHAESTUS_updateZ(check_Z, oldZ)
            newXi <- oldXi + X[[j]] %*% newW[[j]] - newZ
        }
    }
    return(list("w" = newW, "z" = newZ, "xi" = newXi))
}

#' Relaxed Penalized Matrix Decomposition CCA for Two Datasets
#'
#' @description Performs sparse CCA using ADMM framework with relaxed independence
#' assumption. Allows for dependent features within datasets. Based on Suo et al. (2017).
#'
#' @param X_1 First dataset matrix (n x p1), samples in rows, features in columns
#' @param X_2 Second dataset matrix (n x p2), same samples as X_1
#' @param lambda Relaxation parameter controlling step size (typically fixed at 10)
#' @param tauW_1 Sparsity parameter for X_1 (higher = more sparse, typically 0.5-0.9)
#' @param tauW_2 Sparsity parameter for X_2
#' @param initW_1 Optional initial canonical vector for X_1 (default: rep(0.1, p1))
#' @param initW_2 Optional initial canonical vector for X_2 (default: rep(0.1, p2))
#' @param a SCAD shape parameter (default = 3.7, from Fan & Li 2001)
#' @param sd.mu Standard deviation for step size initialization (default = sqrt(2))
#' @param nIter Maximum iterations (default = 1000). Algorithm stops early if converged.
#' @param penalty Penalty function: "LASSO" (default), "SCAD", or "ELASTIC-NET"
#' @param element_wise Use element-wise proximal updates (default = TRUE, recommended)
#' @param tau_EN Secondary tuning parameter for Elastic-Net penalty
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
#' # Run RelPMDCCA
#' result <- HEPHAESTUS_relPMDCCA(X_1, X_2, lambda = 10, tauW_1 = 0.8, tauW_2 = 0.8, nIter = 100)
#'
#' # Check canonical correlation
#' cor(X_1 %*% result$W_1, X_2 %*% result$W_2)
#'
HEPHAESTUS_relPMDCCA <- function(X_1, X_2, lambda, tauW_1, tauW_2, initW_1 = NULL, initW_2 = NULL, a = 3.7, sd.mu = sqrt(2), nIter = 1000, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL, R = 1){
    if (R == 1){
        values = HEPHAESTUS_relPMDCCA_algo(X_1 = X_1, X_2 = X_2, lambda = lambda, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        finalW_1 <- values$w_1
        finalW_2 <- values$w_2
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = HEPHAESTUS_relPMDCCA_algo(X_1 = X_1, X_2 = X_2, lambda = lambda, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        tempW_1 <- values$w_1
        tempW_2 <- values$w_2
        for (i in 1:(R-1)){
            addX_1_1 <- t(tempW_1) %*% t(X_1) %*% X_1
            addX_1_2 <- t(tempW_2) %*% t(X_2) %*% X_1
            addX_2_1 <- t(tempW_2) %*% t(X_1) %*% X_2
            addX_2_2 <- t(tempW_1) %*% t(X_1) %*% X_2
            Xtile_1 <- rbind(X_1, addX_1_1, addX_1_2)
            Xtilde_2 <- rbind(X_2, addX_2_1, addX_2_2)
            values = HEPHAESTUS_relPMDCCA_algo(X_1 = Xtilde_1, X_2 = Xtilde_2, lambda = lambda, tauW_1 = tauW_1, tauW_2 = tauW_2, initW_1 = initW_1, initW_2 = initW_2, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
            tempW_1 <- cbind(tempW_1, values$w_1)
            tempW_2 <- cbind(tempW_2, values$w_2)
        }
        finalW_1 <- tempW_1
        finalW_2 <- tempW_2
    }
    return(list("W_1" = finalW_1, "W_2" = finalW_2))
}

#' Multi-Dataset Relaxed Penalized Matrix Decomposition CCA
#'
#' @description Performs sparse CCA on multiple (>2) datasets simultaneously using
#' ADMM framework. Finds canonical vectors that maximize mutual correlation across
#' all datasets while allowing for dependent features.
#'
#' @param X List of M datasets, where X[[i]] is an (n x p_i) matrix.
#'   All datasets must have the same number of samples (rows)
#' @param lambda Relaxation parameter controlling step size (typically fixed at 10)
#' @param tau List of M sparsity parameters, one per dataset (e.g., list(0.8, 0.8, 0.7))
#' @param a SCAD shape parameter (default = 3.7)
#' @param sd.mu Standard deviation for step size initialization (default = sqrt(2))
#' @param nIter Maximum number of iterations (default = 1000). Stops early if converged.
#' @param penalty Penalty function: "LASSO" (default), "SCAD", or "ELASTIC-NET"
#' @param element_wise Use element-wise proximal updates (default = TRUE)
#' @param tau_EN Secondary tuning parameter for Elastic-Net penalty
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
#' # Run multi.relPMDCCA
#' result <- HEPHAESTUS_multi_relPMDCCA(
#'   X = list(X1, X2, X3),
#'   lambda = 10,
#'   tau = list(0.8, 0.8, 0.8),
#'   nIter = 100
#' )
#'
#' # Check mutual correlations
#' scores <- lapply(1:3, function(i) result$X[[i]] %*% result$W[[i]])
#' cor(scores[[1]], scores[[2]])
#'
HEPHAESTUS_multi_relPMDCCA <- function(X, lambda, tau, a = 3.7, sd.mu = sqrt(2), nIter = 1000, penalty = "LASSO", element_wise = TRUE, tau_EN = NULL, R = 1){
    if (R == 1){
        values = HEPHAESTUS_multi_relPMDCCA_algo(X = X, lambda = lambda, tau = tau, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
        finalW <- values$w
    } else if (R < 1){
        stop("Please provide a positive integer for the number of canonical pairs (R)")
    } else {
        values = HEPHAESTUS_multi_relPMDCCA_algo(X = X, lambda = lambda, tau = tau, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
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
            values = HEPHAESTUS_multi_relPMDCCA_algo(X = Xtilde, lambda = lambda, tau = tau, a = a, sd.mu = sd.mu, nIter = nIter, penalty = penalty, element_wise = element_wise, tau_EN = tau_EN)
            for (j in 1:length(tempW)){
                tempW[[j]] <- cbind(tempW[[j]], values$w[[j]])
            }
        }
        finalW <- tempW
    }
    return(list("W" = finalW))
}
