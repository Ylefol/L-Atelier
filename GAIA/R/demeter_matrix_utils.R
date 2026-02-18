###############################################################################
########### Matrix Utilities ###########
###############################################################################

#' Matrix Power Operator
#'
#' @description Raises a matrix to a given power using eigenvalue decomposition
#'
#' @param x A square matrix
#' @param n Power to raise the matrix to (can be negative for inverse powers)
#'
#' @return Matrix x raised to the power n
#'
#' @examples
#' \dontrun{
#' # Compute inverse square root of a covariance matrix
#' Sigma <- cov(matrix(rnorm(100), 10, 10))
#' Sigma_inv_sqrt <- Sigma %^% (-0.5)
#'
#' }
#' @export
"%^%" <- function(x, n){
    with(eigen(x), vectors %*% (values^n * t(vectors)))
}

#' Matrix Power Function (Function Form)
#'
#' @description Raises a matrix to a given power using eigenvalue decomposition.
#' Function form of the %^% operator for explicit module referencing.
#'
#' @param x A square matrix
#' @param n Power to raise the matrix to (can be negative for inverse powers)
#'
#' @return Matrix x raised to the power n
#'
#' @examples
#' \dontrun{
#' # Compute inverse square root of a covariance matrix
#' Sigma <- cov(matrix(rnorm(100), 10, 10))
#' Sigma_inv_sqrt <- DEMETER_matrix_power(Sigma, -0.5)
#'
#' }
#' @export
DEMETER_matrix_power <- function(x, n){
    with(eigen(x), vectors %*% (values^n * t(vectors)))
}

#' Compute Pseudoinverse of a Matrix
#'
#' @description Computes the Moore-Penrose pseudoinverse using SVD.
#' Handles rank-deficient matrices by truncating near-zero singular values.
#'
#' @param X Input matrix
#' @param tolerance Singular values below this threshold are treated as zero (default: 1e-3)
#'
#' @return Pseudoinverse of X
#'
#' @export
DEMETER_gettingInverse <- function(X, tolerance = 1e-3){
    A <- X
    svdA <- svd(A)
    d <- round(svdA$d, 3)
    or <- order(d, decreasing = TRUE)
    U <- svdA$u[or, ]
    V <- svdA$v[, or]
    d <- d[or]
    l <- length(which(d == 0))
    val <- nrow(A) - l
    S <- diag(svdA$d)
    Splus <- diag(c(1/svdA$d[1:val], rep(0, l)))
    Aplus <- V %*% Splus %*% t(U)
    return(Aplus)
}
