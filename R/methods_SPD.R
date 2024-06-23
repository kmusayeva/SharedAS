#' Returns the matrix of eigenvectors computed from the stepwise estimation of eigenvectors, joint diagonalization based on the FG algorithm, or from the eigendecomposition of the sum of SPD matrices.
#'
#' @param grads Gradient matrix, where the columns are the concatenation of the objective-wise gradients
#' @param n Sample size
#' @param m Number of objectives
#' @param d Number of input dimensions
#' @param method Name of the SPD method to choose from
#' \itemize{
#' \item 'SSPD' to compute the eigenvectors from the sum of SPD matrices,
#' \item 'SEE' to compute the eigenvectors in a stepwise fashion, or
#' \item 'FG' to do simultaneous diagonalization using the Flury and Gautschi algorithm.
#' }
#'
#' @export
#' @return Matrix of eigenvectors (column vectors)
#'
#' @references
#' \itemize{
#' \item{Flury, B. N. and W. Gautschi (1986). An algorithm for simultaneous orthogonal transformation of several positive definite symmetric matrices to nearly diagonal form. SIAM Journal on Scientific and Statistical Computing 7 (1), 169–184}
#' \item{Trendafilov, N. T. (2010). Stepwise estimation of common principal components. Computational Statistics & Data Analysis 54 (12), 3446–3457.}
#' }
#' @examples
#' fn <- problem_car_side_impact
#' d <- 7
#' n <- 1000
#' X <- matrix(runif(d * n), n)
#' Y <- t(apply(X, 1, fn))
#' grads <- t(apply(X, 1, function(x) c(t(numDeriv::jacobian(fn, x = x)))))
#' methods_SPD(grads, n, ncol(Y), d, method = "FG")
#'
methods_SPD <- function(grads, n, m, d, method = c("SSPD", "SEE", "FG")) {
  if (!method %in% c("SSPD", "SEE", "FG")) {
    stop(cat("Possible choices: 'SSPD', 'SEE' or 'FG'\n"))
  }


  if (method == "SSPD") {
    return(sumSPD(grads, n, m, d))
  }

  S <- array(NA, c(d, d, m))

  for (i in 1:m) {
    g <- grads[, ((i - 1) * d + 1):(i * d)]

    S[, , i] <- (t(g) %*% g / n) + diag(1e-12, d)
  }

  if (method == "FG") {
    FG_res <- FG(covmats = S, nvec = rep(n, m), method = "LS") ### it can also be ML, but LS gave better results
  } else {
    FG_res <- stepwisecpc(covmats = S, nvec = rep(n, m))
  }

  return(FG_res$B)
}





#' Computes eigenvectors from the sum of SPD matrices
#' @param grads Matrix of gradients
#' @param n Sample size
#' @param m Number of objectives
#' @param d Number of input dimensions
#' @noRd
#' @return Matrix of eigenvectors
sumSPD <- function(grads, n, m, d) {
  S <- matrix(0, nrow = d, ncol = d)

  for (i in 1:m) {
    g <- grads[, ((i - 1) * d + 1):(i * d)]

    S <- S + t(g) %*% g / n
  }

  return(eigen(S)$vectors)
}
