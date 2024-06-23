#' Computes representative vectors for the Jacobians evaluated at data points and returns the matrix of eigenvectors obtained from the SPD matrix build from such vectors.
#'
#' @param grads Gradient matrix, where the columns are the concatenation of the objective-wise gradients
#' @param n Sample size
#' @param m Number of objectives or function outputs
#' @param d Number of input dimensions
#'
#' @param method The method to compute a reprentative vector, either
#' \itemize{
#' \item 'AG' computes the average of gradient vectors
#' \item 'LP'  for each input variable, computes the rank-1 linear projection of the objective-wise derivatives wrt that variable
#' (based on the leading eigenvector of the SPD matrix of these same derivatives).
#'
#' \item 'MCH' computes the minimum-norm element of the convex hull computed from the Jacobian
#' }
#'
#' @param is_grad_norm If set to \code{TRUE}, performs the L2 normalization of the gradients. This is used only for AG and MCH.
#' @export
#' @return Matrix of eigenvectors (column vectors)
#'
#' @examples
#' fn <- problem_car_side_impact
#' d <- 7
#' n <- 1000
#' X <- matrix(runif(d * n), n)
#' Y <- t(apply(X, 1, fn))
#' grads <- t(apply(X, 1, function(x) c(t(numDeriv::jacobian(fn, x = x)))))
#' methods_grad(grads, n, ncol(Y), d, method = "MCH")
#'
methods_grad <- function(grads, n, m, d, method = c("AG", "LP", "MCH"), is_grad_norm = FALSE) {
  out <- matrix(0, nrow = n, ncol = d)

  if (method == "LP") {
    out <- project_gradients(grads, n, m, d)
  } else {
    for (i in 1:n) {
      J <- matrix(grads[i, ], ncol = d, byrow = T)

      if (is_grad_norm) J <- J / (sqrt(rowSums(J^2)))

      if (method == "MCH") {
        out[i, ] <- min_convex_hull(J)
      } else {
        out[i, ] <- colSums(J) / m
      }
    }
  }


  C <- (t(out) %*% out) / n

  return(eigen(C)$vectors)
}



#' For each input parameter, rank-1 linear projection of the objective-wise derivatives wrt that variable
#' based on the leading eigenvector obtained from the eigendecomposition of the SPD matrix of these same derivatives
#'
#' @param grads Matrix of gradients
#' @param n Sample size
#' @param m Number of objectives
#' @param d Number of input dimensions
#' @noRd
#' @return Matrix  of projected gradients
project_gradients <- function(grads, n, m, d) {
  projected_grads <- matrix(0, nrow = n, ncol = d)

  P <- list()

  for (i in 1:d) {
    P[[i]] <- matrix(0, ncol = m, nrow = n)

    for (j in 1:m) {
      P[[i]][, j] <- grads[, ((j - 1) * d + i)]
    }

    P_dim <- t(P[[i]]) %*% P[[i]] / n

    e <- eigen(P_dim)

    projected_grads[, i] <- P[[i]] %*% e$vectors[, 1]
  }

  return(projected_grads)
}
