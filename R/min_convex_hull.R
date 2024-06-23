#' A recursive algorithm to compute the minimum L2-norm point in a polytope.
#'
#' @param P Matrix for the convex hull based on which to compute the minimum L2-norm element
#'
#' @return Vector representing the minimum L2-norm point of the convex hull of \code{P}
#' @export
#' @references Sekitani, K. and Y. Yamamoto (1993). A recursive algorithm for finding the minimum
#' norm point in a polytope and a pair of closest points in two polytopes. Mathematical Programming 61, 233–249
#'
#' @examples
#' X <- matrix(runif(3 * 10), 10)
#' min_convex_hull(X)
min_convex_hull <- function(P) {
  if (is.null(dim(P))) P <- matrix(P, nrow = 1)

  if (nrow(P) == 1) {
    return(P)
  }

  x <- P[which.min(rowSums(P^2)), ]

  tol <- 1e-10

  while (TRUE) {
    # Step 1
    if (is.null(dim(x))) x <- matrix(x, ncol = length(x))

    p_alpha <- P[which.min(rowSums(sweep(P, MARGIN = 2, x, `*`))), ]

    if (sum(x * x) <= sum(x * p_alpha)) {
      return(x)
    }

    # Step 2
    ind <- apply(P, 1, function(p) sum(x * (p - p_alpha)) < tol)

    Pk <- P[ind, ]

    if (is.null(dim(Pk))) Pk <- matrix(Pk, ncol = length(Pk))

    if (nrow(Pk) == nrow(P)) {
      return(x)
    }

    y <- min_convex_hull(Pk)

    # Step 3

    P_Pk <- P[!apply(P, 1, function(x) any(apply(Pk, 1, function(y) all(x == y)))), ]

    if (is.null(dim(P_Pk))) P_Pk <- matrix(P_Pk, ncol = length(P_Pk))

    p_beta <- P_Pk[which.min(rowSums(sweep(P_Pk, MARGIN = 2, y, `*`))), ]

    if (sum(y * (y - p_beta)) < tol) {
      return(y)
    }

    # Step 4
    P_aux <- P_Pk[apply(P_Pk, 1, function(p) sum((y - x) * (y - p)) > tol & sum(x * (y - p)) != 0), ]

    if (length(P_aux) == 0) {
      return(x)
    }

    if (is.null(dim(P_aux))) P_aux <- matrix(P_aux, ncol = length(P_aux))

    p_lambda <- P_aux[which.min(apply(P_aux, 1, function(p) sum(y * (y - p)) / sum(x * (y - p)))), ]

    if (is.null(dim(p_lambda))) p_lambda <- matrix(p_lambda, ncol = length(p_lambda))

    lam <- sum(x * (p_lambda - y)) / sum((y - x) * (y - p_lambda))

    x <- x + lam * (y - x)
  }
}
