#' Extracts the first \code{j} dominating eigenvectors \code{Wj}, computes the rank-j reconstruction of the data and applies the given test function on the reconstructed data.
#'
#' @param evectors Matrix of eigenvectors
#' @param fn Test function
#' @param X Matrix of input data
#' @param num_dim Number of dominating eigenvectors
#' @param which_distr Data distribution
#' @param norm_params Normalization parameters
#' @param params List of additional parameter to be passed to to the test function
#' @noRd
#' @return List of dominating eigenspace \code{W1} and the function estimations \code{Yest} on the reconstructed data
project_and_estimate <- function(evectors, fn, X, num_dim, which_distr, norm_params, params) {

    W1 <- evectors[, 1:num_dim] #extract important eigenvectors

    X_projected <- X %*% W1 %*% t(W1)

    Yest <- wrapper(X_projected, fn, which_distr, norm_params, params)

    return(list(W1=W1, Yest=Yest))

      }

