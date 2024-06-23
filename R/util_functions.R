#' Computes the jacobian of a vector-valued function
#'
#' @param fn Test function
#' @param X Input data matrix
#' @param distr Data distribution
#' @param norm_params If \code{TRUE}, mean-variance normalizes the function outputs
#' @param params Additional parameters to the test function
#' @noRd
#' @return Matrix where rows contain objective-wise gradients
compute_grad <- function(fn, X, which_distr, norm_params = NULL, params = NULL) {
  grads <- t(apply(X, 1, function(x) c(t(jacobian(wrapper, x = x, fn = fn, which_distr = which_distr, norm_params = norm_params, params = params)))))

  return(grads)
}


#' Maps the input data to the unit hypercube, and if specified, normalizes the outputs
#' @param fn Test function
#' @param x Input matrix or vector
#' @param which_distr Data distribution
#' @param norm_params If \code{TRUE}, mean-variance normalizes the function outputs (logical)
#' @param params additional parameters to the function \code{fn}
#' @noRd
#' @return vector of function outputs.
wrapper <- function(x, fn, which_distr, norm_params = NULL, params = NULL) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1) # if not a matrix, convert to one

  if (which_distr == "norm") {
    x <- sigmoid(x)
  } # map to [0,1]
  else {
    x <- (x + 1) / 2
    x <- pmin(pmax(x, 0), 1)
  }

  if (is.null(params)) {
    res <- t(apply(x, 1, function(x) do.call(fn, list(x))))
  } else {
    res <- t(apply(x, 1, function(x) do.call(fn, list(x, params))))
  }


  if (!is.null(norm_params)) {
    res <- mean_var_norm(res, norm_params$mean, norm_params$std)
  }

  return(res)
}




#' Computes eigenvectors for each method
#'
#' @param methods Vector of methods
#' @param grads Matrix of gradients
#' @param n Numberof samples
#' @param m Number of function outputs
#' @param d Input dimensionality
#'
#' @return Matrix of eigenvectors
#' @noRd
compute_evectors <- function(methods, grads, n, m, d) {
  evectors <- list()


  for (method in methods) {
    # start_time <- Sys.time()

    evectors[[method]] <- switch(method,
      "SSPD" = methods_SPD(grads, n, m, d, method = "SSPD"),
      "SEE" = methods_SPD(grads, n, m, d, method = "SEE"),
      "FG" = methods_SPD(grads, n, m, d, method = "FG"),
      "AG" = methods_grad(grads, n, m, d, method = "AG"),
      "MCH" = methods_grad(grads, n, m, d, method = "MCH"),
      "LP" = methods_grad(grads, n, m, d, method = "LP")
    )

    # end_time <- Sys.time()

    # duration <- as.numeric(end_time - start_time, units = "secs")

    # cat("method", method, ", runtime:", round(duration, 4))
  }



  return(evectors)
}



#' Generates a data according to a given distribution, if \code{X} is not provided, applies the given function to the data, computes the normalization parameters
#'
#' @param fn Test function
#' @param X User provided data matrix, optional
#' @param input_dim Input dimensionality
#' @param num_obj Number of function outputs
#' @param distr Data distribution
#' @param params  List pf additional parameters to the test function
#'
#' @return Matrix of generated data, matrix of function outputs on these data, number of objectives, matrix of normalized outputs and list normalization paramaters
#'
#' @noRd
handle_data <- function(fn, X = NULL, input_dim, num_obj, distr, params = NULL) {
  n <- distr$n

  d <- input_dim

  which_distr <- distr$name

  if (is.null(X)) X <- if (distr$name == "unif") matrix(runif(d * n, min = -1, max = 1), n) else rmvnorm(n = n, mean = rep(0, d), sigma = distr$sigma)

  Y <- wrapper(X, fn, which_distr, norm_params = NULL, params)

  if (!is.null(params$t)) { # penicillin problem

    ind <- which(Y[, 1] == params$t)

    Y <- Y[ind, ]

    X <- X[ind, ]

    num_obj <- ncol(Y) - 1
  }

  # compute the mean and the std of the function outputs for the mean-variance normalization
  Y_mean <- colMeans(Y)
  Y_std <- apply(Y, 2, sd)
  norm_params <- list(mean = Y_mean, std = Y_std)

  Y_norm <- mean_var_norm(Y, Y_mean, Y_std)

  return(list(X = X, Y = Y, num_obj = num_obj, Y_norm = Y_norm, norm_params = norm_params))
}



#' Computes the gradients of the test function.
#'
#' @param fn Test function
#' @param X Matrix of input data
#' @param which_distr Data distribution
#' @param is_norm If \code{TRUE} normalizes when computing the gradients
#' @param norm_params Normalization parameters
#' @param params List pf additional parameters to the test function
#'
#' @noRd
#'
handle_grads <- function(fn, X, which_distr, is_norm = FALSE, norm_params = NULL, params = NULL) {
  grads <- if (is_norm) compute_grad(fn, X, which_distr, norm_params, params) else compute_grad(fn, X, which_distr, NULL, params)

  if (!is.null(params$t)) {
    grads <- grads[, -c(1:ncol(X))]
  } # for the penicillin problem

  return(grads)
}




#' Mean-variance normalizes the given matrix
#'
#' @param A Matrix whose columns are to be normalized
#' @param mean Column-wise mean of the function outputs on the original space
#' @param std Column-wise standard deviation of the function outputs on the original space
#' @noRd
#' @return Mean-variance normalized matrix
mean_var_norm <- function(A, mean, std) {
  if (is.null(dim(A))) A <- matrix(A, nrow = 1)

  A_norm <- t(t(sweep(A, 2, mean)) / std)

  ind <- which(std == 0)

  A_norm[, ind] <- A[, ind]

  return(A_norm)
}




#' Min-max normalizes the given matrix
#'
#' @param A Matrix whose columns are to be normalized
#' @noRd
#' @return Min-max normalized matrix
min_max_norm <- function(A) {
  if ((max(A) - min(A)) != 0) {
    return((A - min(A)) / (max(A) - min(A)))
  }

  return(A)
}



#' Computes the root-mean-squared error (rmse)
#'
#' @param A matrix
#' @param B matrix
#' @noRd
#' @return Vector of column-wise rmses
rmse <- function(A, B) {
  # not needed for already normalized matrices
  # A <- mean_var_norm(A)
  # B <- mean_var_norm(B)

  sqrt(colMeans((A - B)^2))
}


#' Does any additional checking before computing the rmse
#' @param Y function outputs on the original input space
#' @param Yest function outputs on the reconstructed space
#' @noRd
#' @return Vector of objective-wise rmses and their sums.
rmse_wrapper <- function(Y, Yest, params) {
  if (!(is.null(params$t))) { # for penicillin problem, this should be improved because other functions may contain such a parameter
    ind <- which(Yest[, 1] == params$t)
    Yest <- Yest[ind, -1]
    rmses <- rmse(as.matrix(Y[ind, -1]), as.matrix(Yest))
  } else {
    rmses <- rmse(as.matrix(Y), as.matrix(Yest))
  }

  return(c(rmses, sum(rmses)))
}



#' Creates a matrix for bookeeping rmses.
#'
#' @param d Number of input dimensions
#' @param m Number of objectives
#' @param nexps Number of experiments
#' @param methods List of methods
#' @noRd
#' @return List of matrices
create_Result_matrix <- function(start_dim, end_dim, m, nexps, methods) {
  rmses <- matrix(rep(0, nexps * (m + 1)), nrow = nexps)

  colnames(rmses) <- c(sprintf("f%d", 1:m), "sum")

  # L <- vector("list", d)

  L <- list()

  for (i in start_dim:end_dim) L[[i]] <- rmses

  Result <- list()

  for (method in methods) Result[[method]] <- L

  return(Result)
}
