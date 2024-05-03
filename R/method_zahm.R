#' Implements the method of Zahm et al.
#'
#'
#' @param fn Test function
#' @param X Matrix of input data
#' @param distr Data distribution
#' @param grads  Matrix of gradients
#' @param num_obj Number of objectives
#' @param norm_params List containing the mean and the standard deviation of Y (ground truth)
#' @param num_orth_sample Sample size for the orthogonal projection
#' @param seed Random seed
#' @param params List of additional parameters to be passed to the function (to be specified in a list named \code{params})
#'
#' @return Matrix of function estimations for all dimensions
#' @noRd
#' @references
#' Zahm, O., P. G. Constantine, C. Prieur, and Y. M. Marzouk (2020). Gradient-based dimension reduction of multivariate vector-valued functions. SIAM Journal on Scientific Computing 42 (1), A534–A558
Zahm <- function(fn, X, distr, grads, num_obj, norm_params, num_orth_sample = 20, seed=16, params) {

  d <- ncol(X); n <- distr$n

  Result <- list()

  H_hat <- matrix(0, nrow=d, ncol=d) # represents the sum of covariance matrices

  Sigma <- if(distr$name=="norm") distr$sigma else diag(2/sqrt(d), d)

  #calculate the sum of covariance matrices
  for(i in 1:num_obj) {

    g <- grads[,((i-1)*d+1):(i*d)]

    H_hat <- H_hat + t(g) %*% g/n

    }

  Sigma_inv <- if(distr$name=="norm") inv(Sigma) else Sigma

  H_hat_sigma_geigen <- geigen(H_hat, Sigma_inv, symmetric=TRUE)

  idx <- order(H_hat_sigma_geigen$values, decreasing = T)

  H_hat_sigma_geigen$vectors <- H_hat_sigma_geigen$vectors[,idx]

  I <- diag(1, d)

  if (exists(".Random.seed", .GlobalEnv))  oldseed <- .GlobalEnv$.Random.seed

  else oldseed <- NULL

  set.seed(seed)

  #sample_orth <- rmvnorm(n = num_orth_sample, mean = rep(0,d), sigma = Sigma)
  sample_orth <- if(distr$name=="unif") matrix(runif(d*num_orth_sample, min = -1, max = 1), num_orth_sample) else rmvnorm(n=num_orth_sample, mean=rep(0,d), sigma=Sigma)

  for(i in 1:d) {

    P_r <- H_hat_sigma_geigen$vectors[,1:i] %*% t(H_hat_sigma_geigen$vectors[,1:i]) %*% Sigma_inv

    Yest <- 0

    for(iM in 1:num_orth_sample) {

      b <- t((I-P_r)%*%sample_orth[iM,])

      X_approx_gauss <- t(apply(X, 1,  function(x) x%*%P_r + b))

      Yest <- Yest + wrapper(X_approx_gauss, fn, which_distr = distr$name, norm_params, params)

       }

      Yest <- Yest/num_orth_sample

      Result[[i]] <- Yest

     }

    if (!is.null(oldseed)) .GlobalEnv$.Random.seed <- oldseed

    else rm(".Random.seed", envir = .GlobalEnv)

    return(Result)

  }

