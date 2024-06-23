#' Computes shared active subspace for vector-valued functions.
#'
#' @param fn Vector-valued function for which to compute a shared active subspace
#' @param grads Optional, matrix of objective-wise gradients concatenated together
#' @param X Matrix of input data: to provide only if \code{grads} are provided
#' @param distr List indicating the data distribution and its parameters:
#' \itemize{
#' \item \code{name} either "unif" for the uniform or "norm" for a normal (default is the standard normal) distribution
#' \item \code{n} number of samples
#' \item \code{sigma} covariance matrix, if the distribution is set to "norm" (the mean is always zero)
#'  }
#' @param input_dim Dimensionality of the input data
#' @param num_obj Number of outputs of the test function
#' @param methods Vector of method names to compute a shared subspace
#' \itemize{
#' \item AG: average of gradients
#' \item LP: linear projection of gradients
#' \item MCH: minimum convex hull of gradients
#' \item SEE: stepwise estimation of eigenvectors
#' \item SSPD: sum of SPD matrices
#' \item FG: joint diagonalization of SPD matrices
#' \item Zahm: the method of Zahm
#' }
#' @details Depending on the problem, the method of MCH can be costly to compute.
#' @param nexps Number of experiments
#' @param is_norm Logical, indicates whether to normalize the function outputs in the gradient computation (default is TRUE)
#' @param start_dim Optional, start number of the dominating eigenvectors
#' @param end_dim End number of the dominating eigenvectors: To be set if \code{start_dim} is set.
#' @param seed Random seed for reproducibility
#' @param params List of additional parameters for the test function
#'
#' @return List of objective-wise root-mean-square error and its sum for each number of dimensions and for each method, over \code{nexps} number of experiments
#' @export
#'
#' @examples
#'
#' fn <- problem_car_side_impact
#' input_dim <- 7
#' num_obj <- 11
#' nexps <- 5
#' # distr <- list(name="unif", n=10)
#' distr <- list(name = "norm", n = 1000)
#' SharedAS(
#'   fn = fn, distr = distr, input_dim = input_dim, num_obj = num_obj,
#'   methods = c("AG", "LP", "MCH", "SSPD", "SEE", "FG", "Zahm"),
#'   nexps = nexps, is_norm = TRUE, seed = 126, params = NULL
#' )
#'
#'
#' fn <- problem_marinedes
#' input_dim <- 6
#' num_obj <- 3
#' nexps <- 1
#' distr <- list(name = "norm", n = 1000)
#' SharedAS(
#'   fn = fn, distr = distr, input_dim = input_dim, num_obj = num_obj,
#'   methods = c("AG", "LP", "MCH", "SSPD", "SEE", "FG", "Zahm"),
#'   nexps = nexps, is_norm = TRUE, seed = 126, params = NULL
#' )
#'
#'
#' fn <- problem_penicillin
#' input_dim <- 7
#' num_obj <- 5
#' nexps <- 1
#' distr <- list(name = "norm", n = 1000)
#' params <- list(t = 100, returnCST = TRUE)
#' SharedAS(
#'   fn = fn, distr = distr, input_dim = input_dim, num_obj = num_obj,
#'   methods = c("AG", "LP", "MCH", "SSPD", "SEE", "FG", "Zahm"), nexps = nexps,
#'   is_norm = TRUE, seed = 126, params = params
#' )
#'
#'
#' fn <- problem_switch_ripple
#' input_dim <- 8
#' num_obj <- 5 # if returCST=TRUE, then numb_obj=10
#' distr <- list(name = "norm", n = 1000)
#' nexps <- 1
#' params <- list(returnCST = FALSE)
#' SharedAS(
#'   fn = fn, distr = distr, input_dim = input_dim, num_obj = num_obj,
#'   methods = c("AG", "LP", "MCH", "SSPD", "SEE", "FG", "Zahm"), nexps = nexps,
#'   is_norm = TRUE, seed = 126, params = params
#' )
#'
#'
#' fn <- problem_synthetic
#' input_dim <- 3
#' num_obj <- 2
#' nexps <- 30
#' distr <- list(name = "unif", n = 1000)
#' SharedAS(
#'   fn = fn, distr = distr, input_dim = input_dim, num_obj = num_obj,
#'   methods = c("AG", "LP", "MCH", "SSPD", "SEE", "FG", "Zahm"), nexps = nexps,
#'   is_norm = TRUE, seed = 126
#' )
#'
SharedAS <- function(fn, grads = NULL, X = NULL, distr = list(name = "norm", n = 1000), input_dim, num_obj, start_dim = NULL, end_dim = NULL, methods = NULL, nexps = 1, is_norm = TRUE, seed = 126, params = NULL) {
  all_methods <- c("AG", "LP", "MCH", "SSPD", "SEE", "FG", "Zahm")

  if (is.null(methods)) {
    methods <- all_methods
  } else if (sum(!methods %in% all_methods)) stop(cat("Possible choices: ", all_methods, "\n"))

  # if(distr$name=="unif" && "Zahm" %in% methods) {
  #
  #   methods <- methods[which(methods!="Zahm")]
  #
  #   if(length(methods)==0) stop("The method of Zahm applies only to a normal distribution. \n")
  #
  #   cat("The method of Zahm applies only to a normal distribution.\n")
  #
  #   }

  if (!is.null(grads) && is.null(X)) stop("Please specify the X matrix.")

  if (!distr$name %in% c("norm", "unif")) stop("Possible choice for distribution: 'norm' or 'unif'.")

  if (is.null(distr$n)) stop("Please specify the number of samples.")

  if (is.null(input_dim)) stop("Please specify the number of input dimensions.")

  if (is.null(num_obj)) stop("Please specify the number of function outputs.")

  if (distr$name == "norm" && is.null(distr$sigma)) distr$sigma <- diag(1, input_dim)

  if (!is.null(grads) && !is.null(X)) nexps <- 1

  if (is.null(start_dim)) start_dim <- 1

  if (is.null(end_dim)) end_dim <- input_dim

  if (end_dim < start_dim || end_dim > input_dim) stop("Start or end dimension is incorrect.")

  Result <- compute(fn, grads, X, input_dim, num_obj, distr, methods, nexps, is_norm, params, start_dim, end_dim, seed)

  # print(Result)

  saveRDS(Result, paste0(format(Sys.time(), "%X"), ".rds"))
}


#' Computes the eigenspace, rotates the data and evaluates the test function
#'
#' @param fn Vector-valued test function
#' @param grads Matrix of gradients
#' @param X Matrix of input data, optional
#' @param input_dim Dimensionality of the input space
#' @param num_obj Number of function outputs
#' @param distr Data distribution
#' @param methods Vector of method names
#' @param nexps Number of experiments
#' @param is_norm If \code{TRUE}, normalizes when computing the gradients
#' @param seed random seed for reproducibility
#'
#' @return Matrix of RMSEs and their sum for each experiment, dimension and method
#' @noRd
compute <- function(fn, grads, X_given, input_dim, num_obj, distr, methods, nexps, is_norm, params, start_dim, end_dim, seed) {
  which_distr <- distr$name

  Result <- create_Result_matrix(start_dim, end_dim, num_obj, nexps, methods)

  set.seed(seed)

  for (nexp in 1:nexps) {
    cat("Sample: ", nexp, "\n")

    # data part
    dat <- handle_data(fn, X_given, input_dim, num_obj, distr, params)

    X <- dat$X
    d <- ncol(X)
    Y <- dat$Y
    Y_norm <- dat$Y_norm
    num_obj <- dat$num_obj
    n <- nrow(X)

    norm_params <- dat$norm_params

    # grads part
    cat("computing gradients...")

    grads1 <- if (is.null(grads)) handle_grads(fn, X, which_distr, is_norm, norm_params, params) else grads

    cat("done \n")

    # eigenspace part
    cat("computing eigenspaces... ")

    evectors <- compute_evectors(methods, grads1, n, num_obj, d)

    cat("done\n")

    # evaluation part
    Result <- evaluate(Result, nexp, fn, X, grads1, distr, Y_norm, num_obj, start_dim, end_dim, evectors, norm_params, params)
  }

  return(Result)
}




#' Given the eigenvectors, reconstructs the data and evaluates the test function on it.
#'
#' @param Result Matrix of rmses
#' @param nexp The number of iteration in the experiment
#' @param fn Test function
#' @param X Matrix of input data
#' @param grads Matrix of gradients
#' @param distr Data distribution
#' @param start_dim Start dimensionality
#' @param end_dim End dimensionality
#' @param evectors Matrix of eigenvectors
#' @param Y_norm Normalized Y matrix
#' @param norm_params Normalization parameters for the estimated Y values
#' @param params List of additional parameters to the test function
#'
#' @return Matrix of rmses and their sum for each dimension and for each method
#' @noRd
evaluate <- function(Result, nexp, fn, X, grads, distr, Y_norm, num_obj, start_dim, end_dim, evectors, norm_params, params) {
  # if specified, execute the method of Zahm
  if (length(Result$Zahm)) {
    cat("executing the method of Zahm... ")

    # start_time <- Sys.time()
    res_zahm <- Zahm(fn, X, distr, grads, num_obj = num_obj, norm_params, num_orth_sample = 20, seed = 16, params)
    # end_time <- Sys.time()
    # cat("Zahm runtime:", round(end_time - start_time, 4))
    cat("done\n")
  }

  cat("evaluating... ")

  for (num_dim in start_dim:end_dim) {
    if (length(evectors)) {
      for (ei in 1:length(evectors)) {
        res <- project_and_estimate(evectors[[ei]], fn, X, num_dim, which_distr = distr$name, norm_params, params)

        rmses <- rmse_wrapper(Y_norm, res$Yest, params)

        Result[[names(evectors[ei])]][[num_dim]][nexp, ] <- rmses
      }
    }

    if (length(Result$Zahm)) {
      rmses <- rmse_wrapper(Y_norm, res_zahm[[num_dim]], params)

      Result$Zahm[[num_dim]][nexp, ] <- rmses
    }
  }


  cat("done\n")
  return(Result)
}
