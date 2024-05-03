#' Car side-impact problem
#'
#' @param x a vector in [0, 1]^7 corresponding to
#' 1: Thickness of B-Pillar inner
#' 2: Thickness of B-Pillar reinforcement
#' 3: Thickness of floor side inner
#' 4: Thickness of cross members
#' 5: Thickness of door beam
#' 6 : Thickness of door beltline reinforcement
#' 7 : Thickness of roof rail
#' @param params List of additional parameters:
#' \itemize{
#' \item {
#' \code{p} optional vector of values for the random parameter, if \code{fixed = TRUE}, no sampling is done and they are fixed at their mean:
#' 8: Material of B-Pillar inner
#' 9: Material of floor side inner
#' 10: Barrier height
#' 11: Barrier hitting position}
#' \item{\code{mp,sp} mean and standard deviation for p}
#' \item{\code{fixed} Logical}
#'}
#'
#' @return Vector of size 11

#' @export
#'
#' @references
#' Deb, K. et al. “Reliability-Based Optimization Using Evolutionary Algorithms.”
#' Evolutionary Computation, IEEE Transactions on 13.5 (2009): 1054-1074. 2009 Institute of
#' Electrical and Electronics Engineers
#'
#' @examples
#' d <- 7
#' n <- 1000
#' X <- matrix(runif(d*n), n)
#' Y <- t(apply(X, 1, problem_car_side_impact))
#'
#' pairs(Y)
#'
problem_car_side_impact <- function(x, params=list(p = NULL, mp = c(0.345, 0.192, 0, 0), sp = c(0.06, 0.06, 10, 10), fixed = TRUE)) {
  p <- params$p
  mp <- params$mp
  sp <- params$sp
  fixed <- params$fixed

  lower <- c(0.5, 0.45, 0.5, 0.5, 0.875, 0.4, 0.4)
  upper <- c(1.5, 1.35, 1.5, 1.5, 2.625, 1.2, 1.2)

  x <- x * (upper - lower) + lower

  if(fixed) p <- mp else p <- rnorm(4, mean = mp, sd = sp)

  x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; x5 <- x[5]; x6 <- x[6]; x7 <- x[7]
  x8 <- p[1]; x9 <- p[2]; x10 <- p[3]; x11 <- p[4]


  f1 <- 1.98 + 4.9* x1 + 6.67*x2 + 6.98*x3  + 4.01*x4 + 1.78*x5 + 0.00001*x6 + 2.73*x7

  g1 <- 1.16 - 0.3717 * x2 * x4 - 0.00931* x2 * x10 - 0.484 *x3 *x9 + 0.01343 * x6 * x10

  g2 <- 0.261 - 0.0159 * x1 * x2- 0.188 * x1 * x8 - 0.019 * x2 * x7  +
    0.0144 * x3 * x5 + 0.87570001 * x5 * x10 + 0.08045 * x6 * x9 +
    0.00139 * x8 * x11 + 0.00001575 * x10 * x11

  g3 <- 0.214 + 0.00817 * x5 - 0.131 * x1 * x8 - 0.0704 * x1 * x9 + 0.03099 * x2 * x6 -
    0.018 * x2 * x7 + 0.0208 * x3 * x8 + 0.121 * x3 * x9 - 0.00364 * x5 * x6 + 0.0007715 * x5 * x10 -
    0.0005354 * x6 * x10 + 0.00121 * x8 * x11 + 0.00184 * x9 * x10 - 0.018 * x2 * x2

  g4 <- 0.74 - 0.61 * x2 - 0.163* x3 * x8 + 0.001232 *x3 * x10 - 0.166 * x7 * x9 + 0.227  *x2 *x2

  g5 <- 28.98 + 3.818 * x3 - 4.2 * x1 * x2 + 0.0207 * x5 * x10 + 6.63 * x6 * x9 - 7.77 * x7 * x8 + 0.32 * x9 * x10

  g6 <- 33.86 + 2.95 * x3  + 0.1792 * x10 - 5.057 * x1 * x2 - 11 * x2 * x8 - 0.0215 * x5 * x10 - 9.98 * x7 * x8 + 22 *x8 * x9

  g7 <- 46.36 - 9.9 * x2 - 12.9 * x1 * x8  + 0.1107 * x3 * x10

  g8 <- 4.72 - 0.5 * x4 - 0.19 * x2 * x3 - 0.0122 * x4 * x10 + 0.009325 * x6 * x10 + 0.000191 * x11 *  x11

  g9 <- 10.58 - 0.674 * x1 * x2 - 1.95 * x2 * x8 + 0.02054 * x3 * x10 - 0.0198 * x4 * x10 + 0.028 * x6 * x10

  g10 <- 16.45 - 0.489 * x3 * x7 - 0.843 * x5 * x6 + 0.0432 * x9 * x10 - 0.0556 * x9 * x11 - 0.000786 * x11 * x11

  return(c(f1, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10))
}

