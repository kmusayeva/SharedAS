#' Switching ripple test function.
#' Corresponds to the switching ripple suppressor design problem for voltage source inversion in powered system.
#' @param x vector specifying the location where the function is to be evaluated, of size k + 4, k > 1, see Details.
#' @param params List of additional parameters:
#' \itemize{
#' \item \code{returnCST} Optional. If set to \code{TRUE} the last 5 columns in the output correspond to the values of the constraints.
#' }
#' @details
#' Columns of \code{x} correspond to L1, L2, L3, C1, ..., Ck, Cf where k is an arbitrary integer > 1.
#'
#' Parameters of the problems follow Table 2 in (Zhang et al, 2019).
#'
#' @return Vector of objectives. The first \code{k} values are related to the suppression of harmonics while the \code{k+1} one is the sum of inductors.
#' @export
#' @references
#' \itemize{
#' \item{ Zhang, Z., He, C., Ye, J., Xu, J., & Pan, L. (2019). Switching ripple suppressor design of the grid-connected inverters: A perspective of many-objective optimization with constraints handling. Swarm and evolutionary computation, 44, 293-303.}
#' \item{ He, C., Tian, Y., Wang, H., & Jin, Y. (2019). A repository of real-world datasets for data-driven evolutionary multiobjective optimization. Complex & Intelligent Systems, 1-9.}
#' }
#'
problem_switch_ripple <- function(x, params=list(returnCST = FALSE)) {

  lower <- c(0.000110815602836879, 7.83532027529331e-06, 1.29313391262165e-06, rep(1.29313391262165e-06, length(x)-3))
  upper <- c(0.000221631205673759, 0.000783532027529331, 0.000783532027529331, rep(6.46566956310825e-05, length(x)-3))

  x <- x * (upper - lower) + lower
  d <- length(x) # number of variables

  nr <- length(x) - 4 # n in the paper

  ## Define parameters
  Ell <- 400
  Prated <- 65000
  Vdc <- 800
  omega0 <- 100 * pi
  omegasw <- 32000 * pi
  Iref <- 141
  Rlk <- rep(0.005, nr)
  Rb <- Ell^2/Prated

  L1 <- x[1]
  L2 <- x[2]
  L3 <- x[3]
  Ck <- x[4:(d-1)]
  Cf <- x[d]

  Lk <- 1 / (Ck * (1:nr * omegasw)^2)

  GLC <- function(s){

    res <- (L2*s * (L3 * Cf * s^2 + 1) + L3 * s) / (Lk * s + Rlk + 1 / (Ck * s))

    sum(res)

  }

  G <- function(s){
    1 / (L1 * s * GLC(s) + L2 * s * (L3 * Cf * s^2 + 1) + L3 * s)
  }

  fi <- function(i) 20 * log(Mod(G(i * omegasw * 1i)))

  fs <- sapply(1:nr, fi)

  fM <- L1 + L2 + L3 + sum(Lk)

  if(params$returnCST){

    g1 <- sum(Ck) + Cf - 0.05/(Rb * omega0)

    g2 <- L1 + L2 + L3 - 0.1 * Rb / omega0

    # g31 <- 0.2 - 2 * pi * Vdc / (8 * L1 * omegasw * Iref) # inactive returnCST
    # g32 <- 2 * pi* Vdc / (8 * L1 * omegasw * Iref) - 0.4

    g4 <- L2 + L3 - L1

    g51 <- omegasw/2 - sqrt((L1 + L2 + L3) / (L1 * (L2 + L3) * (sum(Ck) + Cf)))

    g52 <- sqrt((L1 + L2 + L3)/(L1 * (L2 + L3) * (sum(Ck) + Cf))) - 3*omegasw / 4

    return(c(fs, fM, g1, g2, g4, g51, g52))

    }

  else {

    return(c(fs, fM))

    }

}

