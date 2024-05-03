#' Penicillin production test function based on the implementation github.com/HarryQL/TuRBO-Penicillin.
#' @param input Matrix of input data
#' @param params List of additional parameters:
#' \itemize{
#' \item \code{t} Reaction time in hours, defaulted to 100
#' \item \code{returnCST} Optional. If \code{TRUE} returns the constraints as part of the output.
#' }
#'
#' @return Vector of objective values: if \code{returnCST} is set to \code{TRUE} a vector of size 5 where the last 3 columns correspond to the constraints.
#'
#' @export
#' @references
#' Liang, Q. and L. Lai (2021). Scalable bayesian optimization accelerates process optimization of penicillin production. In NeurIPS 2021 AI for Science Workshop.
problem_penicillin <- function(input, params=list(t=100, returnCST=FALSE)) {

  lower <-c(60, 0.01, 293, 0.01, 0.01, 400, 5)
  upper <-c(120, 12, 303, 18, 0.5, 600, 7.5)

  input <- input * (upper - lower) + lower

  t <- params$t

  V <- input[1]
  X <- input[2]
  T_ <- input[3]
  S <- input[4]
  F_ <- input[5]
  s_f <- input[6]
  H_ <- input[7]


  ### params
  CO2 <- 0
  P <- 0
  H <- 10^(-H_)

  Y_xs <- 0.45
  Y_ps <- 0.90

  K_1 <- 10^(-10)
  K_2 <- 7 * 10^(-5)
  m_X <- 0.014

  alpha_1 <- 0.143
  alpha_2 <- 4 * 10^(-7)
  alpha_3 <- 10^(-4)
  mu_X <- 0.092
  K_X <- 0.15


  mu_p <- 0.005
  K_p <- 0.0002
  K_I <- 0.10
  p_K <- 0.04
  k_g <- 7 * 10^(3)
  E_g <- 5100
  k_d <- 10^(33)
  E_d <- 50000

  alpha <- 70
  lambd <- 2.5 * 10^(-4)

  T_v <- 273
  T_o <- 373

  p_R <- 1.9872
  V_max <- 180


  for (i in 1:t) {

    F_loss <- V * lambd * (exp(5 * ((T_ - T_o) / (T_v - T_o))) - 1)

    dV_dt <- F_ - F_loss

    mu <- (mu_X / (1 + K_1 / H + H / K_2)) * (S / (K_X * X + S)) * ((k_g * exp(-E_g / (p_R * T_))) - (k_d * exp(-E_d / (p_R * T_))))

    dX_dt <- mu * X - (X / V) * dV_dt

    mu_pp <- mu_p * (S / (K_p + S + S^2 / K_I))

    dS_dt <- - (mu / Y_xs) * X - (mu_pp / Y_ps) * X - m_X * X + F_ * s_f / V - (S / V) * dV_dt

    dP_dt <- (mu_pp * X) - p_K * P - (P / V) * dV_dt

    dCO2_dt <- alpha_1 * dX_dt + alpha_2 * X + alpha_3

    # UPDATE
    P <- P + dP_dt
    V <- V + dV_dt
    X <- X + dX_dt
    S <- S + dS_dt
    CO2 <- CO2 + dCO2_dt
    l_t <- i

    if (V > V_max) {
      # print("Too large V")
      break
    }

    if (S < 0) {
      # print("Too small S")
      break
    }

    if (dP_dt < 10e-12) {
      # print("Converged P")
      break
    }

  }

  if(params$returnCST) return(c(l_t, P, CO2, V, X, S))

  return(c(l_t, P, CO2))

  }

