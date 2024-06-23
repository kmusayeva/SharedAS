#' Sufficient summary plots for Ptest
#'
#' @param fn Test function
#' @param X Input data
#' @param Y Output data
#' @param distr Name of the distribution
#' @param evectors Matrix of eigenvectors
#' @param methods List of methods
#' @param width Width of the plot
#' @param height Height of the plot
#' @param Yout Polygon representing the original image of the test function
#' @noRd
#'
#' @examples
#'
#' fn <- problem_synthetic
#' n <- 1000
#' L <- handle_data(fn = fn, X = NULL, input_dim = 3, num_obj = 2, distr = list(name = "unif", n = 1000))
#' grads <- compute_grad(fn, L$X, which_distr = "unif", norm_params = L$norm_params)
#' methods <- c("AG", "LP", "MCH", "SSPD", "SEE", "FG")
#' evectors <- compute_evectors(methods, grads, n = n, m = 2, d = 3)
#' coordspol <- coords_poly()
#'
#' summary_plot(fn, L$X, L$Y, which_distr = "unif", evectors = evectors, methods = methods, width = 9.5, height = 5, coordspol)
#'
summary_plot <- function(fn, X, Y, which_distr, evectors, methods, width, height, Yout) {
  for (im in 1:length(evectors)) {
    pdf(file = paste0("./figures/", methods[im], ".pdf"), width = width, height = height)

    par(mfrow = c(1, 3))

    W1 <- evectors[[im]][, 1]

    W12 <- evectors[[im]][, 1:2]

    AS1 <- X %*% W1

    AS2 <- X %*% W12

    Yest1 <- wrapper(X %*% W1 %*% t(W1), fn, which_distr)

    Yest12 <- wrapper(X %*% W12 %*% t(W12), fn, which_distr)

    ### f1
    plot(AS1, Y[, 1], xlab = "AS1", ylab = expression(f[1]), cex = 0.3, col = gray(0.5), cex.lab = 1.8, cex.axis = 1.2, xlim = c(-3.2, 4), ylim = c(-13, -6.8))
    points(AS1, Yest1[, 1], xlab = "AS1", ylab = expression(f[1]), cex = 0.2, col = "red")

    ### f2
    plot(AS1, Y[, 2], xlab = "AS1", ylab = expression(f[2]), cex = 0.3, col = gray(0.5), cex.lab = 1.8, cex.axis = 1.2, xlim = c(-3.2, 4), ylim = c(20, 52))
    points(AS1, Yest1[, 2], xlab = "AS1", ylab = expression(f[2]), cex = 0.2, col = "red")


    plot(Yout[, 1], Yout[, 2], type = "n", xlab = expression(f[1]), ylab = expression(f[2]), cex.lab = 1.8, cex.axis = 1.2)

    polygon(x = Yout[, 1], y = Yout[, 2], col = rgb(0, 0, 0.5, 0.1), border = rgb(0, 0, 0.5, 0.1))

    points(Yest12[, 1], Yest12[, 2], col = "green", cex = 0.3)

    points(Yest1[, 1], Yest1[, 2], col = "red", cex = 0.3)

    legend(-9.5, 53,
      legend = c("Original", "1d AS", "2d AS"), col = c(rgb(0, 0, 0.5, 0.1), "red", "green"), pch = rep(19, 3), cex = 1.4,
      box.lty = 0, y.intersp = 1.1
    )

    mtext(paste0(methods[im]), side = 3, line = -3, outer = TRUE)

    dev.off()
  }
}



#' Boundary points for the image space of ptest function
#' @noRd
coords_poly <- function() {
  library(alphahull)
  # Get the envelope of point cloud
  n <- 1e5
  d <- 3
  set.seed(19)
  A <- apply(RRembo::selectA(3, 3), c(1, 2), signif, digits = 2)
  B <- apply(RRembo::selectA(3, 3), c(1, 2), signif, digits = 2)
  X <- matrix(runif(n * d), n)

  params <- list(A = A, B = B, seed = 19)

  Y <- t(apply(X, 1, problem_synthetic, params))
  # plot(Y)

  alpha <- 2
  shap <- ashape(Y, alpha = alpha)
  # plot(shap)
  plot(Y[shap$alpha.extremes, ]) # get the index of the points on the contour

  ## Hard coded version of the contour (keeping the seed above)
  # All ids of extreme points
  ids0 <- c(
    77381, 7837, 4458, 95484, 42569, 83219, 38600, 75734, 53170, 62702, 57568, 7620, 76192, 98391, 58887, 90857, 52245, 649, 59481, 42699, 36679, 40649, 39660,
    3960, 98057, 11072, 40938, 96494, 33900, 98940, 76353, 95734, 1882, 25315, 69237, 51341, 134, 553, 608, 626, 1533, 1814, 1961, 2189, 2689, 2885,
    3443, 3578, 3934, 4381, 5499, 5623, 5870, 5937, 6323, 8142, 8205, 8234, 10279, 11774, 11781, 11961, 12126, 12625, 13471, 14698, 15733, 15837, 16280,
    16707, 16848, 17071, 17176, 17594, 17670, 18059, 19316, 19475, 20284, 20367, 20391, 20621, 21382, 22152, 22511, 24018, 24791, 25382, 25538, 26342, 26351, 27268,
    27665, 29850, 30708, 30958, 31030, 31359, 32537, 34061, 34231, 34682, 36517, 37447, 38653, 38675, 38731, 39398, 39717, 39797, 40470, 40599, 40619, 41323, 41759,
    42327, 42846, 43301, 43421, 44071, 45349, 46941, 47163, 47855, 48041, 48252, 49482, 50005, 50134, 50444, 50859, 51347, 52427, 52787, 56191, 56875, 56973, 57199,
    58112, 58180, 58200, 58433, 59103, 59784, 60106, 60690, 60980, 61029, 61501, 62025, 62055, 62839, 62989, 63609, 64529, 66093, 67040, 67537, 68179, 68339, 69639,
    70199, 70252, 70679, 70743, 70966, 71316, 72019, 73933, 76049, 76154, 76985, 77098, 77138, 77761, 77781, 77963, 78671, 79557, 79579, 80139, 80492, 80609, 81028,
    81080, 81100, 83392, 83685, 84417, 84904, 87894, 88188, 88311, 88653, 88744, 89079, 90216, 90412, 91243, 91782, 92031, 92094, 92253, 93937, 94062, 94181, 94274,
    94869, 95124, 96142, 96415, 96451, 96506, 96701, 98847, 99202, 99390
  )

  ids <- ids0[order(order(Y[ids0, 1]))]

  ids1 <- which(Y[ids, 1] < -9 & Y[ids, 1] > -12 & Y[ids, 2] < 30)
  idstmp <- ids[-ids1]
  plot(Y[ids, ])
  points(Y[ids[ids1], ], col = "red", pch = 20)
  Yout <- Y[ids[ids1], ]
  Yout <- Yout[order(Yout[, 1]), ]

  ids2 <- which(Y[idstmp, 1] > -9 & Y[idstmp, 2] < 23)
  points(Y[idstmp[ids2], ], col = "green", pch = 20)
  Ytmp <- Y[idstmp[ids2], ]
  Ytmp <- Ytmp[order(Ytmp[, 1]), ]
  Ytmp <- Ytmp[c(1:9, 11, 12, 14, 15, 13, 10), ]
  Yout <- rbind(Yout, Ytmp)
  idstmp <- idstmp[-ids2]

  ids3 <- which(Y[idstmp, 1] > -11.4)
  ids3 <- ids3[order(order(Y[ids3, 2]))]
  points(Y[idstmp[ids3], ], col = "orange", pch = 20)
  Ytmp <- Y[idstmp[ids3], ]
  Yout <- rbind(Yout, Ytmp[order(Ytmp[, 2]), ])
  idstmp <- idstmp[-ids3]

  ids4 <- which(Y[idstmp, 2] > 50)
  points(Y[idstmp[ids4], ], col = "blue", pch = 20)
  Ytmp <- Y[idstmp[ids4], ]
  Yout <- rbind(Yout, Ytmp[order(-Ytmp[, 1]), ])
  idstmp <- idstmp[-ids4]

  ids5 <- idstmp[order(order(-Y[idstmp, 2]))]
  points(Y[ids5, ], col = "yellow", pch = 20)
  Ytmp <- Y[ids5, ]
  Yout <- rbind(Yout, Ytmp[order(-Ytmp[, 2]), ])

  return(Yout)
}
