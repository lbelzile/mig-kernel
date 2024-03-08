################################################################
## Local limit theorems for the inverse Gaussian distribution ##
################################################################

## Code written by Frederic Ouimet (December 2022)

require("statmod") # Inverse Gaussian distribution
require("mvtnorm") # Multivariate normal distribution
require("latex2exp") # to write LaTeX symbols

##########################
# Set the path for plots #
##########################
setwd(this.path::here())

##############
# Parameters #
##############

d <- 2

low_mu <- 10
up_mu <- 1000
inc_mu <- 10
ratio_plot_points <- 10

mus <- seq(low_mu, up_mu, inc_mu)
len_mu <- length(mus)

#############
# Functions #
#############

delta <- function (x, mu, beta, xi0, Omega) {
  return (as.numeric((mu * (t(
    beta
  ) %*% xi0)) ^ (-1 / 2)) * solve(t(chol(Omega)), x - mu * xi0))
}

log_ratio_error <- function (x, mu, beta, xi0, Omega, ord) {
  d <- length(x)
  del <- as.matrix(delta(x, mu, beta, xi0, Omega))
  if (ord == 0) {
    # Order 0 correction
    return (log((2 * pi) ^ (-d / 2) * mu * t(beta) %*% xi0 * det(Omega) ^ (-1 / 2) * (t(beta) %*% x) ^ (-(d / 2 + 1)) * exp(-1 / (2 * t(beta) %*% x) * t(x - mu * xi0) %*% solve(Omega, x - mu * xi0)) / mvtnorm::dmvnorm(
      t(x),
      mean = t(mu * xi0),
      sigma = as.numeric(mu * (t(beta) %*% xi0)) * Omega
    )
    ))
  } else if (ord == 1) {
    # Order 1 correction
    return (
      log((2 * pi) ^ (-d / 2) * mu * t(beta) %*% xi0 * det(Omega) ^ (-1 / 2) * (t(beta) %*% x) ^ (-(d / 2 + 1)) * exp(-1 / (2 * t(beta) %*% x) * t(x - mu * xi0) %*% solve(Omega, x - mu * xi0)) / mvtnorm::dmvnorm(
        t(x),
        mean = t(mu * xi0),
        sigma = as.numeric(mu * (t(beta) %*% xi0)) * Omega
      )
      ) + ((d + 2) / 2 - t(del) %*% del / 2) * (t(beta) %*% t(chol(Omega)) %*% del / sqrt(mu * t(beta) %*% xi0))
    )
  } else {
    # ord = 2 : Order 2 correction
    return (
      log((2 * pi) ^ (-d / 2) * mu * t(beta) %*% xi0 * det(Omega) ^ (-1 / 2) * (t(beta) %*% x) ^ (-(d / 2 + 1)) * exp(-1 / (2 * t(beta) %*% x) * t(x - mu * xi0) %*% solve(Omega, x - mu * xi0)) / mvtnorm::dmvnorm(
        t(x),
        mean = t(mu * xi0),
        sigma = as.numeric(mu * (t(beta) %*% xi0)) * Omega
      )
      ) + ((d + 2) / 2 - t(del) %*% del / 2) * (t(beta) %*% t(chol(Omega)) %*% del / sqrt(mu * t(beta) %*% xi0)) - ((d + 2) / 4 - t(del) %*% del / 2) * (t(beta) %*% t(chol(Omega)) %*% del / sqrt(mu * t(beta) %*% xi0)) ^ 2
    )
  }
}

#############
# Main code #
#############

# The choice of beta and xi0 below implies t(beta) %*% xi0 = 1
beta = as.matrix(c(1 / 2, 1 / 2))
xi0 = as.matrix(c(1, 1))

for (a1 in 2:5) {
  for (a2 in 2:5) {
    for (bb in 1:1) {

      start_time <- Sys.time()

      Omega = matrix(c(a1, bb, bb, a2), nrow = d, ncol = d)

      max_error_vec <- array(NaN, dim = c(len_mu, 3))

      for (i in 1:len_mu) {
        musi <- mus[i]
        lst <- vector(mode = "list", length = d)
        for (k in 1:d) {
          lst[[k]] <- seq(musi - sqrt(musi), musi + sqrt(musi), 1)
        }
        hypercube <-
          expand.grid(lst) # hypercube of d dimensions centered around r
        len_hyper <-
          dim(hypercube)[1] # there was an error here in the previous version of the code
        error_mat <- array(NaN, dim = c(len_mu, len_hyper, 3))
        for (j in 1:len_hyper) {
          x_hyp_j <- as.numeric(hypercube[j,])
          error_mat[i, j, 1] <-
            abs(log_ratio_error(as.matrix(x_hyp_j), musi, beta, xi0, Omega, 0))
          error_mat[i, j, 2] <-
            abs(log_ratio_error(as.matrix(x_hyp_j), musi, beta, xi0, Omega, 1))
          error_mat[i, j, 3] <-
            abs(log_ratio_error(as.matrix(x_hyp_j), musi, beta, xi0, Omega, 2))
        }
        max_error_vec[i, 1] <- max(error_mat[i, , 1])
        max_error_vec[i, 2] <- max(error_mat[i, , 2])
        max_error_vec[i, 3] <- max(error_mat[i, , 3])
      }

      mypath <-
        paste0("../_manuscript/figures/loglog_errors_plot_",
               a1, "_", a2, "_", bb, ".pdf")
      pdf(file = mypath,
          width = 5,
          height = 5)
      par(
        bg = "white",
        # mgp = c(2, 0.8, 0),
        mar = c(3.2,3.2,1,1),
        cex = 1.3
      )
      mus_plot_points <-
        seq(low_mu, up_mu, inc_mu * ratio_plot_points)
      mus_index_plot_points <- seq(1, len_mu, ratio_plot_points)
      plot(
        mus_plot_points,
        1 / max_error_vec[mus_index_plot_points, 1],
        log = "xy",
        ylim = c(1, 1 / min(max_error_vec[, 3])),
        col = 4,
        pch = 15,
        xlab = "",
        ylab = "",
        cex.lab = 1.2,
        cex.axis = 1.2,
        yaxt = "n"
      )
      axis(2, at = c(1e1, 1e2, 1e3,1e4, 1e5),
           cex = 1.2,
           cex.axis = 1.2,
           labels = c(10, "",expression(10^3), "", expression(10^5))) #c(1, 10, expression(10^2), expression(10^3), expression(10^4),                      expression(10^5)))
      points(mus_plot_points,
             1 / max_error_vec[mus_index_plot_points, 2],
             col = 3,
             pch = 16)
      points(mus_plot_points,
             1 / max_error_vec[mus_index_plot_points, 3],
             col = 2,
             pch = 17)
      lines(mus, 1 / max_error_vec[, 1], col = 4, lwd = 2)
      lines(mus, 1 / max_error_vec[, 2], col = 3, lwd = 2)
      lines(mus, 1 / max_error_vec[, 3], col = 2, lwd = 2)
      legend(
        "topleft",
        legend = c(TeX("$1 / E_0$"), TeX("$1 / E_1$"), TeX("$1 / E_2$")),
        pch = c(15, 16, 17),
        lty = c(1, 1, 1),
        col = c(4, 3, 2),
        pt.cex = 1,
        cex = 1,
        bty = "n"
      )
      dev.off()
      mypath <-
        paste0("../_manuscript/figures/error_exponents_plot_",
               a1, "_", a2, "_", bb, ".pdf")
      pdf(file = mypath,
          width = 5,
          height = 5)
      par(
        # bg = "white",
        # mgp = c(2, 0.8, 0),
        mar = c(3.2, 3.2, 1, 1),
        cex = 1.2
      )
      mus_plot_points <-
        seq(low_mu, up_mu, inc_mu * ratio_plot_points)
      mus_index_plot_points <- seq(1, len_mu, ratio_plot_points)
      log_mu <- log(as.matrix(mus))
      log_mu_plot_points <- log(as.matrix(mus_plot_points))
      plot(
        mus_plot_points,
        -log(max_error_vec[mus_index_plot_points, 1]) / log_mu_plot_points,
        log = "x",
        ylim = c(-0.6, 2.1),
        col = 4,
        pch = 15,
        xlab = "",
        ylab = "",
        cex.lab = 1.3,
        cex.axis = 1.3
      )
      points(
        mus_plot_points,
        log(max_error_vec[mus_index_plot_points, 2]) / -log_mu_plot_points,
        col = 3,
        pch = 16
      )
      points(
        mus_plot_points,
        log(max_error_vec[mus_index_plot_points, 3]) / -log_mu_plot_points,
        col = 2,
        pch = 17
      )
      lines(mus,
            -log(max_error_vec[, 1]) / log_mu,
            col = 4,
            lwd = 2)
      lines(mus,
            -log(max_error_vec[, 2]) / log_mu,
            col = 3,
            lwd = 2)
      lines(mus,
            -log(max_error_vec[, 3]) / log_mu,
            col = 2,
            lwd = 2)
      legend(
        "bottomright",
        legend = c(
          TeX("$-\\log(E_0) / \\log(\\mu)$"),
          TeX("$-\\log(E_1) / \\log(\\mu)$"),
          TeX("$-\\log(E_2) / \\log(\\mu)$")
        ),
        pch = c(15, 16, 17),
        lty = c(1, 1, 1),
        col = c(4, 3, 2),
        pt.cex = 1,
        cex = 1,
        bty = "n"
      )
      dev.off()
      end_time <- Sys.time()
      print(end_time - start_time)
    }
  }
}
