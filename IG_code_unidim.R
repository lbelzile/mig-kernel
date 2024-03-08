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

low_mu <- 10
up_mu <- 100000
inc_mu <- 100
ratio_plot_points <- 100

mus <- seq(low_mu, up_mu, inc_mu)
len_mu <- length(mus)

#############
# Functions #
#############

delta <- function (x, mu) {
  return ((x - mu) / sqrt(mu))
}

log_ratio_error <- function (x, mu, ord) {
  if (ord == 0) { # Order 0 correction
    return (log(dinvgauss(x, mean = mu, shape = mu ^ 2) / dnorm(delta(x, mu), 0, 1) * sqrt(mu)))
  } else if (ord == 1) { # Order 1 correction
    return (log(dinvgauss(x, mean = mu, shape = mu ^ 2) / dnorm(delta(x, mu), 0, 1) * sqrt(mu)) + (3 / 2 - delta(x, mu) ^ 2 / 2) * delta(x, mu) / sqrt(mu))
  } else { # ord = 2 : Order 2 correction
    return (log(dinvgauss(x, mean = mu, shape = mu ^ 2) / dnorm(delta(x, mu), 0, 1) * sqrt(mu)) + (3 / 2 - delta(x, mu) ^ 2 / 2) * delta(x, mu) / sqrt(mu) - (3 / 4 - delta(x, mu) ^ 2 / 2) * delta(x, mu) ^ 2 / mu)
  }
}

#############
# Main code #
#############

start_time <- Sys.time()

max_error_vec <- array(NaN, dim = c(len_mu, 3))

for (i in 1:len_mu) {
  musi <- mus[i]
  lst <- seq(musi - sqrt(musi), musi + sqrt(musi), 1)
  len_x <- length(lst)
  error_mat <- array(NaN, dim = c(len_mu, len_x, 3))
  for (j in 1:len_x) {
    x_j <- as.numeric(lst[j])
    error_mat[i, j, 1] <- abs(log_ratio_error(x_j, musi, 0))
    error_mat[i, j, 2] <- abs(log_ratio_error(x_j, musi, 1))
    error_mat[i, j, 3] <- abs(log_ratio_error(x_j, musi, 2))
  }
  max_error_vec[i, 1] <- max(error_mat[i, , 1])
  max_error_vec[i, 2] <- max(error_mat[i, , 2])
  max_error_vec[i, 3] <- max(error_mat[i, , 3])
}

mypath <- "../_manuscript/figures/errors_plot_unidim.pdf"
pdf(file = mypath, width = 9, height = 4)
par(#bg = "white", pty = "s",
    mgp = c(2, 0.8, 0), mar = c(3,3, 0.1, 0.5), cex = 1.2, mfrow = c(1,2))
mus_plot_points <- seq(low_mu, up_mu, inc_mu * ratio_plot_points)
mus_index_plot_points <- seq(1, len_mu, ratio_plot_points)
plot(mus_plot_points, 1 / max_error_vec[mus_index_plot_points, 1], log = "xy", ylim = c(1, 1 / min(max_error_vec[, 3])),
     col = 4, pch = 15, xlab = "", ylab = "", #cex.lab = 1.3, cex.axis = 1.3,
     xaxt = "n", yaxt = "n")
axis(1, at = c(10, 1e3, 1e5),
     labels = c("10",expression(10^3), expression(10^5)))
axis(2, at = c(1e0, 1e2, 1e4, 1e6, 1e8),
     labels = c(1, expression(10^2), expression(10^4),expression(10^6), expression(10^8)))
points(mus_plot_points, 1 / max_error_vec[mus_index_plot_points, 2], col = 3, pch = 16)
points(mus_plot_points, 1 / max_error_vec[mus_index_plot_points, 3], col = 2, pch = 17)
lines(mus, 1 / max_error_vec[, 1], col = 4, lwd = 2)
lines(mus, 1 / max_error_vec[, 2], col = 3, lwd = 2)
lines(mus, 1 / max_error_vec[, 3], col = 2, lwd = 2)
legend("topleft", legend = c(TeX("$1 / E_0$"), TeX("$1 / E_1$"), TeX("$1 / E_2$")),bty = "n",
       pch = c(15, 16, 17), lty = c(1, 1, 1), col = c(4, 3, 2), pt.cex = 1, cex = 1)
# dev.off()
# mypath <- "../_manuscript/figures/error_exponents_plot_unidim.pdf"
# par(bg = "white", pty = "s", mgp = c(2, 0.8, 0), mar = c(3.1, 0.8, 1.2, 0.1), cex = 1.3)
mus_plot_points <- seq(low_mu, up_mu, inc_mu * ratio_plot_points)
mus_index_plot_points <- seq(1, len_mu, ratio_plot_points)
log_mu <- log(as.matrix(mus))
log_mu_plot_points <- log(as.matrix(mus_plot_points))
plot(mus_plot_points, -log(max_error_vec[mus_index_plot_points, 1]) / log_mu_plot_points, log = "x", ylim = c(-0.6, 2.1),
     col = 4, pch = 15, xlab = "", ylab = "", #cex.lab = 1.3, cex.axis = 1.3,
     xaxt = "n")
axis(1, at = c(10, 1e3, 1e5),
     labels = c("10",expression(10^3), expression(10^5)))
points(mus_plot_points, -log(max_error_vec[mus_index_plot_points, 2]) / log_mu_plot_points, col = 3, pch = 16)
points(mus_plot_points, -log(max_error_vec[mus_index_plot_points, 3]) / log_mu_plot_points, col = 2, pch = 17)
lines(mus, -log(max_error_vec[, 1]) / log_mu, col = 4, lwd = 2)
lines(mus, -log(max_error_vec[, 2]) / log_mu, col = 3, lwd = 2)
lines(mus, -log(max_error_vec[, 3]) / log_mu, col = 2, lwd = 2)
legend("bottomright", legend = c(TeX("$-\\log(E_0) / \\log(\\mu)$"), TeX("$-\\log(E_1) / \\log(\\mu)$"), TeX("$-\\log(E_2) / \\log(\\mu)$")),
       pch = c(15, 16, 17), lty = c(1, 1, 1), col = c(4, 3, 2), pt.cex = 1, cex = 1, bty = "n")
dev.off()

end_time <- Sys.time()
print(end_time - start_time)
