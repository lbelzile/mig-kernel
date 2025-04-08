setwd(this.path::here())
# Case 1: all entries of beta are positive
library(mig)
source("utils.R")

pdf(file = "../_manuscript/figures/simulation_2d_densities.pdf", width = 8, height = 8)
par(mfrow = c(2,2), mar = c(4.5,4.5,1,1))
xg <- seq(-5, 20, length.out = 101)
grid <- as.matrix(expand.grid(xg,xg))
beta <- rep(1, 2)
d1 <- dsim(x = grid, model = 1, beta = beta, log = TRUE)
d2 <- dsim(x = grid, model = 2, beta = beta, log = TRUE)
image(xg, xg, z = matrix(d1, nrow = length(xg)),
      breaks = max(d1, na.rm = TRUE) + sort(seq(0, -20, by = -2)),
      col = viridis::viridis(n = 10, alpha = 0.8),
      xlab = expression(X[1]),
      ylab = expression(X[2]))
contour(xg, xg,
        z = matrix(d1, nrow = length(xg)),
        levels = max(d1, na.rm = TRUE) + sort(seq(0, -20, by = -2)),
        labels = seq(-20, 0, by = 2),
        add = TRUE)
text(min(xg), min(xg), labels = expression(F[1]),
     cex = 2.5, adj = c(-0.2,-0.2))
abline(a = 0, b = -1)

image(xg, xg,
      z = matrix(d2, nrow = length(xg)),
      breaks = max(d2, na.rm = TRUE) + sort(seq(0, -20, by = -2)),
      col = viridis::viridis(n = 10, alpha = 0.8),
      xlab = expression(X[1]),
      ylab = expression(X[2]))
contour(xg, xg,
        z = matrix(d2, nrow = length(xg)),
        levels = max(d2, na.rm = TRUE) + sort(seq(0, -20, by = -2)),
        labels = seq(-20, 0, by = 2),
        add = TRUE)
text(min(xg), min(xg), labels = expression(F[2]),
     cex = 2.5, adj = c(-0.2,-0.2))
abline(a = 0, b = -1)


xg <- seq(-1, 6, length.out = 101)
grid <- as.matrix(expand.grid(xg,xg))
d3 <- dsim(x = grid, model = 3, beta = beta, log = TRUE)
image(xg, xg,
      z = matrix(d3, nrow = length(xg)),
      breaks = max(d3, na.rm = TRUE) + c(-2000, sort(seq(0, -20, by = -2))),
      col = viridis::viridis(n = 11, alpha = 0.8),
      xlab = expression(X[1]),
      ylab = expression(X[2]))
contour(xg, xg,  z = matrix(d3, nrow = length(xg)),
        levels = max(d3, na.rm = TRUE) + sort(seq(0, -20, by = -2)),
        labels = seq(-20, 0, by = 2),
        add = TRUE)
text(min(xg), min(xg), labels = expression(F[3]),
     cex = 2.5, adj = c(-0.2,-0.2))
abline(a = 0, b = -1)

xg <- seq(-5, 20, length.out = 101)
grid <- as.matrix(expand.grid(xg, xg))
d4 <- dsim(x = grid, model = 4, beta = beta, log = TRUE)
d4 <- ifelse(grid %*% beta < 0, NA, d4)
image(xg, xg,
      z = matrix(d4, nrow = length(xg)),
      breaks = max(d4, na.rm = TRUE) + sort(c(-1000, 0,-0.25,-0.5,-1,seq(-2, -20, by = -2))),
      col = viridis::viridis(n = 14, alpha = 0.8),
      xlab = expression(X[1]),
      ylab = expression(X[2]))
contour(xg, xg,
        z = matrix(d4, nrow = length(xg)),
        levels = max(d4, na.rm = TRUE) + rev(c(0,-0.5,-1,seq(-2, -20, by = -2))),
        labels = c(seq(-20, -2, by = 2), -1,-0.5,0),
        add = TRUE)
text(min(xg), min(xg), labels = expression(F[4]),
     cex = 2.5, adj = c(-0.2,-0.2))
abline(a = 0, b = -1)
dev.off()
