setwd(this.path::here())
library(mig)
library(ggplot2)
library(patchwork)
set.seed(202503)

beta <- rep(1,2L)
xs <- seq(-2.5, 5, length.out = 101)
grid <- as.matrix(expand.grid(xs, xs))
grid <- grid[c(grid %*% beta) > 0,]



xdat <- mig::rtellipt(
  n = 100,
  beta = beta,
  mu = c(0,0),
  sigma = diag(c(1,2)),
  df = 5)
H_mig <- mig::kdens_bandwidth(
  x = xdat,
  beta = beta,
  family = "mig",
  type = "full",
  transformation = "none",
  method = "lcv")
kdens_mig <- mig_kdens(
  x = xdat,
  newdata = grid,
  Omega = H_mig,
  beta = beta,
  log = TRUE)

lcv_mig <- c(mig:::mig_loo(x = xdat, Omega = H_mig, beta = beta))
H_tnorm <- mig::kdens_bandwidth(
  x = xdat,
  beta = beta,
  family = "tnorm",
  type = "full",
  transformation = "none",
  method = "lcv")
lcv_tnorm <- c(mig:::tnorm_loo(x = xdat,
                               Omega = H_tnorm,
                               beta = beta))
kdens_tnorm <-tnorm_kdens(
  x = xdat,
  newdata = grid,
  Sigma = H_tnorm,
  beta = beta,
  log = TRUE)

g1 <- ggplot() +
  geom_polygon(data = data.frame(x = c(-2.5, -2.5, 2.5),
                                 y = c(-2.5, 2.5, -2.5)),
               mapping = aes(x = x, y = y),
               fill = "grey90",
               alpha = 0.5) +
  geom_abline(slope = -1, intercept = 0) +
  geom_contour(data = data.frame(x1 = grid[,1],
                                 x2 = grid[,2],
                                 z = kdens_mig),
               mapping = aes(x = x1, y = x2, z = z),
               bins = 9L,
               col = "grey") +
  geom_point(data = data.frame(x1 = xdat[,1],
                               x2 = xdat[,2],
                               weight = lcv_mig),
             mapping = aes(x = x1, y = x2, color = weight)) +
  scale_x_continuous(limits = c(-2.5, 5),
                     expand = expansion()) +
  scale_y_continuous(limits = c(-2.5, 5),
                     expand = expansion()) +
  scale_color_continuous(type = "viridis",
                         limits = range(c(lcv_mig, lcv_tnorm))) +
  labs(x = expression(x[1]), y = expression(x[2]),
       color = "likelihood cross validation score",
       subtitle = "MIG kernel") +
  theme_classic() +
  theme(legend.position = "bottom")

g2 <- ggplot() +
  geom_contour(data = data.frame(x1 = grid[,1],
                                 x2 = grid[,2],
                                 z = kdens_tnorm),
               mapping = aes(x = x1, y = x2, z = z),
               bins = 9L,
               col = "grey") +
  geom_polygon(data = data.frame(x = c(-2.5, -2.5, 2.5),
                                 y = c(-2.5, 2.5, -2.5)),
               mapping = aes(x = x, y = y),
               fill = "grey90",
               alpha = 0.5) +
  geom_abline(slope = -1, intercept = 0) +
  geom_point(data = data.frame(
    x1 = xdat[,1],
    x2 = xdat[,2],
    weight = lcv_tnorm),
       mapping = aes(x = x1, y = x2, color = weight)) +
  scale_x_continuous(limits = c(-2.5, 5),
                     expand = expansion()) +
  scale_y_continuous(limits = c(-2.5, 5),
                     expand = expansion()) +
  scale_color_continuous(type = "viridis",
                         limits = range(c(lcv_mig, lcv_tnorm))) +
  labs(x = expression(x[1]), y = expression(x[2]),
       color = "likelihood cross validation score",
       subtitle = "truncated Gaussian kernel") +
  theme_classic() +
  theme(legend.position = "bottom")

g1 + g2 +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom', 
        legend.margin = margin(0,0,0,0),
        plot.margin = unit(c(1,2,1,1), "mm"))
ggsave(filename = "../_manuscript/figures/LCV_kernel.pdf",
       width = 8,
       height = 4.3)
