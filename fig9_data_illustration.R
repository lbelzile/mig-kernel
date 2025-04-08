setwd(this.path::here())
set.seed(2024)
library(ggplot2)
library(patchwork)
library(revdbayes)
# remotes::install_github("lbelzile/TruncatedNormal")
# remotes::install_github("lbelzile/mig")
library(mig)
nobs <- 1e3L
data(geomagnetic, package = "mig")
data <- geomagnetic
thresh <- 220
data <- data[data > thresh] - thresh
maxx <- max(data)
grimshaw_gp_mle(x = data)
# Simulate exact draws from the posterior distribution
post_samp <- revdbayes::rpost(
   n = 1e5L,
   model = "gp",
   data = data,
   thresh = 0,
   prior = set_prior(prior = "mdi", model = "gp")
)$sim_vals
colnames(post_samp) <- c("scale", "shape")
# Use 1000 first observations for smoothing
post <- post_samp[seq_len(nobs),]
# Compute maximum likelihood estimator and
beta <- c(1, maxx)
mle <- mig::fit_mig(x = post, beta = c(1, maxx), "mle")
xs <- seq(20, 160, by = 1)
ys <- seq(-0.5, 1, length.out = 101)
z <- expand.grid(x = xs, y = ys)
z$logd <- dmig(
   as.matrix(z),
   xi = mle$xi,
   Omega = mle$Omega,
   beta = beta,
   log = TRUE
)
# Create a bivariate scatter plot of posterior draws
# with contour curves for the MIG density
g1 <- ggplot() +
   geom_point(
      data = data.frame(scale = post[, 1], shape = post[, 2]),
      mapping = aes(x = scale, y = shape),
      col = "grey"
   ) +
   geom_contour(
      data = z,
      color = "black",
      alpha = 0.5,
      mapping = aes(
         x = x,
         y = y,
         z = 2 * (logd - max(logd))
      ),
      breaks = -qchisq(seq(0.05, 0.95, by = 0.10), df = 2)
   ) +
   geom_abline(intercept = 0, slope = -1 / maxx) +
   scale_y_continuous(limits = c(-0.5, 1),
                      breaks = seq(-0.5, 1, by = 0.5),
                      labels = c("-0.5","0","0.5","1"), expand = c(0, 0)) +
   scale_x_continuous(limits = c(20, 175), expand = c(0, 0)) +
   labs(x = expression(sigma[u]), y = expression(xi)) +
   theme_minimal() +
   theme(axis.title.y = element_text(
      angle = 0,
      vjust = 1,
      hjust = 1
   ),
   text=element_text(size=16), #change font size of all text
   axis.text=element_text(size=15), #change font size of axis text
   axis.title=element_text(size=18), #change font size of axis titles
   plot.title=element_text(size=15), #change font size of plot title
   legend.text=element_text(size=20), #change font size of legend text
   legend.title=element_text(size=20))

# Compute bandwidth for different distributions
bw1 <- kdens_bandwidth(
   x = post,
   beta = beta,
   family = "mig",
   type = "full",
   method = "lcv",
   buffer = 0
)
bw2 <- kdens_bandwidth(
   x = post,
   beta = beta,
   family = "mig",
   type = "full",
   method = "amise",
   approx = "mig",
   N = 1e4L)
# Check correlation between parameters in bandwidth
cov2cor(bw1)[1,2] - cov2cor(bw2)[1,2]
bw3 <- kdens_bandwidth(
   x = post,
   beta = beta,
   method = "lcv",
   family = "tnorm",
   type = "full",
   N = 1e4L
)
bw4 <- kdens_bandwidth(
   x = post,
   beta = beta,
   method = "lcv",
   family = "hsgauss",
   type = "diag",
   N = 1e4L
)

z$mixtdens1 <- mig_kdens(
   x = post,
   newdata = as.matrix(z)[, 1:2],
   Omega = bw1,
   beta = beta,
   log = TRUE
)
z$mixtdens2 <- mig_kdens(
   x = post,
   newdata = as.matrix(z)[, 1:2],
   Omega = bw2,
   beta = beta,
   log = TRUE
)
z$mixtdens3 <- tnorm_kdens(
   x = post,
   newdata = as.matrix(z)[, 1:2],
   Sigma = bw3,
   beta = beta,
   log = TRUE
)
z$mixtdens4 <- hsgauss_kdens(
   x = post,
   newdata = as.matrix(z)[, 1:2],
   Sigma = bw4,
   beta = beta,
   log = TRUE
)

g2 <- ggplot() +
   geom_hex(
      data = as.data.frame(post_samp),
      mapping = aes(x = scale, y = shape),
      bins = 100
   ) +
   scale_fill_gradient(low = "grey90", high = "grey10") +
   geom_contour(
      data = z,
      mapping = aes(x = x, y = y, z = mixtdens1),
      breaks = max(z$mixtdens1) - qchisq(seq(0.05, 0.95, by = 0.10), df = 2),
      col = "black"
   ) +
   geom_contour(
      data = z,
      mapping = aes(x = x, y = y, z = mixtdens3),
      breaks = max(z$mixtdens3) - qchisq(seq(0.05, 0.95, by = 0.10), df = 2),
      col = "grey50"
   ) +
   # geom_contour(
   #    data = z,
   #    mapping = aes(x = x, y = y, z = mixtdens3),
   #    breaks = -qchisq(seq(0.05, 0.95, by = 0.10), df = 2),
   #    linetype = "dashed",
   #    col = "black") +
   # geom_contour(
   #    data = z,
   #    mapping = aes(x = x, y = y, z = mixtdens4),
   #    breaks = -qchisq(seq(0.05, 0.95, by = 0.10), df = 2),
   #    linetype = "dashed",
#    col = "grey50") +
geom_abline(intercept = 0, slope = -1 / maxx) +
   labs(x = expression(sigma[u]), y = expression(xi)) +
   scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, by = 0.5),
                      labels = c("-0.5","0","0.5","1"), expand = c(0, 0)) +
   scale_x_continuous(limits = c(20, 175), expand = c(0, 0)) +
   theme_minimal() +
   theme(axis.title.y = element_text(
      angle = 0,
      vjust = 1,
      hjust = 1
   ),
   legend.position = "none",
   text=element_text(size=16), #change font size of all text
   axis.text=element_text(size=15), #change font size of axis text
   axis.title=element_text(size=18), #change font size of axis titles
   plot.title=element_text(size=15), #change font size of plot title
   legend.text=element_text(size=20), #change font size of legend text
   legend.title=element_text(size=20))

g1 + g2

ggsave("../_manuscript/figures/application.pdf",
       width = 10,
       height = 5)

g1 / g2
ggsave("../_manuscript/figures/application_tp.pdf",
       width = 5,
       height = 10)
