setwd(this.path::here())
library(simsalapar)
library(dplyr)
library(ggplot2)
library(patchwork)


# Timing reported in the text
round(
  t(apply(
    results[-7, "timing", , "500", "3", ],
    1,
    quantile,
    probs = c(0.25, 0.5, 0.75),
    na.rm = TRUE
  )),
  1
)
round(
  t(apply(
    results[-7, "timing", , "1000", "3", ],
    1,
    quantile,
    probs = c(0.25, 0.5, 0.75),
    na.rm = TRUE
  )),
  1
)
# The big differences in timing are due to (a) the mig plug-in for mig (vs kernel for others) + (b) diagonal matrix for hsgauss
round(
  data.frame(
    "N.1000" = apply(
      results[, "timing", , "1000", "3", ],
      1,
      mean,
      na.rm = TRUE
    ),
    "N.500" = apply(results[, "timing", , "500", "3", ], 1, mean, na.rm = TRUE)
  ),
  1
)
# Maximum time
max(apply(results[, "timing", , "1000", "3", ], 1, max, na.rm = TRUE))
# All less than 1 minute for the bandwidth

set.seed(2024)
nrep <- tail(dim(results), 1)

# Relabel results to get consistent naming convention
dimnames(results)$bandwidth[c(4, 8)] <- c("amise-spher-mig", "normal-full-mig")
# Transform array to data frame for plotting
r_df <- array2DF(results, responseName = "coef") |>
  dplyr::filter(bandwidth != "normal-full-mig") |>
  dplyr::mutate(
    kernel = factor(
      stringr::str_split_i(bandwidth, pattern = "-", 3),
      labels = c("trans. Gauss", "MIG", "trunc. Gaussian")
    ),
    transfo = factor(stringr::str_split_i(bandwidth, pattern = "-", 2)),
    method = factor(stringr::str_split_i(bandwidth, pattern = "-", 1)),
    nobs = factor(nobs, levels = c("250", "500", "1000")),
    bandwidth = factor(
      bandwidth,
      levels = dimnames(results)$bandwidth[c(4, 2, 1, 3, 5:7, 9:11)],
      labels = LETTERS[1:10]
    ),
    dim = factor(dim),
    model = factor(model, labels = paste0("F", 1:4))
  ) |>
  dplyr::rename(dimension = dim)
# Boxplot of results
ggplot(
  data = r_df |> dplyr::filter(criterion == "rmise"),
  mapping = aes(x = bandwidth, y = coef, fill = kernel, color = nobs)
) +
  #ggdist::stat_slabinterval(p_limits = c(0, 0.95),) +
  geom_boxplot(outliers = FALSE, alpha = 0.5) +
  facet_grid(
    rows = vars(model),
    cols = vars(dimension),
    labeller = "label_both",
    scales = "free_y",
    axes = "all_x"
  ) +
  scale_color_grey() +
  MetBrewer::scale_fill_met_d(name = "Hiroshige") +
  labs(
    x = "bandwidth estimator",
    y = "",
    subtitle = "root mean integrated squared error",
    color = "sample size",
    fill = "kernel"
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      drop0trailing = TRUE,
      style_negative = "minus"
    ),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(filename = "../_manuscript/figures/rmise.pdf", width = 9, height = 15)

ggplot(
  data = r_df |> dplyr::filter(criterion == "kldiv"),
  mapping = aes(
    x = bandwidth,
    y = log1p(pmax(0, coef)),
    fill = kernel,
    color = nobs
  )
) +
  # ggdist::stat_slabinterval() +
  geom_boxplot(outliers = FALSE, alpha = 0.5) +
  facet_grid(
    rows = vars(model),
    cols = vars(dimension),
    labeller = "label_both",
    scales = "free_y",
    axes = "all_x"
  ) +
  scale_color_grey() +
  scale_y_continuous(
    # limits = c(0, NA),
    # transform = "log",
    breaks = scales::breaks_pretty(),
    labels = scales::label_number(
      drop0trailing = TRUE,
      style_negative = "minus"
    ),
    #expand = expansion(mult = c(0, 0.01))
  ) +
  MetBrewer::scale_fill_met_d(name = "Hiroshige") +
  labs(
    x = "bandwidth estimator",
    y = "",
    subtitle = "Kullback-Leibler divergence (log scale)",
    color = "sample size",
    fill = "kernel"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(filename = "../_manuscript/figures/kldiv.pdf", width = 9, height = 15)


# Boundary mise
results_bmise <- readRDS("simulation_results_R1_bmise.rds")
dimnames(results_bmise)$bandwidth[2] <- c("amise-spher-mig")
r_df_bmise <- array2DF(results_bmise, responseName = "coef") |>
  dplyr::mutate(
    kernel = factor(
      stringr::str_split_i(bandwidth, pattern = "-", 3),
      labels = c("trans. Gauss", "MIG", "trunc. Gaussian")
    ),
    transfo = factor(stringr::str_split_i(bandwidth, pattern = "-", 2)),
    method = factor(stringr::str_split_i(bandwidth, pattern = "-", 1)),
    nobs = factor(nobs, levels = c("250", "500", "1000")),
    bandwidth = factor(
      bandwidth,
      levels = dimnames(results_bmise)$bandwidth[c(2, 1, 3:10)],
      labels = LETTERS[1:10]
    ),
    dim = factor(dim),
    model = factor(model, labels = paste0("F", 1:4))
  ) |>
  dplyr::rename(dimension = dim)
ggplot(
  data = r_df_bmise |> dplyr::filter(criterion == "bmise", !is.na(coef)),
  mapping = aes(x = bandwidth, y = coef, fill = kernel, color = nobs)
) +
  #ggdist::stat_slabinterval(p_limits = c(0, 0.95),) +
  geom_boxplot(outliers = FALSE, alpha = 0.5) +
  facet_wrap(
    facets = vars(model, dimension),
    nrow = 4,
    ncol = 2,
    # rows = vars(model),
    # cols = vars(dimension),
    labeller = labeller(
      model = ~ paste("Model: ", .),
      dimension = ~ paste0("d=", .),
      .multi_line = FALSE
    ),
    scales = "free",
    axes = "all"
  ) +
  scale_color_grey() +
  MetBrewer::scale_fill_met_d(name = "Hiroshige") +
  labs(
    x = "bandwidth estimator",
    y = "",
    subtitle = "Boundary root mean integrated squared error",
    color = "sample size",
    fill = "kernel"
  ) +
  scale_y_continuous(
    labels = scales::label_number(drop0trailing = TRUE),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(filename = "../_manuscript/figures/rbmise.pdf", width = 9, height = 15)
