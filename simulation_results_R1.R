setwd(paste0(this.path::here(), "/results"))
files <- list.files(recursive = FALSE, pattern = "*.rds")

results <- array(dim = c(12, 5, 4, 3, 2, 1000))
dimnames(results) <- list(
  bandwidth = c(
    "lcv-full-mig",
    "amise-full-mig",
    "rlcv-full-mig-bad",
    "rlcv-full-mig",
    "amise-isotropic-spherical",
    "lcv-full-tn",
    "rlcv-full-tn",
    "lscv-full-tn",
    "normal-mig",
    "lcv-spher-hsgauss",
    "rlcv-spher-hsgauss",
    "lscv-spher-hsgauss"
  ),
  criterion = c("rwmise", "rmise", "rbmise", "kldiv", "timing"),
  model = 1:4,
  nobs = c(250, 500, 1000),
  dim = c(2, 3),
  rep = 1:1000
)

# n d
# 1 250 2
# 2 500 2
# 3 250 3
# 4 500 3
library(dplyr)
str <- stringr::str_split(files, '_', simplify = TRUE)
batch <- case_when(
  str[, 4] == "HS" ~ 4,
  str[, 4] == "TN" ~ 3,
  str[, 4] == "MIG2" ~ 2,
  .default = 1
)
mod <- as.integer(ifelse(nchar(str[, 4]) == 1L, str[, 4], str[, 5]))
b <- as.integer(stringr::str_split(
  ifelse(nchar(str[, 4]) == 1L, str[, 5], str[, 6]),
  ".r",
  simplify = TRUE
)[, 1])
dd <- dim(results)[4:5]
ncombo <- prod(dd)
se_nobs <- rep(seq_len(dd[1]), length.out = ncombo)
se_dim <- rep(seq_len(dd[2]), each = dd[1])
pos <- case_when(
  batch == 1 ~ list(1:3),
  batch == 2 ~ list(4:5),
  batch == 3 ~ list(6:9),
  batch == 4 ~ list(10:12)
)
for (i in seq_along(files)) {
  test <- readRDS(file = files[i])
  if (length(test) == ncombo) {
    for (j in 1:ncombo) {
      results[pos[[i]], , mod[i], se_nobs[j], se_dim[j], b[i]] <- test[[
        j
      ]]$value
    }
  }
}
results <- results[-3, , , , , ]
setwd(this.path::here())
saveRDS(results, file = "simulation_results_R1.rds", version = 2)
