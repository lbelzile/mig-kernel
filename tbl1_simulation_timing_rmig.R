library(mig)
library(bench)
# Dimension of the vector
dims <- c(2:5, 10L)
set.seed(2024)
results <- press(
  d = c(2L,3L,5L,10L),
  {
  beta <- exp(rnorm(n = d, mean = 0, sd = 1))
  xi <- beta + rexp(n = d, rate = 0.5)
  Omega <- rWishart(n = 1,
                  df = d + 2L,
                  Sigma = matrix(0.5, d, d) + diag(rep(0.5, d)))[,,1]
  bench::mark(
    invsim = rmig(n = 100L, xi = xi, beta = beta, Omega = Omega, method = "invsim"),
    bm = rmig(n = 100L, xi = xi, beta = beta, Omega = Omega, method = "bm"),
    check = FALSE,
    iterations = 25L,
    time_unit = "ms")
  }
)
save(results, file = "timing_rmig.RData", version = 2)

matres <- matrix(results$median, nrow = 2)
colnames(matres) <- c(2L,3L,5L,10L)
rownames(matres) <- c("exact","Minami")
knitr::kable(matres, digits = 2, format = "latex", booktabs = TRUE)
