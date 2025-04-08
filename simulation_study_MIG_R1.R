setwd(this.path::here())
source("utils.R")
library(mig)
library(simsalapar)
nobs <- c(250L, 500L, 1000L)
modelid <- 1:4
dim <- 2:3

ncombo <- length(nobs) * length(modelid) * length(dim)

# Nsim <- 1000L
ncores <- 6L
block <- 1L #5L
Nsim <- 1
N <- 1000
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = Nsim),
  n = list(type = "grid", value = nobs),
  model = list(type = "grid", value = modelid),
  d = list(type = "grid", value = dim),
  method = list(
    type = "inner",
    value = c(2, 6, 1, 3, 12:14, 9:11)
  ),
  criterion = list(
    type = "inner",
    value = c("rmise", "kldiv", "rwmise", "rbmise", "timing")
  )
)
# Move to local directory
if (!dir.exists("results")) {
  dir.create("results")
}
setwd("results")
name <- "res-mig_mod_"
for (i in seq(1, N, by = Nsim)) {
  if (!paste0(name, "_", i, ".rds") %in% list.files()) {
    result <- doLapply(
      vList = varList,
      doAL = FALSE,
      sfile = paste0(name, "_", i, ".rds"),
      keepSeed = FALSE,
      seed = i + 1000,
      doOne = doOne
    )
  }
}
setwd(this.path::here())
