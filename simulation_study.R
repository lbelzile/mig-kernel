setwd(this.path::here())
source("utils.R")
library(mig)
library(simsalapar)

nrep <- 1e3L
nobs <- c(250L, 500L, 1000L)
modelid <- 1:4
dim <- 2:4

# Nsim <- 1000L
ncores <- 25L
block <- 40L
Nsim <- ncores*block
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = Nsim),
  n = list(type = "grid",  value = nobs),
  model = list(type = "grid",  value = modelid),
  d = list(type = "grid", value = dim),
  method = list(type = "inner", value = 1:8),
  criterion = list(type = "inner", value = c("RMISE","KLD"))
)
doOne <- function(n, model, d, method, criterion){
  beta <- rep(1, d)
  samp <- rsim(n = n, beta = beta, model = model)
  bw1 <- try(mig_kdens_bandwidth(
    x = samp,
    beta = beta,
    type = "full",
    method = "lcv"), silent = TRUE)
  bw2 <- try(mig_kdens_bandwidth(
    x = samp,
    beta = beta,
    type = "full",
    method = "amise",
    approx = "tnorm",
    buffer = 0), silent = TRUE)
  bw3 <- try(mig_kdens_bandwidth(
     x = samp,
     beta = beta,
     type = "full",
     method = "amise",
     approx = "mig",
     buffer = 0), silent = TRUE)
  bw4 <- try(mig_kdens_bandwidth(
    x = samp,
    beta = beta,
    method = "lcv",
    transformation = "spherical",
    type = "isotropic"), silent = TRUE)
  bw5 <- try(mig_kdens_bandwidth(
     x = samp,
     beta = beta,
     method = "amise",
     transformation = "spherical",
     approx = "tnorm",
     type = "isotropic"), silent = TRUE)
  bw6 <- try(mig_kdens_bandwidth(
     x = samp,
     beta = beta,
     method = "amise",
     transformation = "spherical",
     approx = "mig",
     type = "isotropic"), silent = TRUE)
  bwnorm <- diag(apply(samp, 2, sd)) * exp((log(4) - log(d + 2) - log(n))/(d+4))
  criteria <- matrix(NA, nrow = 8, ncol = 2)
  if(is.matrix(bw1) & !isTRUE(any(is.na(bw1)))){
   criteria[1,] <- metrics(x = samp, model = model, bandwidth = bw1, beta = beta)
  }
  if(is.matrix(bw2) & !isTRUE(any(is.na(bw2)))){
    criteria[2,] <- metrics(x = samp, model = model, bandwidth = bw2, beta = beta)
  }
  if(is.matrix(bw3) & !isTRUE(any(is.na(bw3)))){
    criteria[3,] <- metrics(x = samp, model = model, bandwidth = bw3, beta = beta)
  }
  if(is.matrix(bw4) & !isTRUE(any(is.na(bw4)))){
    criteria[4,] <- metrics(x = samp, model = model, bandwidth = bw4, beta = beta)
  }
  if(is.matrix(bw5) & !isTRUE(any(is.na(bw5)))){
     criteria[5,] <- metrics(x = samp, model = model, bandwidth = bw5, beta = beta)
  }
  if(is.matrix(bw6) & !isTRUE(any(is.na(bw6)))){
     criteria[6,] <- metrics(x = samp, model = model, bandwidth = bw6, beta = beta)
  }
  criteria[7,] <- metrics(x = samp, model = model, bandwidth = bwnorm,
                          beta = beta, kernel = "tnorm")
  criteria[8,] <- metrics(x = samp, model = model, bandwidth = bwnorm, beta = beta)
  criteria
}
