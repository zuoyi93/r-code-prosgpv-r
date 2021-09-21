# -----------------------
# GET SIMULATION RESULTS
# -----------------------

library("ProSGPV")
library("BeSS")
library("SIS")
library("survival")
library("pROC")
library("doParallel")

source("utils.R")

# -----------------------
# FIGURES 4
# -----------------------

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

# gaussian
out.main.gaussian <- foreach(
  p = seq(100, 1000, 50),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    p = p, n = 100, s = 10, beta.min = 1, beta.max = 2,
    intercept = 0, family = "gaussian"
  )
}

# logistic
out.main.logistic <- foreach(
  n = seq(32, 320, 16),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 16, s = 6, beta.min = 0.4, beta.max = 1.2, rho = 0.6,
    intercept = 0, family = "binomial", gvif = T
  )
}

# poisson
out.main.poisson <- foreach(
  p = seq(40, 400, 20),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = 40, p = p, s = 4, beta.min = 0.2, beta.max = 0.5,
    intercept = 2, family = "poisson"
  )
}

# cox
out.main.cox <- foreach(
  n = seq(80, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 40, s = 20, beta.min = 0.3, beta.max = 1,
    intercept = 0, family = "cox"
  )
}

# run time
# gaussian
out.run.gaussian <- foreach(
  p = seq(100, 1000, 50),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  run.time.many(
    p = p, n = 100, s = 10, beta.min = 1, beta.max = 2,
    intercept = 0, family = "gaussian"
  )
}

# logistic
out.run.logistic <- foreach(
  n = seq(32, 320, 16),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  run.time.many(
    n = n, p = 16, s = 6, beta.min = 0.4, beta.max = 1.2, rho = 0.6,
    intercept = 0, family = "binomial", gvif = T
  )
}

# poisson
out.run.poisson <- foreach(
  p = seq(40, 400, 20),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  run.time.many(
    n = 40, p = p, s = 4, beta.min = 0.2, beta.max = 0.5,
    intercept = 2, family = "poisson"
  )
}

# cox
out.run.cox <- foreach(
  n = seq(80, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  run.time.many(
    n = n, p = 40, s = 20, beta.min = 0.3, beta.max = 1,
    intercept = 0, family = "cox"
  )
}

save(out.main.gaussian, out.main.logistic, out.main.poisson, out.main.cox,
     out.run.gaussian, out.run.logistic, out.run.poisson, out.run.cox,
     file = "data/sim.data.RData"
)
