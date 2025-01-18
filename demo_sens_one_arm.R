rm(list = ls())
source("utils_sens_one_arm.R")
library(rootSolve)
library(dplyr)
library(ggplot2)
library(mice)
library(glmnet)
library(nleqslv)
library(SuperLearner)
library(flexmix)
# useful functions
expit <- function(x) {
  exp(x) / {
    1 + exp(x)
  }
}
seed <- 1
# parameters for data generation
delta_S <- 0
delta_R0 <- 0
delta_R1 <- 0
# sensitivity parameters for model estimation
gamma_S <- 0
gamma_R0 <- 0
gamma_R1 <- 0
B <- 50
# nuisance models are correctly specified
OR <- 1
PS <- 1
# sample sizes for single-arm trial and external controls
N_R <- 200
N_E <- 500
# generate simulated data
sim_list <- getdata_sens(
  seed = seed,
  N_R = N_R,
  N_E = N_E,
  delta_S = delta_S,
  delta_R1 = delta_R1,
  delta_R0 = delta_R0
)
sim_data <- sim_list$data
tau.truth <- sim_list$ate_truth

# two sensitivity parameters
## gamma_S: EC compatibility
## gamma_R: MNAR for single-arm trial and external controls
est.list <- SensDR(sim_data,
  gamma_S = gamma_S,
  gamma_R0 = gamma_R0,
  gamma_R1 = gamma_R1,
  OR = OR, PS = PS,
  # number of Gaussian mixtures for modeling conditional density
  k = 1
)

# bootstrap for variance estimation
est.boot.list <- lapply(1:B, function(seed) {
  set.seed(seed)
  N <- nrow(sim_data)
  idx <- sample(N, replace = TRUE)

  SensDR(sim_data[idx, ],
    gamma_S = gamma_S,
    gamma_R0 = gamma_R0,
    gamma_R1 = gamma_R1
  )$est$tau
})
est.list$est$var_boot <- do.call(rbind, est.boot.list) %>%
  apply(2, var)
est.list$est <- rbind(est.list$est,
  tau.truth = c(tau.truth, NA)
)
# estimates
est.list$est
# calibrated sensitivity parameters
est.list$gamma.t_star
