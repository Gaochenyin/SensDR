getdata_sens <- function(seed, N_R, N_E,
                         delta_S = 0,
                         delta_R1 = 0,
                         delta_R0 = 0) {
  set.seed(seed)
  # N_R <- N_E <- 5e5
  N <- N_R + N_E
  # original covariates
  X <- cbind(
    replicate(4, rnorm(N, mean = 0.25)),
    rnorm(N)
  )

  d.X <- ncol(X)
  Z <- {
    X**2 + 2 * sin(X) - 1.5
  } / sqrt(2)
  Z[, 5] <- X[, 5]

  Y1 <- rnorm(N, mean = 3 * Z %*% rep(1 / 6, 5))
  Y0 <- rnorm(N, mean = 2 * Z %*% rep(1 / 6, 5))

  # generate data origin
  alpha.opt <- nleqslv::nleqslv(c(0), fn = function(alpha) {
    mean(expit(alpha + Z %*% rep(0.1, 5) + delta_S * Y0) * N) - N_R
  })$x
  S <- rbinom(N,
    size = 1,
    prob = expit(alpha.opt + Z %*% rep(0.1, 5) + delta_S * Y0)
  )

  ## S = 1
  # transformed covariates
  X_R <- X[S == 1, ]
  Z_R <- Z[S == 1, ]
  Y1_R <- Y1[S == 1]
  Y0_R <- Y0[S == 1]
  # generate treatment (all treated)
  A_R <- 1 # rbinom(sum(S), size = 1, prob = expit(Z_R[,-5]%*%rep(0.1, 4)))
  # generate missing
  alphaR1.opt <- nleqslv::nleqslv(c(0), fn = function(alpha) {
    mean(expit(alpha + (2 * A_R - 1) * Z_R %*% rep(1 / 6, 5) + delta_R1 * Y1_R)) - 1 / 2
  })$x

  prob_R1.X <- expit(alphaR1.opt +
    (2 * A_R - 1) * Z_R %*% rep(1 / 6, 5) + delta_R1 * Y1_R)
  R_R <- 1 - rbinom(sum(S),
    size = 1,
    prob = prob_R1.X
  )
  # generate outcomes
  Y_R <- Y1_R * R_R + 0 * (1 - R_R)


  ## S = 0
  # transformed covariates
  X_E <- X[S == 0, ]
  Z_E <- Z[S == 0, ]
  Y0_E <- Y0[S == 0]
  # generate missingness for external resources
  alphaR0.opt <- nleqslv::nleqslv(c(0), fn = function(alpha) {
    mean(expit(alpha + Z_E %*% rep(1 / 6, 5) + delta_R0 * Y0_E)) - 1 / 2
  })$x

  prob_E0.X <- expit(alphaR0.opt + Z_E %*% rep(1 / 6, 5) + delta_R0 * Y0_E)
  prob_R0.X <- expit(alphaR0.opt + Z_R %*% rep(1 / 6, 5) + delta_R0 * Y0_R)

  R_E <- 1 - rbinom(sum(S == 0),
    size = 1,
    prob = prob_E0.X
  )
  # generate outcomes (same at first)
  Y_E <- Y0_E * R_E +
    0 * (1 - R_E)
  # use correct model
  sim_data <- data.frame(
    S = c(
      rep(1, nrow(X_R)),
      rep(0, nrow(X_E))
    ),
    X = rbind(X_R, X_E),
    Z = rbind(Z_R, Z_E),
    R = c(R_R, R_E),
    Y = c(Y_R, Y_E)
  )

  # Y1: S = 1, R = 1
  mu1.R1X <- predict(
    lm(Y ~ .,
      data = data.frame(
        Y = Y1_R[R_R == 1],
        Z = Z_R[R_R == 1, ]
      )
    ),
    newdata = data.frame(Z = Z_R)
  )

  # Y1: S = 1, R = 0
  mu1.R0X <- predict(
    lm(Y ~ .,
      data = data.frame(
        Y = Y1_R[R_R == 0],
        Z = Z_R[R_R == 0, ]
      )
    ),
    newdata = data.frame(Z = Z_R)
  )

  # Y0: S = 1, R = 1
  mu0.R1X <- predict(
    lm(Y ~ .,
      data = data.frame(
        Y = Y0_R[R_R == 1],
        Z = Z_R[R_R == 1, ]
      )
    ),
    newdata = data.frame(Z = Z_R)
  )
  # Y0: S = 1, R= 0
  mu0.R0X <- predict(
    lm(Y ~ .,
      data = data.frame(
        Y = Y0_R[R_R == 0],
        Z = Z_R[R_R == 0, ]
      )
    ),
    newdata = data.frame(Z = Z_R)
  )

  ate_truth <- mean(prob_R1.X * mu1.R1X + (1 - prob_R1.X) * mu1.R0X -
    prob_R0.X * mu0.R1X - (1 - prob_R0.X) * mu0.R0X)

  return(list(
    data = sim_data,
    ate_truth = ate_truth
  ))
}


SensDR <- function(sim_data,
                   gamma_S = 0,
                   gamma_R0 = 0,
                   gamma_R1 = 0,
                   OR = TRUE, PS = TRUE,
                   predictors.OR = NULL,
                   predictors.PS = NULL,
                   k = 1 # number of Gaussian mixtures
) {
  N_R <- sum(sim_data$S)
  N_E <- sum(1 - sim_data$S)
  N <- N_R + N_E

  if (is.null(predictors.OR)) {
    if (OR) {
      predictors.OR <- paste0("Z.", 1:5)
    } else {
      predictors.OR <- paste0("X.", 1:5)
    }
  }

  if (is.null(predictors.PS)) {
    if (PS) {
      predictors.PS <- paste0("Z.", 1:5)
    } else {
      predictors.PS <- paste0("X.", 1:5)
    }
  }

  # construct the estimating function
  ## pi_S(X), pi_R(X), tau
  # estimate nuisance functions by GMM
  m1.model <- flexmix::flexmix(as.formula(paste0("Y~", paste0(predictors.OR, collapse = "+"))),
    data = sim_data[sim_data$S == 1 & sim_data$R == 1, ],
    k = k,
    model = list(FLXMRglm(family = "gaussian")),
    control = list(iter = 100)
  )

  m0.model <- flexmix::flexmix(as.formula(paste0("Y~", paste0(predictors.OR, collapse = "+"))),
    data = sim_data[sim_data$S == 0 & sim_data$R == 1, ],
    k = k,
    model = list(FLXMRglm(family = "gaussian")),
    control = list(iter = 100)
  )

  fit.S1 <- glm(paste0("S~", paste(predictors.PS, collapse = "+")),
    family = "binomial",
    data = sim_data
  )

  fit.R1 <- glm(paste0("R~", paste(predictors.PS, collapse = "+")),
    family = "binomial",
    data = sim_data,
    subset = S == 1
  )

  fit.R0 <- glm(paste0("R~", paste(predictors.PS, collapse = "+")),
    family = "binomial",
    data = sim_data,
    subset = S == 0
  )

  S1.X <- predict(fit.S1,
    newdata = sim_data[, predictors.PS],
    type = "response"
  )
  R1.X <- predict(fit.R1,
    newdata = sim_data[, predictors.PS],
    type = "response"
  )
  R0.X <- predict(fit.R0,
    newdata = sim_data[, predictors.PS],
    type = "response"
  )


  # predict the outcomes under treatment 0 and 1
  mu0.X.all <- predict(m0.model, newdata = sim_data)
  mu1.X.all <- predict(m1.model, newdata = sim_data)

  pi.m1 <- flexmix::prior(m1.model)
  pi.m0 <- flexmix::prior(m0.model)

  sigma.m1 <- parameters(m1.model)["sigma", ]
  sigma.m0 <- parameters(m0.model)["sigma", ]

  sim_data$mu1.X <- apply(
    mapply(
      function(mu, pi) {
        mu * pi
      },
      mu = mu1.X.all, pi = pi.m1
    ),
    1, mean
  )

  sim_data$mu0.X <- apply(
    mapply(
      function(mu, pi) {
        mu * pi
      },
      mu = mu0.X.all, pi = pi.m0
    ),
    1, mean
  )
  q_R1.X <- (1 - R1.X) / R1.X
  q_R0.X <- (1 - R0.X) / R0.X
  q_S.X <- S1.X / (1 - S1.X)

  ## sensitivity analysis, s(y) = y * phi(y)
  # gamma_0 <- 0
  # compute the distribution by Gaussian distribution
  den_R0_S <- mapply(function(mu, pi, sigma) {
    pi * exp(mu * (gamma_R0 + gamma_S) + sigma**2 * (gamma_R0 + gamma_S)**2 / 2)
  }, pi = pi.m0, mu = mu0.X.all, sigma = sigma.m0) %>% apply(1, sum)

  den_R0 <- mapply(function(mu, pi, sigma) {
    pi * exp(mu * (gamma_R0) + sigma**2 * (gamma_R0)**2 / 2)
  }, pi = pi.m0, mu = mu0.X.all, sigma = sigma.m0) %>% apply(1, sum)


  den_S <- mapply(function(mu, pi, sigma) {
    pi * exp(mu * (gamma_S) + sigma**2 * (gamma_S)**2 / 2)
  }, pi = pi.m0, mu = mu0.X.all, sigma = sigma.m0) %>% apply(1, sum)

  den_R1 <- mapply(function(mu, pi, sigma) {
    pi * exp(mu * (gamma_R1) + sigma**2 * (gamma_R1)**2 / 2)
  }, pi = pi.m1, mu = mu1.X.all, sigma = sigma.m1) %>% apply(1, sum)

  ## compute the expectation for the numerator


  num_R0_S <- mapply(function(mu, pi, sigma) {
    pi * (mu + (gamma_R0 + gamma_S) * sigma**2) *
      exp(mu * (gamma_R0 + gamma_S) + sigma**2 * (gamma_R0 + gamma_S)**2 / 2)
  }, pi = pi.m0, mu = mu0.X.all, sigma = sigma.m0) %>% apply(1, sum)


  num_R0 <- mapply(function(mu, pi, sigma) {
    pi * (mu + (gamma_R0) * sigma**2) *
      exp(mu * (gamma_R0) + sigma**2 * (gamma_R0)**2 / 2)
  }, pi = pi.m0, mu = mu0.X.all, sigma = sigma.m0) %>% apply(1, sum)

  num_S <- mapply(function(mu, pi, sigma) {
    pi * (mu + (gamma_S) * sigma**2) *
      exp(mu * (gamma_S) + sigma**2 * (gamma_S)**2 / 2)
  }, pi = pi.m0, mu = mu0.X.all, sigma = sigma.m0) %>% apply(1, sum)

  num_R1 <- mapply(function(mu, pi, sigma) {
    pi * (mu + (gamma_R1) * sigma**2) *
      exp(mu * (gamma_R1) + sigma**2 * (gamma_R1)**2 / 2)
  }, pi = pi.m1, mu = mu1.X.all, sigma = sigma.m1) %>% apply(1, sum)

  # compute the adjusted outcomes
  ## Y0
  sim_data$Y0.adj <- sim_data$Y *
    exp(gamma_R0 * sim_data$Y) / den_R0
  sim_data$Y0.adj.E <- num_R0 / den_R0
  sim_data$Y0.adj2.E <- num_R0 / den_R0^2 * exp(gamma_R0 * sim_data$Y)

  ## Y1
  sim_data$Y1.adj <- sim_data$Y *
    exp(gamma_R1 * sim_data$Y) / den_R1
  sim_data$Y1.adj.E <- num_R1 / den_R1
  sim_data$Y1.adj2.E <- num_R1 / den_R1^2 * exp(gamma_R1 * sim_data$Y)


  # other necessary terms
  sim_data$B_S <- num_S
  sim_data$A_S <- den_S
  sim_data$B_R0 <- num_R0
  sim_data$A_R0 <- den_R0
  sim_data$B_R0_S <- num_R0_S
  sim_data$A_R0_S <- den_R0_S

  sim_data$num_mu0.t2t <- with(
    sim_data,
    R0.X * B_S * A_R0 + (1 - R0.X) * B_R0_S
  )
  sim_data$den_mu0.t2t <- with(
    sim_data,
    R0.X * A_S * A_R0 + (1 - R0.X) * A_R0_S
  )

  # naive method: weighting
  tau.OM <- (with(
    sim_data,
    sum(
      S * Y1.adj.E -
        S * num_mu0.t2t / den_mu0.t2t,
      na.rm = TRUE
    )
  )) / N_R


  tau.PS <- (with(
    sim_data,
    sum(
      S * R * Y / R1.X -
        (1 - S) * R / R0.X * q_S.X * Y,
      na.rm = TRUE
    )
  )) / N_R


  ## Jump2Reference with two-parameter framework
  ## by EIF
  tau.hat.J2R <-
    sum(
      with(
        sim_data,
        (
          S * R1.X * mu1.X +
            S * (R * Y - R1.X * mu1.X) -
            S * R *
              {
                num_mu0.t2t
              } / {
                den_mu0.t2t
              } -

            (1 - S) * q_S.X * R1.X *
              {
                (R - R0.X) * (B_S * A_R0 - B_R0_S) +
                  R / R0.X * (R0.X * A_R0 * (Y * exp(gamma_S * Y) - B_S) +
                    R0.X * B_S * exp(gamma_R0 * Y) +
                    (1 - R0.X) * Y * exp((gamma_R0 + gamma_S) * Y))
              } / den_mu0.t2t +

            (1 - S) * q_S.X * R1.X *
              {
                (R - R0.X) * (A_S * A_R0 - A_R0_S) +
                  R / R0.X * (R0.X * A_S * exp(gamma_R0 * Y) +
                    R0.X * A_R0 * (exp(gamma_S * Y) - A_S) +
                    (1 - R0.X) * exp((gamma_R0 + gamma_S) * Y))
              } *
              num_mu0.t2t / den_mu0.t2t^2
        )
      ),
      na.rm = TRUE
    ) / N_R
  ## three-parameter tilting framework
  tau.hat.t2t <- sum(
    with(
      sim_data,
      (S * R1.X * mu1.X + S * {
        R * Y - R1.X * mu1.X
      } +
        S * (1 - R) * Y1.adj.E +
        S * R * q_R1.X * {
          Y1.adj - Y1.adj2.E
        } -
        S *
          {
            num_mu0.t2t
          } / {
            den_mu0.t2t
          } -
        (1 - S) * q_S.X *
          {
            (R - R0.X) * (B_S * A_R0 - B_R0_S) +
              R / R0.X * (R0.X * A_R0 * (Y * exp(gamma_S * Y) - B_S) +
                R0.X * B_S * exp(gamma_R0 * Y) +
                (1 - R0.X) * Y * exp((gamma_R0 + gamma_S) * Y))
          } / den_mu0.t2t +
        (1 - S) * q_S.X *
          {
            (R - R0.X) * (A_S * A_R0 - A_R0_S) +
              R / R0.X * (R0.X * A_S * exp(gamma_R0 * Y) +
                R0.X * A_R0 * (exp(gamma_S * Y) - A_S) +
                (1 - R0.X) * exp((gamma_R0 + gamma_S) * Y))
          } *
          num_mu0.t2t / den_mu0.t2t^2)
    ),
    na.rm = TRUE
  ) / N_R

  # calibration for tilting function
  ## for sensitivity parameter gamma_R0 (R|S=0)
  fit.R0.all <- glm(as.formula(paste0("(1-R)~", paste0(predictors.PS, collapse = "+"))),
    data = sim_data, subset = S == 0,
    family = "binomial"
  )
  var.m.R0 <- var(predict(fit.R0.all))
  var.m.R0_j <- sapply(1:length(predictors.PS), function(idx) {
    fit.R0 <- glm(as.formula(paste0("(1-R)~", paste0(predictors.PS[-idx], collapse = "+"))),
      data = sim_data, subset = S == 0,
      family = "binomial"
    )
    var(predict(fit.R0))
  })
  rho.R0.star <- max((var.m.R0 - var.m.R0_j) / (var.m.R0 + pi^2 / 3))

  ## for sensitivity parameter gamma_R1 (R|S=1)
  fit.R1.all <- glm(as.formula(paste0("(1-R)~", paste0(predictors.PS, collapse = "+"))),
    data = sim_data, family = "binomial", subset = S == 1
  )
  var.m.R1 <- var(predict(fit.R1.all))

  var.m.R1_j <- sapply(1:length(predictors.PS), function(idx) {
    fit.R1 <- glm(as.formula(paste0("(1-R)~", paste0(predictors.PS[-idx], collapse = "+"))),
      data = sim_data, family = "binomial", subset = S == 1
    )
    var(predict(fit.R1))
  })

  rho.R1.star <- max((var.m.R1 - var.m.R1_j) / (var.m.R1 + pi^2 / 3))

  ## for sensitivity parameter gamma_S
  fit.S.all <- glm(as.formula(paste0("S~", paste0(predictors.PS, collapse = "+"))),
    data = sim_data, family = "binomial"
  )
  var.m.S <- var(predict(fit.S.all))

  var.m.S_j <- sapply(1:length(predictors.PS), function(idx) {
    fit.S <- glm(as.formula(paste0("S~", paste0(predictors.PS[-idx], collapse = "+"))),
      data = sim_data, family = "binomial"
    )
    var(predict(fit.S))
  })

  rho.S.star <- max((var.m.S - var.m.S_j) / (var.m.S + pi^2 / 3))

  # found the conditional residual error of Y
  fit.Y0 <- glm(as.formula(paste0("Y~", paste0(predictors.OR, collapse = "+"))),
    data = sim_data,
    subset = S == 0
  )
  sigma.Y0 <- summary(fit.Y0)$dispersion

  fit.Y1 <- glm(as.formula(paste0("Y~", paste0(predictors.OR, collapse = "+"))),
    data = sim_data,
    subset = S == 1
  )
  sigma.Y1 <- summary(fit.Y1)$dispersion

  fit.Y <- glm(as.formula(paste0("Y~S+", paste0(predictors.OR, collapse = "+"))),
    data = sim_data
  )
  sigma.Y <- summary(fit.Y)$dispersion

  gamma_R0.t.star <- nleqslv::nleqslv(
    x = 0,
    fn = function(gamma_R) {
      sigma.Y0^2 * gamma_R^2 /
        {
          var.m.R0 + pi^2 / 3 + sigma.Y0^2 * gamma_R^2
        } -
        rho.R0.star
    }
  )$x

  gamma_R1.t.star <- nleqslv::nleqslv(
    x = 0,
    fn = function(gamma_R) {
      sigma.Y1^2 * gamma_R^2 /
        {
          var.m.R1 + pi^2 / 3 + sigma.Y1^2 * gamma_R^2
        } -
        rho.R1.star
    }
  )$x

  gamma_S.t.star <- nleqslv::nleqslv(
    x = 0,
    fn = function(gamma_S) {
      sigma.Y^2 * gamma_S^2 /
        {
          var.m.S + pi^2 / 3 + sigma.Y^2 * gamma_S^2
        } -
        rho.S.star
    }
  )$x

  list(
    est = data.frame(tau = c(
      tau.OM = tau.OM,
      tau.PS = tau.PS,
      tau.hat.J2R = tau.hat.J2R,
      tau.hat.t2t = tau.hat.t2t
    )),
    gamma = c(
      gamma_R1 = gamma_R1, gamma_R0 = gamma_R0,
      gamma_S = gamma_S
    ),
    gamma.t_star = c(
      gamma_R1 = gamma_R1.t.star,
      gamma_R0 = gamma_R0.t.star,
      gamma_S.t.star = gamma_S.t.star
    )
  )
}
