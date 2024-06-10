
# install.packages(
#   './',
#   repos = NULL,
#   type  = 'source'
# )

SLURM_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')

# SLURM_ID = 1

library(instrument)

set.seed(as.numeric(SLURM_ID))

# Simulate MIRT parameters given model type
# type - type of mirt model to simulate
# type = c("irt", "mirt", "soirt", "birt")
# "irt"
# "mirt"
# "soirt"
# "birt"
sim_mirt_pars = function() {

  # n = number of observations (sample size)
  n = 400 * 3

  # d = number of dimensions (1st order)
  d = 3

  # j is number of total questions (20 q's per dimension)
  j = 20*d

  lambda_ind = rep(1:3, each = 20) # change this for the bigger model

  # number of response categories
  ncat = 3

  # number of categories per item
  ncategi = c(rep(ncat, j))

  # maximum number of categories
  ncateg_max = max(ncategi)

  # fixed effect variables
  ff = as.matrix(data.frame(x1 = sample(c(0, 1), n, TRUE),
                            x2 = sample(c(0, 1), n, TRUE)))

  # alpha - discrimination parameter (slope from factor analysis)
  alpha = matrix(0, d, j)

  # define the dominant alphas for each dimension
  # alpha_dominant = list(1:5, 6:10, 11:15)
  alpha_dominant = list(1:20, 21:40, 41:60)

  # Dominant alphas follow Unif(1.7, 3.0)
  # non-dominant alphas follow Unif(0.2, 1.0)
  for(dd in 1:d) {
    alpha[dd, alpha_dominant[[dd]]] = sort(runif(length(alpha_dominant[[dd]]), 1.5, 4))
    alpha[dd, setdiff(unlist(alpha_dominant), alpha_dominant[[dd]])] = rep(0.0, length(unlist(alpha_dominant[-dd]))) #sort(runif(length(unlist(alpha_dominant[-dd])), 0.0, 0.4))
  }

  # delta - difficulty parameters
  delta = matrix(nrow = j, ncol = ncateg_max - 1)

  # deltas are ordered and follow N(0, 1)
  for(jj in 1:j) {
    # delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 2.5))
    dval = c()
    for(d_pos in 1:(ncategi[jj] - 1)) {
      dval = c(dval, rnorm(1, -1 + d_pos, 2.5))
    }
    delta[jj, 1:(ncategi[jj]-1)] = sort(dval)
  }

  # first level of deltas are zero for estimability
  delta = cbind(0, delta)

  # theta parameters - ability at individual level
  theta = matrix(0, nrow = n, ncol = d)
  theta_g = matrix(0, nrow = n, ncol = 1)

  # ability follows N(0, 1)
  for(dd in 1:d) {
    theta[, dd] = rnorm(n, 0, 0.1) # change the s.d. of this variable?
  }
  theta_g[, 1] = rnorm(n, 0, 1)

  # beta = c(0.7, -0.5)
  beta = c(0.0, 0.0)

  lambda = c(2, 1.2, -1.8)

  # two correlated random effects (one intercept, one slope)
  n_per_level = 10
  Lz = n / n_per_level
  # sigma = matrix(c(4.0, 0.0, 0.0, 3.8), nrow = 2)
  sigma = 1
  # int_slope = mvtnorm::rmvnorm(Lz, mean = c(0, 0), sigma = sigma)
  ints = rnorm(Lz, 0, sigma)
  # cov(int_slope)
  z_c_int = matrix(rep(diag(Lz), each = n_per_level), nrow = n)
  # zc_is = z_c_int %*% int_slope[, 1]
  # z_c_slope = z_c_int
  # z_c_slope[z_c_slope == 1] = rnorm(n)
  # zc_is = zc_is + z_c_slope %*% int_slope[, 2]
  z_c = z_c_int
  z = data.frame(
      school = rep(paste0("school", 1:Lz), each = n_per_level)
    # pred = z_c_slope[z_c_slope != 0]
    )


  #
  # start_index = 1
  # beta_dstart = NULL
  # beta_dend = NULL

  return(
    list(
      n = n, d = d, j = j, ncat = ncat, ncategi = ncategi,
      ncateg_max = ncateg_max, alpha = alpha, alpha_dominant = alpha_dominant, delta = delta,
      ff = ff, beta = beta, theta = theta, lambda = lambda, theta_g = theta_g,
      lambda_ind = lambda_ind, z = z, z_c = z_c, ints = ints
    )
  )

}

# Given the set of parameters we generated, sample a new data set
# pars = sim_mirt_pars()
sim_mirt_data = function(type, pars) {
  # pars = sim_mirt_pars()
  # pull values from sim_mirt_pars() output
  n = pars$n
  d = pars$d
  j = pars$j
  ncat = pars$ncat
  ncategi = pars$ncategi
  ncateg_max = pars$ncateg_max
  alpha = pars$alpha
  alpha_dominant = pars$alpha_dominant
  delta = pars$delta
  theta = pars$theta
  theta_g = pars$theta_g
  ff = pars$ff
  beta = pars$beta
  lambda = pars$lambda
  lambda_ind = pars$lambda_ind
  z = pars$z
  z_c = pars$z_c
  ints = pars$ints

  # set.seed(seed_replication_level)

  # sample data set
  data = matrix(0, nrow = n, ncol = j)

  # theta_g = 2.0 * theta_g / sd(theta_g)

  z_cm = z_c %*% ints

  # rescale theta values
  for(i in 1:n) {
    tr = 0.0 #z_c %*% ints  #as.vector(beta %*% ff[i, ])
    # tr = 0
    theta[i, ] = theta[i, ] + tr
  }

  # theta mean = 0, sd = 1
  # theta = apply(theta, 2, \(x) {(x - mean(x)) / (sd(x))})

  # for individual i and question j, calculate the probability for responding at
  # each level of the response
  nu = array(0, dim = c(n, j, ncateg_max))
  for(i in 1:n) {
    for(jj in 1:j) {
      # regression equation
      nu[i, jj, ] = sum((t(alpha[lambda_ind[jj], jj]) * (theta[i, lambda_ind[jj]] + lambda[lambda_ind[jj]]*theta_g[i] + z_cm[i,1] ))) - (delta[jj, 1:ncategi[jj]])
    }
  }

  # theta_g = apply(theta_g, 2, \(x) {(x - mean(x)) / (sd(x))})

  # nu = (nu / sd(nu)) * 1.5

  for(i in 1:j) {
    nu[, i, ] = (nu[, i, ] / sd(nu[, i, ])) * 1.5
  }

  sd(nu)

  for(i in 1:n) {
    for(jj in 1:j) {
      # inverse logit
      prb = (1 / (1 + exp(-(nu[i, jj, ]))))
      # first level is 1
      prb[1] = 1.0
      # last level is 0
      prb = c(prb, 0)
      # difference between the two levels
      prb = prb[-length(prb)] - prb[2:length(prb)]
      # sample given probability vector
      data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
    }
  }

  # preview the data
  apply(data, 2, table)

  hist(apply(data, 1, sum))

  # remove empty categories
  collapse_categories = function(x) {
    return(apply(x, 2, \(y){ match(y, sort(unique(y))) }))
  }

  # apply the remove gaps function
  data = collapse_categories(data)

  # preview the data
  apply(data, 2, table)

  # bind sampled data and predictor variables
  data = cbind(data, ff, z)

  # assign variable names to data
  colnames(data) = c(paste0('x', 1:j), paste0('z', 1:ncol(ff)), paste0('school', 1:ncol(z)))

  # store simulated data
  sim_data = list(alpha = alpha, delta = delta, theta = theta, theta_g = theta_g, beta = beta)

  # store data for fitting the model
  fit_data = list(data = data)

  # structure the alpha's in long format for comparison with estiamtes
  alpha = data.table::data.table(true = as.vector(sim_data$alpha))
  # alpha = data.table::data.table(true = log(as.vector(sim_data$alpha)))

  # update alpha parameter names for later use
  a_names = c()
  for(jj in 1:j) {
    for(i in 1:d) {
      a_names = c(a_names, paste0("alpha[", i, ",", jj, "]"))
    }
  }

  # assign new names to long parameter data
  alpha[, parameter := a_names]

  # cor(
  #     sim_data$alpha[1, ],
  #     alpha[
  #       grep('alpha', parameter),
  #     ][
  #       grep('\\[1,', parameter), true
  #     ]
  #   )

  # structure the theta's in long format for comparison with estiamtes
  theta = data.table::data.table(true = as.vector(sim_data$theta))

  # update alpha parameter names for later use
  t_names = c()
  for(dd in 1:d) {
    for(i in 1:n) {
      t_names = c(t_names, paste0("theta[", i, ",", dd, "]"))
    }
  }

  # assign new names to long parameter data
  theta[, parameter := t_names]

  theta_g = data.table::data.table(true = as.vector(sim_data$theta_g))

  t_names = c()
  for(i in 1:n) {
      t_names = c(t_names, paste0("theta_g[", i, "]"))
  }

  theta_g[, parameter := t_names]

  # structure the delta's in long format for comparison with estiamtes
  delta = data.table::data.table(true = as.vector(sim_data$delta[, -1]))

  # update alpha parameter names for later use
  d_names = c()
  for(ncat in 1:(ncateg_max-1)) {
    for(jj in 1:j) {
      d_names = c(d_names, paste0("delta[", jj, ",", ncat, "]"))
    }
  }

  # assign new names to long parameter data
  delta[, parameter := d_names]

  # bind true parameter values together
  sim_data$true = data.table::rbindlist(list(alpha, delta, theta, theta_g))

  sim_data$ints = ints

  # cor(
  #     sim_data$alpha[3, ],
  #     sim_data$true[
  #       grep('alpha', parameter),
  #     ][
  #       grep('\\[3,', parameter), true
  #     ]
  #   )

  # list of return values
  return(list(sim_data = sim_data, fit_data = fit_data))

}

sim_data = sim_mirt_data(pars = sim_mirt_pars())

model = 'theta1 = c(1:20)
         theta2 = c(21:40)
         theta3 = c(41:60)
         thetag = theta1 + theta2 + theta3
         thetag ~ (1|school1)'

# model = 'theta1 = c(1:30)
#          theta2 = c(1:30)
#          theta3 = c(1:30)'

fit =
  instrument::instrument(
    data          = sim_data$fit_data[['data']],
    model         = model,
    iter_warmup   = 100,
    iter_sampling = 100,
    seed          = 12345,
    fweights = rep(1, nrow(sim_data$fit_data[['data']]))
  )

# debugging without devtools
library(stringr)
sourceDir = function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('./instrument2/R/')
# data = sim_data$fit_data[['data']]
# model = model
# iter_warmup = 10
# iter_sampling = 20
# seed = 12345

# evaluate the model by pulling different parameters and correlating them, looking at them

# sfit = summary.instrumentObj(fit)

# sfit = object$stanfit$summary()

library(stringr)
library(rstan)
library(tidyverse)

draws = summary.instrumentObj(fit)

output =
  list(
    sim_data = sim_data,
    fit      = fit,
    draws    = draws
  )

save(
  output,
  file = paste0('hoirt_rand_', SLURM_ID, '.rda')
)

#
# # go to summary.instrumentObj
#
# draws[
#   grep('alpha', parameter),
# ][
#   , mean
# ] %>%
#   matrix(., byrow = FALSE, nrow = 90/3) %>%
#   as.data.frame() %>%
#   mutate_all(\(x) { ifelse(x < 0.4, '--', as.character(round(x, digits = 1)))})
#
# sim_data$sim_data$alpha %>%
#   t() %>%
#   matrix(., byrow = FALSE, nrow = 90/3) %>%
#   as.data.frame() %>%
#   mutate_all(\(x) { ifelse(x < 0.4, '--', as.character(round(x, digits = 1)))})
#
# alpha_est =
#   draws[
#     grep('alpha', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = FALSE, nrow = 90/3)
#
# alpha_true =
#   sim_data$sim_data$alpha %>%
#   t()
#
# # cor(
# #   alpha_est[,1],
# #   exp(log(alpha_true[,1] + 1)) - 1
# # )
#
#
#
# # here right now!
#
# mean(exp(alpha_est) -1  - alpha_true)
# plot(
#   exp(alpha_est[,1]) - 1,
#   alpha_true[,3]
# )
# plot(
#   exp(alpha_est[,2]) - 1,
#   alpha_true[,2]
# )
# plot(
#   exp(alpha_est[,3]) - 1,
#   alpha_true[,1]
# )
#
# sqrt(
#   sum(
#     (exp(alpha_est[,3]) - 1 - (alpha_true[,1]))^2
#   )/30
# )
#
# sqrt(sum((alpha_est[,1] - (exp(log(alpha_true[,3] + 1)) - 1))^2)/30)
#
#
# plot(
#   alpha_est[,3],
#   exp((log(alpha_true[,2] + 1)) ) - 1.0
# )
# plot(
#   alpha_est[,2],
#   exp((log(alpha_true[,1] + 1)) ) - 1.0
# )
# plot(
#   alpha_est[,1],
#   exp((log(alpha_true[,3] + 1)) ) - 1.0
# )
#
# mean(alpha_est[,1] - (exp(log(alpha_true[,3] + 1)) - 1))
#
# cbind(alpha_est[,1], (exp(log(alpha_true[,2] + 1)) - 1)) %>% View()
#
# sqrt(sum((alpha_est[,1] - (exp(log(alpha_true[,3] + 1)) - 1))^2)/30)
# sqrt(sum((alpha_est[,2] - (exp(log(alpha_true[,3] + 1)) - 1.0))^2)/30)
# sqrt(sum((alpha_est[,3] - (exp(log(alpha_true[,1] + 1)) - 1))^2)/30)
#
# mean((alpha_est - 1) - (exp(log(alpha_true + 1)) - 1))
#
# mean(alpha_est - alpha_true)
#
# cor(
#   alpha_est[,1],
#   exp(log(alpha_true[,3] + 1.0))
# )
#
# cor(
#   exp(exp(alpha_est[,2])),
#   exp(log(alpha_true[,2] + 1.0))
# )
#
# cor(
#   exp(exp(alpha_est[,3])),
#   exp(log(alpha_true[,1] + 1.0))
# )
#
# plot(
#   alpha_est[,1],
#   log(alpha_true[,1] + 1)
# )
#
# plot(
#   exp(alpha_est[,2]),
#   exp(log(alpha_true[,2] + 1))
# )
#
# plot(
#   exp(alpha_est[,3]),
#   log(alpha_true[,3] + 1)
# )
#
# plot(
#   alpha_est[,2],
#   alpha_true[,1]
# )
#
# plot(
#   alpha_est[,2],
#   alpha_true[,2]
# )
#
# plot(
#   alpha_est[,3],
#   alpha_true[,3]
# )
#
# mean(alpha_est - alpha_true)
#
# plot(
#   alpha_est[,1],
#   log(alpha_true[,1] + 1.0)
# )
#
# plot(
#   alpha_est[,2],
#   log(alpha_true[,2] + 1.0)
# )
#
# plot(
#   exp(alpha_est[,3]) - 1.0,
#   alpha_true[,3]
# )
#
# cor(
#   exp(alpha_est[,1]) - 1.0,
#   alpha_true[,1]
# )
#
# plot(
#   exp(alpha_est[,2]) - 1.0,
#   alpha_true[,2]
# )
#
# mean(exp(alpha_est[,2]) - 1.0 - alpha_true[,2])
#
# cor(
#   exp(alpha_est[,3]) - 1.0,
#   alpha_true[,3]
# )
#
# mean(alpha_est - log(alpha_true + 1.0))
#
# cor(
#   alpha_est[,1],
#   log(alpha_true[,1] + 1)
# )
#
# delta_est =
#   draws[
#     grep('delta', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = TRUE, nrow = 90/3)
#
# delta_true = sim_data$sim_data$delta[,-1]
#
# cor(
#   delta_est[,1],
#   delta_true[,1]
# )
#
# plot(
#   delta_est[,1],
#   delta_true[,1]
# )
#
# cor(
#   delta_est[-1,2],
#   delta_true[-1,2]
# )
#
# plot(
#   delta_est[,2],
#   delta_true[,2]
# )
#
# cor(
#   delta_est[,2],
#   delta_true[,2]
# )
#
# mean(delta_est - delta_true)
#
# theta_est =
#   draws[
#     grep('theta', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = TRUE, nrow = 600)
#
# theta_true =
#   sim_data$sim_data$theta
#
# cor(
#   theta_est[,3],
#   theta_true[,2]
# )
#
# draws[
#   grep('beta', parameter),
# ][
#   , mean
# ]
#
# sim_data$sim_data$beta
#
#
#
#
#
#
#
#
#
# # theta_g
#
#
# theta_g_est =
#   draws[
#     grep('theta\\[', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = TRUE, nrow = 600)
#
# cor(theta_g_est, sim_data$sim_data$theta_g)
# plot(theta_g_est, sim_data$sim_data$theta_g)
#
# lambda_est =
#   draws[
#     grep('lambda_identify\\[', parameter),
#   ][
#     , mean
#   ]
#
# theta_est =
#   draws[
#     grep('theta_resid\\[', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = TRUE, nrow = 600)
#
# lambda_est*sd(theta_g_est) /
#   apply(
#     MARGIN = 2,
#     theta_g_est %*% lambda_est + theta_est,
#     \(x) {
#       sd(x)
#     }
#   )
#
# draws[
#   grep('lambda_identify\\[', parameter),
# ]
#
# draws[
#   grep('sig_thetag_reg', parameter),
# ]
#
# cor(
#  theta_est[,2],
#  sim_data$sim_data$theta[,2]
# )
#
# cor(
#   theta_g_est*lambda_est[1] + theta_est[,1],
#   sim_data$sim_data$theta_g*(0.7) + sim_data$sim_data$theta[,1]
#   )
#
# cor(
#   theta_g_est*lambda_est[2] + theta_est[,2],
#   sim_data$sim_data$theta_g*(0.5) + sim_data$sim_data$theta[,2]
# )
#
# cor(
#   theta_g_est*lambda_est[3] + theta_est[,3],
#   sim_data$sim_data$theta_g*(-0.4) + sim_data$sim_data$theta[,3]
# )
#
# c(0.7, 0.5, -0.4)
#
# # why is the negative not getting picked up???
#
# # alpha
# alpha_est =
#   draws[
#     grep('alpha', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = FALSE, nrow = 90/3)
#
# alpha_true =
#   sim_data$sim_data$alpha %>%
#   t()
#
# draws[
#   grep('alpha\\[3,', parameter),
# ][
#   , mean
# ]
#
# # check out the alpha parameters
#
# # alpha parameters look wrong!
# cor(alpha_est[1:10,1], alpha_true[1:10,1])
# cor(alpha_est[11:20,2], alpha_true[11:20,2])
# cor(alpha_est[21:30,3], alpha_true[21:30,3])
#
#
# plot(alpha_est[1:10,1], alpha_true[1:10,1])
# plot(alpha_est[11:20,2], alpha_true[11:20,2])
# plot(alpha_est[21:30,3], alpha_true[21:30,3])
#
# plot(alpha_est[,1], alpha_true[,1])
# plot(alpha_est[,2], alpha_true[,2])
# plot(alpha_est[,3], alpha_true[,3])
#
# cor(
#   c(alpha_est[1:10,1], alpha_est[11:20,2], alpha_est[21:30,3]),
#   c(alpha_true[1:10,1], alpha_true[11:20,2], alpha_true[21:30,3])
# )
#
#
# delta_est =
#   draws[
#     grep('delta', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = TRUE, nrow = 90/3)
#
# delta_true = sim_data$sim_data$delta[,-1]
#
# cor(
#   delta_est[,1],
#   delta_true[,1]
# )
#
# plot(
#   delta_est[,1],
#   delta_true[,1]
# )
#
# cor(
#   delta_est[-1,2],
#   delta_true[-1,2]
# )
#
# plot(
#   delta_est[,2],
#   delta_true[,2]
# )
#
# cor(
#   delta_est[,2],
#   delta_true[,2]
# )
#
#
# # delta
# delta_est =
#   draws[
#     grep('delta', parameter),
#   ][
#     , mean
#   ] %>%
#   matrix(., byrow = TRUE, nrow = 90/3)
#
# delta_true = sim_data$sim_data$delta[,-1]
#
# cor(
#   delta_est[,1],
#   delta_true[,1]
# )
#
# plot(
#   delta_est[,1],
#   delta_true[,1]
# )
#
# cor(
#   delta_est[-1,2],
#   delta_true[-1,2]
# )
#
# plot(
#   delta_est[,2],
#   delta_true[,2]
# )
#
# cor(
#   delta_est[,2],
#   delta_true[,2]
# )
#
#
#
#
#
#
#
#
# plot(
#   sim_data$sim_data$ints,
#   draws[
#     grep('zeta\\[', parameter),
#   ][
#     , mean
#   ]
# )
#
# cor(
#   sim_data$sim_data$ints,
#   draws[
#     grep('zeta\\[', parameter),
#   ][
#     , mean
#   ]
# )
#
#
#
# draws[
#   grep('zeta_l_sd\\[', parameter),
# ]
#
#
#
#
#
#
#
#
