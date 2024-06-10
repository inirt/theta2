
install.packages(
  './',
  repos = NULL,
  type  = 'source'
  )

library(instrument)

# Simulate MIRT parameters given model type
# type - type of mirt model to simulate
# type = c("irt", "mirt", "soirt", "birt")
# "irt"
# "mirt"
# "soirt"
# "birt"
sim_mirt_pars = function() {

  # n = number of observations (sample size)
  n = 200 * 3

  # d = number of dimensions (1st order)
  d = 3

  # j is number of total questions (20 q's per dimension)
  j = 10*d

  # number of response categories
  ncat = 3

  # number of categories per item
  ncategi = c(rep(ncat, j))

  # maximum number of categories
  ncateg_max = max(ncategi)

  # alpha - discrimination parameter (slope from factor analysis)
  alpha = matrix(0, d, j)

  # define the dominant alphas for each dimension
  # alpha_dominant = list(1:5, 6:10, 11:15)
  alpha_dominant = list(1:10, 11:20, 21:30)

  # Dominant alphas follow Unif(1.7, 3.0)
  # non-dominant alphas follow Unif(0.2, 1.0)
  for(dd in 1:d) {
    alpha[dd, alpha_dominant[[dd]]] = sort(runif(length(alpha_dominant[[dd]]), 1.5, 4))
    alpha[dd, setdiff(unlist(alpha_dominant), alpha_dominant[[dd]])] = sort(runif(length(unlist(alpha_dominant[-dd])), 0.0, 0.4))
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

  # ability follows N(0, 1)
  for(dd in 1:d) {
    theta[, dd] = rnorm(n, 0, 1)
  }

  #
  # start_index = 1
  # beta_dstart = NULL
  # beta_dend = NULL

  return(
    list(
      n = n, d = d, j = j, ncat = ncat, ncategi = ncategi,
      ncateg_max = ncateg_max, alpha = alpha, alpha_dominant = alpha_dominant, delta = delta,
      theta = theta
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
  # set.seed(seed_replication_level)

  # sample data set
  data = matrix(0, nrow = n, ncol = j)

  # rescale theta values

  # theta mean = 0, sd = 1
  theta = apply(theta, 2, \(x) {(x - mean(x)) / (sd(x))})

  # for individual i and question j, calculate the probability for responding at
  # each level of the response
  nu = array(0, dim = c(n, j, ncateg_max))
  for(i in 1:n) {
    for(jj in 1:j) {
      # regression equation
      nu[i, jj, ] = sum((t(alpha[, jj]) %*% (theta[i, ]))) - (delta[jj, 1:ncategi[jj]])
    }
  }

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

  # assign variable names to data
  colnames(data) = c(paste0('x', 1:j))

  # store simulated data
  sim_data = list(alpha = alpha, delta = delta, theta = theta)

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
  sim_data$true = data.table::rbindlist(list(alpha, delta, theta))

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

model = 'theta1 = c(1:30)
         theta2 = c(1:30)
         theta3 = c(1:30)'

fit =
  instrument::instrument(
    data = sim_data$fit_data[['data']],
    model = model,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 12345
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
sourceDir('./R/')
data = sim_data$fit_data[['data']]
model = model
iter_warmup = 10
iter_sampling = 20
seed = 12345

# evaluate the model by pulling different parameters and correlating them, looking at them

# sfit = summary.instrumentObj(fit)

# sfit = object$stanfit$summary()

library(stringr)
library(rstan)
library(tidyverse)

draws = summary.instrumentObj(fit)

# go to summary.instrumentObj

draws[
    grep('alpha', parameter),
  ][
    , mean
  ] %>%
  matrix(., byrow = FALSE, nrow = 90/3) %>%
  as.data.frame() %>%
  mutate_all(\(x) { ifelse(x < 0.4, '--', as.character(round(x, digits = 1)))})

sim_data$sim_data$alpha %>%
  t() %>%
  matrix(., byrow = FALSE, nrow = 90/3) %>%
  as.data.frame() %>%
  mutate_all(\(x) { ifelse(x < 0.4, '--', as.character(round(x, digits = 1)))})

alpha_est =
  draws[
    grep('alpha', parameter),
  ][
    , mean
  ] %>%
  matrix(., byrow = FALSE, nrow = 90/3)

alpha_true =
  sim_data$sim_data$alpha %>%
  t()

cor(
  alpha_est[,3],
  alpha_true[,2]
)

delta_est =
  draws[
    grep('delta', parameter),
  ][
    , mean
  ] %>%
  matrix(., byrow = TRUE, nrow = 90/3)

delta_true = sim_data$sim_data$delta[,-1]

cor(
  delta_est[,1],
  delta_true[,1]
)

cor(
  delta_est[,2],
  delta_true[,2]
)

theta_est =
  draws[
    grep('theta', parameter),
  ][
    , mean
  ] %>%
  matrix(., byrow = TRUE, nrow = 600)

theta_true =
  sim_data$sim_data$theta

cor(
  theta_est[,2],
  theta_true[,3]
)















