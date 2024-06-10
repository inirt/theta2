library(tidyverse)
library(instrument)
library(gtools)
library(stringr)

# combine runs from two batches and get first 100 with sequential numbering
# 'instrument_n400_p60_2024_04_19_v2/' %>% 
#   list.files() %>% 
#   grep(pattern = 'mirt_fixed_', .) %>% 
#   (
#     \(x) {
#       ('instrument_n400_p60_2024_04_19_v2/' %>% 
#         list.files())[x] %>% 
#         length()
#     }
#   )(.) %>% 
#   (
#     \(x) {
#       file.rename(
#         list.files(path = 'instrument_n400_p60_2024_04_19_v2/', 
#                    pattern = 'mirt_fixed_*',
#                    full.names = TRUE), 
#         paste0("instrument_n400_p60_2024_04_19_v2/mirt_fixed_", 1:x, '.rda')
#       )
#     }
#   )(.)


# file.rename(
#   list.files(path = './', 
#              pattern = 'mirt_fixed_*',
#              full.names = FALSE),
#   paste0('instrument_n400_p60_2024_04_19_v2/',
#          list.files(path = './', 
#                     pattern = 'mirt_fixed_*',
#                     full.names = FALSE)
#          )
# )

# 100-76 # remainder

# 'instrument_n400_p60_2024_04_19_v2_batch2/' %>%
#   list.files() %>%
#   grep(pattern = 'mirt_fixed_*', .) %>%
#   (
#     \(x) {
#       ('instrument_n400_p60_2024_04_19_v2_batch2/' %>%
#         list.files())[x] %>%
#         length()
#     }
#   )(.)
# 
# file.rename(
#   list.files(path = 'instrument_n400_p60_2024_04_19_v2_batch2/', 
#              pattern = 'mirt_fixed_*',
#              full.names = TRUE),
#   paste0("instrument_n400_p60_2024_04_19_v2/mirt_fixed_", 1:30 + 76, '.rda')
#   )

dt = data.frame()

for(j in 1:100) {
  jobid = j
  
  f_name = 
    paste0(
      'instrument_n400_p60_2024_04_19_v2/mirt_fixed_',
      jobid,
      '.rda'
    )
  
  if(file.exists(f_name)) {
    load(
      f_name
    )
  } else {
    next
  }
  
  draws = output$draws
  sim_data = output$sim_data
  
  alpha_est =
    draws[
      grep('alpha', parameter),
    ][
      , mean
    ] %>%
    matrix(., byrow = FALSE, nrow = 60)
  
  alpha_true =
    sim_data$sim_data$alpha %>%
    t()
  
  perms = permutations(3, 3, 1:3)
  diff = rep(0, 6)
  
  for(i in 1:6) {
    estimates = alpha_est
    truth = alpha_true
    ord = perms[i, ]
    estimates = 
      rbind(
        estimates[,ord[1], drop = FALSE],
        estimates[-1,ord[2], drop = FALSE],
        estimates[-c(1:2),ord[3], drop = FALSE]
      )
    truth = 
      rbind(
        truth[,1, drop = FALSE],
        truth[-1,2, drop = FALSE],
        truth[-c(1:2),3, drop = FALSE]
      )
    diff[i] = cor(estimates, truth)
  }
  
  alpha_est = 
    alpha_est[, perms[which.max(diff),]]
  
  alpha_estimates = 
    rbind(
      alpha_est[,1, drop = FALSE],
      alpha_est[-1,2, drop = FALSE],
      alpha_est[-c(1:2),3, drop = FALSE]
    )
  
  alpha_truth = 
    rbind(
      alpha_true[,1, drop = FALSE],
      alpha_true[-1,2, drop = FALSE],
      alpha_true[-c(1:2),3, drop = FALSE]
    )
  
  bias = alpha_estimates - alpha_truth
  
  difference = 
    tibble(
      parm = 'alpha',
      bias = bias[,1]
    )
  
  
  delta_est =
    draws[
      grep('delta', parameter),
    ][
      , mean
    ] %>%
    matrix(., byrow = TRUE, nrow = 60)
  
  delta_true = sim_data$sim_data$delta[, -1]
  
  bias = as.vector(delta_est) - as.vector(delta_true)
  
  difference = 
    bind_rows(
      difference,
      tibble(
        parm = 'delta',
        bias = bias
      )
    )
  
  
  theta_est =
    draws[
      grep('theta', parameter),
    ][
      , mean
    ] %>%
    matrix(., byrow = TRUE, nrow = 600)
  
  theta_est = 
    theta_est[, perms[which.max(diff),]]
  
  theta_true =
    sim_data$sim_data$theta
  
  bias = as.vector(theta_est) - as.vector(theta_true)
  
  difference = 
    bind_rows(
      difference,
      tibble(
        parm = 'theta',
        bias = bias
      )
    )
  
  beta_est = 
    draws[
      grep('beta', parameter),
    ][
      , mean
    ]
  
  beta_est = 
    beta_est[beta_est != 0]
  
  theta_sd = 
    rep(apply(theta_est, 2, sd), each = 2)
  
  beta_est = 
    beta_est / theta_sd
  # 
  # beta_est = 
  #   matrix(beta_est, byrow = FALSE, nrow = 2)
  
  beta_true = 
    rep(sim_data$sim_data$beta, 3)
  
  bias = beta_est - beta_true
  
  difference = 
    bind_rows(
      difference,
      tibble(
        parm = 'beta',
        bias = bias
      )
    )
  
  difference$jobid = jobid
  
  dt = bind_rows(dt, difference)
  
  rm(output, draws, sim_data, alpha_est, alpha_true, perms, diff,
     estimates, truth, ord, bias, difference, delta_est, delta_true,
     beta_est, beta_true, theta_sd, alpha_estimates, alpha_truth,
     theta_est, theta_true, i)
}

write_csv(dt, 'instrument_n400_p60_2024_04_19_v2/sim_out_n400_p60.csv')

dtw = 
  tibble(
    parm = 
      dt %>% 
      filter(jobid == 1) %>% 
      pull(parm)
  )

for(i in 1:100) {
  c_name = paste0('job', i)
  
  dtw = 
    bind_cols(
      dtw,
      tibble(
        {{c_name}} := 
            dt %>% 
            filter(jobid == i) %>% 
            pull(bias)
      )
    )
}
  
mean(dtw[1,-1] %>% as.numeric())

dtw %>% 
  filter(parm == 'alpha') %>% 
  select(-parm) %>% 
  rowMeans() %>% 
  mean()

dtw %>% 
  filter(parm == 'alpha') %>% 
  select(-parm) %>% 
  as.matrix() %>% 
  (
    \(x) {
      y = rep(0, nrow(x))
      for(i in 1:nrow(x)) {
        y[i] = sqrt(mean(x[i, ]^2))
      }
      return(y)
    }
  )(.) %>% 
  mean()
  

dtw %>% 
  filter(parm == 'delta') %>% 
  select(-parm) %>% 
  as.matrix() %>% 
  (
    \(x) {
      x[x > 1e5] = NA
      y = rep(0, nrow(x))
      for(i in 1:nrow(x)) {
        y[i] = mean(x[i, ], na.rm = TRUE)
      }
      y
    }
  )(.) %>% 
  mean()

dtw %>% 
  filter(parm == 'delta') %>% 
  select(-parm) %>% 
  as.matrix() %>% 
  (
    \(x) {
      x[x > 1e5] = NA
      y = rep(0, nrow(x))
      for(i in 1:nrow(x)) {
        y[i] = sqrt(mean(x[i, ]^2, na.rm = TRUE))
      }
      y
    }
  )(.) %>% 
  mean()


dtw %>% 
  filter(parm == 'theta') %>% 
  select(-parm) %>% 
  rowMeans() %>% 
  mean()

dtw %>% 
  filter(parm == 'theta') %>% 
  select(-parm) %>% 
  as.matrix() %>% 
  (
    \(x) {
      y = rep(0, nrow(x))
      for(i in 1:nrow(x)) {
        y[i] = sqrt(mean(x[i, ]^2))
      }
      return(y)
    }
  )(.) %>% 
  mean()





dtw %>% 
  filter(parm == 'beta') %>% 
  select(-parm) %>% 
  rowMeans() %>% 
  mean()

dtw %>% 
  filter(parm == 'beta') %>% 
  select(-parm) %>% 
  as.matrix() %>% 
  (
    \(x) {
      y = rep(0, nrow(x))
      for(i in 1:nrow(x)) {
        y[i] = sqrt(mean(x[i, ]^2))
      }
      return(y)
    }
  )(.) %>% 
  mean()

