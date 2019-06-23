

# Simulation setup ----
library(trialr)
library(dplyr)
library(gtools)

get_data_func <- function() {
  N <- 50
  sigma <- 1
  delta1 <- -0.356;
  mu <- c(0.5 * delta1, delta1)
  Sigma = matrix(c(0.5 * sigma^2, 0.5 * sigma^2, 0.5 * sigma^2, sigma^2), 
                 ncol = 2)
  alphaD <- -1.5
  gammaD <- 0
  y <- MASS::mvrnorm(n = N, mu, Sigma)
  z0 <- runif(N, min = 5, max = 10)
  z1 <- exp(y[, 1]) * z0
  z2 <- exp(y[, 2]) * z0
  d1 <- rbinom(N, size = 1, prob = inv.logit(alphaD + gammaD * z0))
  d2 <- rbinom(N, size = 1, prob = inv.logit(alphaD + gammaD * z1))
  tumour_size = data.frame(z0, z1, z2) # cm
  non_shrinkage_failure <- data.frame(d1, d2)
  list(tumour_size, non_shrinkage_failure)
}

fit_model_func_inf_priors <- function(data) {
  priors <- list(alpha_mean = 0, alpha_sd = 0.1,
                 beta_mean = 0, beta_sd = 0.1,
                 gamma_mean = 0, gamma_sd = 0.1,
                 sigma_mean = 0, sigma_sd = 0.5,
                 omega_lkj_eta = 1,
                 alpha_d1_mean = 0, alpha_d1_sd = 0.5,
                 gamma_d1_mean = 0, gamma_d1_sd = 0.25,
                 alpha_d2_mean = 0, alpha_d2_sd = 0.5,
                 gamma_d2_mean = 0, gamma_d2_sd = 0.25)
  fit <- stan_augbin(tumour_size = data[[1]], 
                     non_shrinkage_failure = data[[2]], 
                     prior_params = priors, 
                     model = '2t-1a', 
                     cores = 4, refresh = 0)
  fit
}

fit_model_func_diffuse_priors <- function(data) {
  priors <- list(alpha_mean = 0, alpha_sd = 1,
                 beta_mean = 0, beta_sd = 1,
                 gamma_mean = 0, gamma_sd = 1,
                 sigma_mean = 0, sigma_sd = 1,
                 omega_lkj_eta = 1,
                 alpha_d1_mean = 0, alpha_d1_sd = 1,
                 gamma_d1_mean = 0, gamma_d1_sd = 1,
                 alpha_d2_mean = 0, alpha_d2_sd = 1,
                 gamma_d2_mean = 0, gamma_d2_sd = 1)
  fit <- stan_augbin(tumour_size = data[[1]], 
                     non_shrinkage_failure = data[[2]], 
                     prior_params = priors, 
                     model = '2t-1a', 
                     cores = 4, refresh = 0)
  fit
}

summarise_func <- function(data, fit) {
  augbin_analysis <- predict(fit, y2_upper = log(0.7))
  augbin_analysis %>% 
    summarise(mean(prob_success)) %>% .[[1]] -> augbin_mean_prob_success
  augbin_analysis %>% 
    summarise(mean(ci_width)) %>% .[[1]] -> augbin_mean_ci_width
  binary_analysis <- binary_prob_success(fit, methods = 'exact')
  
  list(
    augbin_prob_success_lower = augbin_analysis$lower, 
    augbin_prob_success = augbin_analysis$prob_success, 
    augbin_mean_prob_success = augbin_mean_prob_success,
    augbin_prob_success_upper = augbin_analysis$upper,
    augbin_ci_width = augbin_analysis$ci_width, 
    augbin_mean_ci_width = augbin_mean_ci_width,
    bin_prob_success_lower = binary_analysis$lower,
    bin_prob_success = binary_analysis$mean,
    bin_prob_success_upper = binary_analysis$upper,
    bin_ci_width = binary_analysis$ci_width
  )
}

# Run simulation ----

# Batch 1
save_func1 <- function(x) {
  saveRDS(object = x, file = 'Batch1_v2.rds')
}
set.seed(123)
sims <- trialr_simulate(N = 1000, 
                        get_data_func = get_data_func, 
                        fit_model_func = fit_model_func_inf_priors, 
                        summarise_func = summarise_func,
                        num_logs = 20, 
                        save_func = save_func1, num_saves = 5)


# Batch 2
save_func2 <- function(x) {
  saveRDS(object = x, file = 'Batch2.rds')
}

set.seed(123)
sims2 <- trialr_simulate(N = 1000, 
                         get_data_func = get_data_func, 
                         fit_model_func = fit_model_func_diffuse_priors, 
                         summarise_func = summarise_func,
                         num_logs = 20,
                         save_func = save_func2, num_saves = 5)




# Results ----
library(purrr)
truth <- 0.334  # From Table 1 of Wason & Seaman

## Batch 1
# sims1 <- readRDS('Batch1.rds')
sims %>% map_dbl('augbin_mean_ci_width') %>% summary() # Mean = 17.73%
sims %>% map_dbl('bin_ci_width') %>% summary() # Mean = 27.09%
sims %>% map_dbl('augbin_mean_prob_success') %>% summary() # Mean = 28.90%
sims %>% map_dbl('bin_prob_success') %>% summary() # Mean = 32.99%
# Prob(Coverage)
# Bayesian AugBin method
sims %>% 
  map_dfc(.f = function(x) truth < x$augbin_prob_success_upper & 
            truth > x$augbin_prob_success_lower) %>% 
  unlist() %>% mean # 71.70%
# Binary method
sims %>% 
  map_lgl(.f = function(x) truth < x$bin_prob_success_upper & 
            truth > x$bin_prob_success_lower) %>% 
  mean # 96.60%

## Batch 2
# sims2 <- readRDS('Batch2.rds')
sims2 %>% map_dbl('augbin_mean_ci_width') %>% summary() # Mean = 21.98%
sims2 %>% map_dbl('bin_ci_width') %>% summary() # Mean = 27.09%
sims2 %>% map_dbl('augbin_mean_prob_success') %>% summary() # Mean = 31.88%
sims2 %>% map_dbl('bin_prob_success') %>% summary() # Mean = 32.99%
# Prob(Coverage)
# Bayesian AugBin method
sims2 %>% map_dfc(.f = function(x) truth < x$augbin_prob_success_upper & 
                    truth > x$augbin_prob_success_lower) %>% 
  unlist() %>% mean # 92.73%
# Binary method
sims2 %>% 
  map_lgl(.f = function(x) truth < x$bin_prob_success_upper & 
            truth > x$bin_prob_success_lower) %>% mean # 96.60%

