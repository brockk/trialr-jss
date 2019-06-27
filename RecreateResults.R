
library(trialr)


# CRM ----
skeleton <- c(0.05, 0.12, 0.25, 0.40, 0.55)
target <- 0.25
fit <- stan_crm(outcome_str = '3N 5N 5T 3N 4N',
                skeleton = skeleton, target = target,
                model = 'logistic', a0 = 3,
                beta_mean = 0, beta_sd = sqrt(1.34),
                seed = 123, refresh = 0)

fit



library(dplyr)
library(tidybayes)
fit %>% 
  gather_draws(prob_tox[dose]) %>% 
  head




# Figure 1
fit %>% 
  gather_draws(prob_tox[dose]) %>% 
  group_by(.draw) %>% 
  summarise(mtd = dose[which.min(abs(.value - target))]) -> mtd_candidates
library(ggplot2)
fit %>% 
  gather_draws(prob_tox[dose]) %>% 
  left_join(mtd_candidates, by = '.draw') %>% 
  filter(.draw <= 200) %>% 
  ggplot(aes(x = dose, y = .value, group = .draw)) +
  geom_line(aes(col = as.factor(mtd)), alpha = 0.5) + 
  geom_hline(yintercept = target, col = 'red', linetype = 'dashed') + 
  labs(title = 'The identify of the MTD is still shrouded in mystery', 
       y = 'Prob(DLT)', col = 'MTD') +
  theme(legend.position = 'bottom')


mtd_candidates %>% 
  count(mtd) %>% 
  mutate(prob_mtd = n / sum(n))


outcomes <- '2NN 3TN'
skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
target <- 0.25
paths1 <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric', 
                   cohort_sizes = c(3, 3), previous_outcomes = outcomes,
                   beta_sd = 1, refresh = 0)
paths2 <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric',
                   cohort_sizes = c(3, 3), previous_outcomes = outcomes,
                   user_dose_func = function(x) {
                     careful_escalation(x, tox_threshold = target + 0.1, 
                                        certainty_threshold = 0.7,
                                        reference_dose = 1)
                   }, beta_sd = 1, seed = 123, refresh = 0)



library(tibble)
paths2_df <- as_tibble(paths2) 
paths2_df %>% head()



# Figure 2
library(tibble)
library(dplyr)
paths2_df <- as_tibble(paths2)

library(DiagrammeR)
# DiagrammeR requires data.frames of nodes and edges to create a graph.
paths2_df %>%
  transmute(id = .node,
            type = NA,
            label = case_when(
              is.na(next_dose) ~ 'Stop',
              TRUE ~ next_dose %>% as.character()),
            shape = 'circle',
            fillcolor = case_when(
              next_dose == 1 ~ 'slategrey',
              next_dose == 2 ~ 'skyblue1',
              next_dose == 3 ~ 'royalblue1',
              next_dose == 4 ~ 'orchid4',
              next_dose == 5 ~ 'royalblue4',
              is.na(next_dose) ~ 'red'
            )
  ) -> nodes_df
paths2_df %>% 
  filter(!is.na(.parent)) %>% 
  select(from = .parent, to = .node, label = outcomes) %>% 
  mutate(rel = "leading_to") -> edges_df
graph <- create_graph(nodes_df = nodes_df, edges_df = edges_df)
render_graph(graph)



spread_paths(paths2_df %>% select(-fit, -parent_fit, -dose_index))


fit <-stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
               doses_given = c(3, 3, 3, 3),
               tox = c(0, 0, 0, 0),
               weights = c(73, 66, 35, 28) / 126,
               model = 'empiric', beta_sd = sqrt(1.34), 
               seed = 123, refresh = 0)
fit$recommended_dose


# EffTox ----
library(purrr)
library(tidyr)

fit_thall_2014 <- function(outcomes, alpha_sd, beta_sd, gamma_sd, zeta_sd) {
  stan_efftox(outcome_str = outcomes, 
              real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
              efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
              p_e = 0.1, p_t = 0.1, eff0 = 0.5, tox1 = 0.65,
              eff_star = 0.7, tox_star = 0.25,
              alpha_mean = -7.9593, alpha_sd = alpha_sd,
              beta_mean = 1.5482, beta_sd = beta_sd,
              gamma_mean = 0.7367, gamma_sd = gamma_sd,
              zeta_mean = 3.4181, zeta_sd = zeta_sd,
              eta_mean = 0, eta_sd = 0.2,
              psi_mean = 0, psi_sd = 1, 
              seed = 123, refresh = 0)
}


# Figure 3
expand.grid(
  alpha_sd = c(2, 3, 4), 
  beta_sd = c(2, 3, 4), 
  gamma_sd = 2.5423,
  zeta_sd = 2.4406
) %>% 
  as_tibble() %>% 
  mutate(fit = pmap(list(alpha_sd, beta_sd, gamma_sd, zeta_sd), 
                    fit_thall_2014, outcomes = ''),
         series = rownames(.)) -> prior_fits1

prior_fits1 %>% 
  mutate(
    dose = map(fit, 'dose_indices'),
    prob_obd = map(fit, 'prob_obd'),
    entropy = map_dbl(fit, 'entropy')
  ) %>% 
  select(-fit) %>% 
  unnest %>% 
  ggplot(aes(x = dose, y = prob_obd, fill = entropy)) + 
  geom_col() + 
  facet_grid(~ alpha_sd ~ beta_sd) + 
  labs(y = 'Prob(dose has maximal utility)', 
       title = 'Prior SDs of alpha (rows) and beta (columns)') + 
  theme(legend.position = 'bottom')


# Figure 4
expand.grid(alpha_sd = 3.5487, 
            beta_sd = 3.5018,
            gamma_sd = c(2, 3, 4),
            zeta_sd = c(2, 3, 4)
) %>% 
  as_tibble() %>% 
  mutate(fit = pmap(list(alpha_sd, beta_sd, gamma_sd, zeta_sd), 
                    fit_thall_2014, outcomes = ''),
         series = rownames(.)) -> prior_fits2

prior_fits2 %>% 
  mutate(
    dose = map(fit, 'dose_indices'),
    prob_obd = map(fit, 'prob_obd'),
    entropy = map_dbl(fit, 'entropy')
  ) %>% 
  select(-fit) %>% 
  unnest %>% 
  ggplot(aes(x = dose, y = prob_obd, fill = entropy)) + 
  geom_col() + 
  facet_grid(~ gamma_sd ~ zeta_sd) + 
  labs(y = 'Prob(dose has maximal utility)', 
       title = 'Effect of prior SD of gamma (rows) and zeta (columns) on identity of optimal dose') + 
  theme(legend.position = 'bottom')


outcomes <- '1NNN 2ENN'

fit <- stan_efftox(outcome_str = outcomes,
                   real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                   efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                   p_e = 0.1, p_t = 0.1, eff0 = 0.5, tox1 = 0.65,
                   eff_star = 0.7, tox_star = 0.25,
                   alpha_mean = -7.9593, alpha_sd = 3.5487,
                   beta_mean = 1.5482, beta_sd = 3.5018,
                   gamma_mean = 0.7367, gamma_sd = 2.5423,
                   zeta_mean = 3.4181, zeta_sd = 2.4406,
                   eta_mean = 0, eta_sd = 0.2,
                   psi_mean = 0, psi_sd = 1, 
                   seed = 123, refresh = 0)
fit


# Figure 5
efftox_contour_plot(fit)


efftox_superiority(fit)


paths <- efftox_dtps(cohort_sizes = c(3), previous_outcomes = outcomes, 
                     real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                     efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                     p_e = 0.1, p_t = 0.1, eff0 = 0.5, tox1 = 0.65,
                     eff_star = 0.7, tox_star = 0.25,
                     alpha_mean = -7.9593, alpha_sd = 3.5487,
                     beta_mean = 1.5482, beta_sd = 3.5018,
                     gamma_mean = 0.7367, gamma_sd = 2.5423,
                     zeta_mean = 3.4181, zeta_sd = 2.4406,
                     eta_mean = 0, eta_sd = 0.2,
                     psi_mean = 0, psi_sd = 1, 
                     next_dose = 3, seed = 123, refresh = 0)


library(tidyr)
library(purrr)

as_tibble(paths) %>% 
  filter(.depth > 0) %>% 
  mutate(prob_obd = map(fit, 'prob_obd'), 
         parent_prob_obd = map(parent_fit, 'prob_obd')) %>% 
  select(outcomes, dose_index, prob_obd, parent_prob_obd) %>% 
  unnest %>% 
  mutate(prob_obd_delta = prob_obd - parent_prob_obd) %>% 
  filter(dose_index == 5)


# Figure 6 - LHS
paths_df <- as_tibble(paths)
paths_df %>%
  transmute(id = .node,
            type = NA,
            label = case_when(
              is.na(next_dose) ~ 'Stop',
              TRUE ~ next_dose %>% as.character()),
            shape = 'circle',
            fillcolor = case_when(
              next_dose == 1 ~ 'slategrey',
              next_dose == 2 ~ 'skyblue1',
              next_dose == 3 ~ 'royalblue1',
              next_dose == 4 ~ 'orchid4',
              next_dose == 5 ~ 'royalblue4',
              is.na(next_dose) ~ 'red'
            )
  ) -> nodes_df

paths_df %>% 
  filter(!is.na(.parent)) %>% 
  select(from = .parent, to = .node, label = outcomes) %>% 
  mutate(rel = "leading_to") -> edges_df

graph <- create_graph(nodes_df = nodes_df, edges_df = edges_df)
render_graph(graph)


# Figure 6 - RHS
paths <- efftox_dtps(cohort_sizes = c(1, 1), previous_outcomes = outcomes, 
                     next_dose = 3, real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                     efficacy_hurdle = 0.5, toxicity_hurdle = 0.3,
                     p_e = 0.1, p_t = 0.1, eff0 = 0.5, tox1 = 0.65,
                     eff_star = 0.7, tox_star = 0.25,
                     alpha_mean = -7.9593, alpha_sd = 3.5487,
                     beta_mean = 1.5482, beta_sd = 3.5018,
                     gamma_mean = 0.7367, gamma_sd = 2.5423,
                     zeta_mean = 3.4181, zeta_sd = 2.4406,
                     eta_mean = 0, eta_sd = 0.2,
                     psi_mean = 0, psi_sd = 1, 
                     seed = 123, refresh = 2000)
paths_df <- as_tibble(paths)

paths_df %>%
  transmute(id = .node,
            type = NA,
            label = case_when(
              is.na(next_dose) ~ 'Stop',
              TRUE ~ next_dose %>% as.character()),
            shape = 'circle',
            fillcolor = case_when(
              next_dose == 1 ~ 'slategrey',
              next_dose == 2 ~ 'skyblue1',
              next_dose == 3 ~ 'royalblue1',
              next_dose == 4 ~ 'orchid4',
              next_dose == 5 ~ 'royalblue4',
              is.na(next_dose) ~ 'red'
            )
  ) -> nodes_df

paths_df %>% 
  filter(!is.na(.parent)) %>% 
  select(from = .parent, to = .node, label = outcomes) %>% 
  mutate(rel = "leading_to") -> edges_df

graph <- create_graph(nodes_df = nodes_df, edges_df = edges_df)
render_graph(graph)



# Augmented Binary method ----
informative_priors <- list(alpha_mean = 0, alpha_sd = 0.1,
                           beta_mean = 0, beta_sd = 0.1,
                           gamma_mean = 0, gamma_sd = 0.1,
                           sigma_mean = 0, sigma_sd = 0.5,
                           omega_lkj_eta = 1,
                           alpha_d1_mean = 0, alpha_d1_sd = 0.5,
                           gamma_d1_mean = 0, gamma_d1_sd = 0.25,
                           alpha_d2_mean = 0, alpha_d2_sd = 0.5,
                           gamma_d2_mean = 0, gamma_d2_sd = 0.25)

diffuse_priors <- list(alpha_mean = 0, alpha_sd = 1,
                       beta_mean = 0, beta_sd = 1,
                       gamma_mean = 0, gamma_sd = 1,
                       sigma_mean = 0, sigma_sd = 1,
                       omega_lkj_eta = 1,
                       alpha_d1_mean = 0, alpha_d1_sd = 1,
                       gamma_d1_mean = 0, gamma_d1_sd = 1,
                       alpha_d2_mean = 0, alpha_d2_sd = 1,
                       gamma_d2_mean = 0, gamma_d2_sd = 1)



set.seed(123)
diffuse_prior_pred_data <- do.call(prior_predictive_augbin_2t_1a, 
                                   append(diffuse_priors, list(num_samps = 1000)))
inf_prior_pred_data <- do.call(prior_predictive_augbin_2t_1a, 
                               append(informative_priors, 
                                      list(num_samps = 1000)))

library(stringr)
library(tidyr)

# Figure 7
bind_rows(
  diffuse_prior_pred_data %>% mutate(Prior = 'Diffuse'),
  inf_prior_pred_data %>% mutate(Prior = 'Informative')
) %>% 
  select(id, Prior, y0, y1, y2) %>% 
  gather(assessment, y, -id, -Prior) %>% 
  mutate(time = str_extract(assessment, '\\d+') %>% as.integer()) %>% 
  ggplot(aes(x = time, y = y)) + 
  geom_line(aes(group = id), alpha = 0.3, col = 'darkgray') + 
  geom_hline(yintercept = log(0.7), col = 'orange', linetype = 'dashed') + 
  scale_x_continuous(breaks = 0:2) + 
  labs(y = 'log(tumour size ratio)') + 
  facet_wrap(~ Prior)



# Figure 8
bind_rows(
  diffuse_prior_pred_data %>% mutate(Prior = 'Diffuse'),
  inf_prior_pred_data %>% mutate(Prior = 'Informative')
) %>% 
  ggplot(aes(x = z0, y = prob_d1)) + 
  geom_point() + geom_smooth(method = 'gam') +
  labs(x = 'Baseline tumour size (cm)', y = 'Prob(Non-shrikage failure at time 1)') + 
  facet_wrap(~ Prior)



N <- 50
sigma <- 1
delta1 <- -0.356
mu <- c(0.5 * delta1, delta1)
Sigma = matrix(c(0.5 * sigma^2, 0.5 * sigma^2, 
                 0.5 * sigma^2, sigma^2), ncol = 2)
alphaD <- -1.5
gammaD <- 0



set.seed(123456)
y <- MASS::mvrnorm(n = N, mu, Sigma)
z0 <- runif(N, min = 5, max = 10)
z1 <- exp(y[, 1]) * z0
z2 <- exp(y[, 2]) * z0
d1 <- rbinom(N, size = 1, prob = gtools::inv.logit(alphaD + gammaD * z0))
d2 <- rbinom(N, size = 1, prob = gtools::inv.logit(alphaD + gammaD * z1))
tumour_size = data.frame(z0, z1, z2) # cm
non_shrinkage_failure <- data.frame(d1, d2)



fit_diffuse <- stan_augbin(tumour_size, non_shrinkage_failure, 
                           model = '2t-1a', prior_params = diffuse_priors, 
                           seed = 123, refresh = 0)
fit_inf <- stan_augbin(tumour_size, non_shrinkage_failure, 
                       model = '2t-1a', prior_params = informative_priors, 
                       seed = 123, refresh = 0)



pred_diffuse <- predict(fit_diffuse, y2_upper = log(0.7))
pred_diffuse %>% head()



pred_binary <- binary_prob_success(fit_diffuse, y2_upper = log(0.7), 
                                   methods = 'exact')
pred_binary



mean(pred_diffuse$ci_width) / pred_binary$ci_width - 1



# Figure 9
pred_inf <- predict(fit_inf, y2_upper = log(0.7))

pred_diffuse %>% 
  ggplot(aes(x = z0, y = prob_success, col = z1)) + 
  geom_point() + 
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  geom_hline(yintercept = pred_binary$mean, col = 'orange', linetype = 'dashed') +
  geom_hline(yintercept = pred_binary$lower, col = 'orange', linetype = 'dotted') +
  geom_hline(yintercept = pred_binary$upper, col = 'orange', linetype = 'dotted') +
  ylim(0, 0.75) + 
  labs(x = 'Baseline tumour size (cm)', y = 'Prob(Success)', 
       col = 'Interim size (cm)', title = 'Diffuse') + 
  theme(legend.position = 'bottom') -> aubgin_diffuse_plot

pred_inf %>% 
  ggplot(aes(x = z0, y = prob_success, col = z1)) + 
  geom_point() + 
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  geom_hline(yintercept = pred_binary$mean, col = 'orange', linetype = 'dashed') +
  geom_hline(yintercept = pred_binary$lower, col = 'orange', linetype = 'dotted') +
  geom_hline(yintercept = pred_binary$upper, col = 'orange', linetype = 'dotted') +
  ylim(0, 0.75) + 
  labs(x = 'Baseline tumour size (cm)', y = 'Prob(Success)', 
       col = 'Interim size (cm)', title = 'Informative') + 
  theme(legend.position = 'bottom') -> aubgin_inf_plot

gridExtra::grid.arrange(aubgin_diffuse_plot, aubgin_inf_plot, ncol = 2)



predict(fit_diffuse, newdata = data.frame(z0 = 5:10, z1 = 4:9), 
        y2_upper = log(0.7))




