
# outcomes <- '1NNN 2NTN 2NNN'
# outcomes <- '1NNN 2NTN'
outcomes <- '2NN 3TN'
skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
target <- 0.25

mod0 <- stan_crm(outcome_str = outcomes, skeleton = skeleton, target = target, 
                 model = 'empiric', beta_sd = 1, seed = 123, refresh = 0)
mod0$recommended_dose  # 2
careful_escalation(mod0, tox_threshold = 0.35, certainty_threshold = 0.7)  # 2

paths1 <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric',
                   cohort_sizes = c(3, 3), previous_outcomes = outcomes,
                   beta_sd = 1, seed = 123, refresh = 0)
paths1[[1]]
library(tibble)
as_tibble(paths1) %>% print(n=100)

paths2 <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric',
                   cohort_sizes = c(3, 3), previous_outcomes = outcomes,
                   user_dose_func = function(x) {
                     careful_escalation(x, tox_threshold = target + 0.1, 
                                        certainty_threshold = 0.7,
                                        reference_dose = 1)
                   }, beta_sd = 1, seed = 123, refresh = 0)
paths2[[1]]
as_tibble(paths2) %>% print(n=100)


library(dplyr)
paths2_df <- as_tibble(paths2)
paths2_df %>% print(n = 100)

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
