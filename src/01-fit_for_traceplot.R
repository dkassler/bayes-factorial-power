setwd(here::here())
library(ProjectTemplate)
load.project()

library(multidplyr)

#' Make data frame of simulation scenarios

study_designs <- list(
  design_frame(2, 2, 2, 2), 
  design_frame(3, 4, 3, 2),
  design_frame(4, 3, 5, 3),
  design_frame(5, 5, 5, 2)
)

fitfuns <- list(
  'bayes' = partial(fit_bayes, model = model, chains = 3, iter = 1000,
                    control = list(adapt_delta = 0.9, max_treedepth = 25))
)

scenarios <- cross_df(list(
  study_design = study_designs,
  N = c(1000, 3500, 5000, 10000),
  fitfun = 'bayes'
)) %>% 
  unnest(study_design)

#' Run simulation to get MDE for each

sims <- scenarios %>% 
  replicate(n = 1, expr = ., simplify = F) %>% 
  bind_rows(.id = 'i')

closeAllConnections()
rm(list = 'cl')

cl <- create_cluster() %>% 
  cluster_copy(fitfuns)
cluster_eval(cl, {
  library(ProjectTemplate)
  setwd(here::here())
  load.project()
})
set_default_cluster(cl)

sims <- sims %>% 
  group_by(i, nA, nB, nC, nD) %>% 
  mutate(coef = list(gen_true_coef(nA[[1]], nB[[1]], nC[[1]], nD[[1]]))) %>% 
  group_by(N, add = TRUE) %>% 
  mutate(data = list(gen_study_data(coef[[1]], N[[1]]))) %>% 
  partition(i, nA, nB, nC, nD, N) %>% 
  group_by(fitfun, add = TRUE) %>% 
  mutate(fit = list(quietly(fitfuns[[fitfun[[1]]]])(data[[1]], coef[[1]]))) %>% 
  collect()

saveRDS(sims, 'cache/sims_single.Rds')
