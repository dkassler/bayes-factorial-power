setwd(here::here())
library(ProjectTemplate)
load.project()

library(multidplyr)

#' Make data frame of simulation scenarios

study_designs <- list(
  design_frame(2, 2, 2, 2), 
  design_frame(3, 3, 5, 4),
  design_frame(4, 3, 3, 2),
  design_frame(5, 5, 5, 2)
)

sigfuns <- list(
  'bayes_arms' = sig_bayes_arms,
  'bayes_main' = sig_bayes_main,
  'freq_arms'  = sig_freq_arms,
  'freq_main'  = sig_freq_main 
)

fitfuns <- list(
  'bayes' = partial(fit_bayes, model = model, chains = 1, iter = 1000,
                    control = list(adapt_delta = 0.9, max_treedepth = 25)),
  'freq_inter' = fit_freq_inter,
  'freq_main' = fit_freq_main
)

cfgs <- list(
  fit_sig_type(T, 'arms', 'bayes', 'bayes_arms'),
  fit_sig_type(T, 'main', 'bayes', 'bayes_main'),
  fit_sig_type(F, 'arms', 'freq',  'freq_arms'),
  fit_sig_type(F, 'main', 'freq',  'freq_main')
)

scenarios <- cross_df(list(
  study_design = study_designs,
  N = c(350, 1000, 3500, 5000, 10000),
  cfg = cfgs
)) %>% 
  unnest(study_design, cfg) 

#' Run simulation to get MDE for each

sims <- scenarios %>% 
  replicate(n = 1000, expr = ., simplify = F) %>% 
  bind_rows(.id = 'i')

closeAllConnections()
cl <- create_cluster() %>% 
  cluster_copy(sigfuns)
cluster_eval(cl, {
  library(ProjectTemplate)
  setwd(here::here())
  load.project()
})
set_default_cluster(cl)

sims1 <- sims %>% 
  group_by(i, nA, nB, nC, nD) %>% 
  mutate(coef = list(gen_true_coef(nA[[1]], nB[[1]], nC[[1]], nD[[1]]))) %>% 
  group_by(N, add = TRUE) %>% 
  mutate(data = list(gen_study_data(coef[[1]], N[[1]]))) %>% 
  partition(i, nA, nB, nC, nD, N) %>% 
  mutate(fit = list(quietly(fitfuns[[fitfun[[1]]]])(data[[1]], coef[[1]])))

sims2 <- sims1 %>% 
  group_by(sigfun, add = TRUE) %>% 
  mutate(PPs = list(quietly(sigfuns[[sigfun[[1]]]])(coef[[1]], fit[[1]]$result)))

sims3 <- sims2 %>%
  mutate(MDE = map(PPs, ~ quietly(get_mde)(.x$result))) %>% 
  collect()

sims4 <- sims3 %>% 
  mutate(fit.warnings = map(fit, 'warnings'),
         PPs.warnings = map(PPs, 'warnings'),
         MDE.warnings = map(MDE, 'warnings')) %>% 
  mutate(MDE = map_dbl(MDE, 'result')) %>% 
  mutate(hmc_diagnostics = 
           map(fit, ~ get_sampler_params(.x[['result']], inc_warmup = F) %>% 
                 do.call(rbind, .) %>% colMeans %>% as.data.frame.list)) %>% 
  select(-fit, -PPs, -data) %>% 
  ungroup

saveRDS(sims4, 'save/mde_full.Rds')
