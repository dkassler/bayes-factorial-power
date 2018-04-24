source(here::here('PowerCalc/loadproject.R'))
library(multidplyr)

#' Make data frame of simulation scenarios

study_designs <- list(
  design_frame(2, 2, 2, 2), 
  design_frame(3, 3, 10, 2),
  design_frame(4, 3, 3, 2),
  design_frame(5, 5, 5, 2)
)

# sigfuns <- list(
#   bayes_arms = partial(sigfun.bayes_arms, .lazy = FALSE, model = model, 
#                        control = list(adapt_delta = 0.99, max_treedepth = 20),
#                        chains = 1)
# )
sigfuns <- list(
  'bayes_arms' = sig.bayes_arms,
  'bayes_main' = sig.bayes_main
)

scenarios <- cross_df(list(
  study_design = study_designs,
  N = c(350, 1000, 3500, 5000, 10000),
  # N = c(350, 10000),
  sigfun = names(sigfuns)
)) %>% 
  unnest(study_design) 

#' Run simulation to get MDE for each

sims <- scenarios %>% 
  replicate(n = 1000, expr = ., simplify = F) %>% 
  bind_rows(.id = 'i')

# cl <- makeCluster(1)
# doSNOW::registerDoSNOW(cl)
# on.exit(stopCluster(cl), add = T)

# foreach_args <- c(
#   as.list(sims),
#   .errorhandling = 'pass',
#   .packages = 'tidyverse'
# )
# 
# MDEs <- do.call(foreach, foreach_args) %do% {
#   sigfun <- sigfuns[[sigfun]]
#   simMDE(nA, nB, nC, nD, N, sigfun)
# }
# 
# sims$MDE <- MDEs



# Scratch -----------------------------------------------------------------

closeAllConnections()
cl <- create_cluster(63) %>% 
  # cluster_copy(fit_bayes) %>% 
  cluster_copy(sigfuns) #%>% 
  # cluster_copy(get_mde) %>% 
  # cluster_copy(model)
cluster_eval(cl, source(here::here('PowerCalc/loadproject.R')))
set_default_cluster(cl)

sims1 <- sims %>% 
  group_by(i, nA, nB, nC, nD) %>% 
  mutate(coef = list(gen_true_coef(nA[[1]], nB[[1]], nC[[1]], nD[[1]]))) %>% 
  group_by(N, add = TRUE) %>% 
  mutate(data = list(gen_study_data(coef[[1]], N[[1]]))) %>% 
  partition(i, nA, nB, nC, nD, N) %>% 
  mutate(fit = list(quietly(fit_bayes)(data[[1]], coef[[1]], model, chains = 1, iter = 1000,
                              control = list(adapt_delta = 0.9, max_treedepth = 10))))

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

saveRDS(sims4, 'save/mde.1K.Rds')
