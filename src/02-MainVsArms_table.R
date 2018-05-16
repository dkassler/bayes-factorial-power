# sims <- readRDS(
#   file.path("N:", "Project", "50381_PCI-W", "DC1", "Factorial Design paper", 
#             "Data", "mde_full.Rds")
# )

# sims <- readRDS("../mde_full.Rds")

tab <- sims %>% 
  mutate(nArms = nA * nB * nC * nD) %>% 
  group_by(sigfun, N, nArms) %>% 
  summarize(
    mean_MDE = mean(MDE, na.rm = TRUE),
    median_MDE = median(replace_na(MDE, Inf)),
    sd_MDE = sd(MDE, na.rm = TRUE),
    pct_na = mean(is.na(MDE))
  ) %>% 
  mutate(median_MDE = ifelse(pct_na > 1/3, NA, median_MDE),
         sd_MDE     = ifelse(pct_na > 1/3, NA, sd_MDE),
         mean_MDE   = ifelse(pct_na > 1/3, NA, mean_MDE)) %>% 
  filter(N > 350) %>% 
  select(sigfun, N, nArms, median_MDE) %>% 
  spread(key = sigfun, value = median_MDE)
