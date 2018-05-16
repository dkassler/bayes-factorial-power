sims <- readRDS("../mde_full.Rds") %>% 
  mutate(nArms = factor(nA * nB * nC * nD))