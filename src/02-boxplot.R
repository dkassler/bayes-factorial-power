
sims <- readRDS("../mde_full.Rds")

# Boxplot -----------------------------------------------------------------

boxplot_data <- sims %>% 
  filter(sigfun == 'bayes_arms') %>%
  mutate(nArms = factor(nA * nB * nC * nD)) %>% 
  replace_na(list(MDE = 1))

ggplot(data = boxplot_data) +
  # facet_wrap(~sigfun) + 
  geom_boxplot(aes(x = factor(N), fill = nArms, y = MDE)) +
  theme_bw() +
  ggsave('graphs/boxplot.png')