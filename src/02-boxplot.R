
sims <- readRDS("../mde_full.Rds")

# Boxplot -----------------------------------------------------------------

boxplot_data <- sims %>% 
  filter(sigfun == 'bayes_arms') %>%
  mutate(nArms = factor(nA * nB * nC * nD)) %>% 
  replace_na(list(MDE = 1))

ggplot(data = boxplot_data, aes(x = factor(N), y = MDE, fill = 
                                  factor(nArms, levels=c(16, 72, 180, 250),
                                         labels=c('2x2x2x2 = 16 arms',
                                                  '3x4x3x2 = 72 arms',
                                                  '4x3x5x3 = 180 arms',
                                                  '5x5x5x2 = 250 arms')))) +
  geom_boxplot() +
  guides(fill=guide_legend(title='')) +
  labs(x='Sample Size') +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1), 
                     labels = c('0.00', '0.25', '0.50', '0.75', expression(phantom(x) >= '1.00'))) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggsave('graphs/boxplot.png', width = 6.5, height = 5)
