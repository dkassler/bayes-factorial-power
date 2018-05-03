sims <- readRDS(
  file.path("N:", "Project", "50381_PCI-W", "DC1", "Factorial Design paper", 
            "Data", "mde_full.Rds")
)

boxplot_data <- sims %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>% 
  group_by(N, nArms, sigfun) %>% 
  mutate(pct_na = mean(is.na(MDE))) %>% 
  ungroup() %>% 
  replace_na(list(MDE = 1))

ggplot(data = boxplot_data) +
  facet_grid(N ~ nArms) +
  geom_boxplot(aes(x = sigfun, y = MDE, color = pct_na > 1/3)) +
  scale_x_discrete(limits = c('bayes_arms', 'freq_arms', 
                              'bayes_main', 'freq_main'),
                   breaks = c('bayes_arms', 'bayes_main',
                              'freq_arms', 'freq_main'),
                   labels = c('Bayes Arms', 'Bayes Main Eff.',
                              'Freq. Arms', 'Freq. Main Eff.')) +
  scale_color_manual(name = "More than 1/3 NA", breaks = c(FALSE, TRUE),
                     values = c('black', 'brown')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Significance Method")

ggsave('graphs/boxplot_arms_main.png')
  
