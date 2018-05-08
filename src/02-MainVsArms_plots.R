setwd(here::here())
library(ProjectTemplate)
load.project()

sims <- readRDS(
  file.path("N:", "Project", "50381_PCI-W", "DC1", "Factorial Design paper", 
            "Data", "mde_full.Rds")
)


# Facetted Boxplot --------------------------------------------------------

boxplot_data <- sims %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>% 
  group_by(N, nArms, sigfun) %>% 
  mutate(pct_na = mean(is.na(MDE))) %>% 
  ungroup() %>% 
  replace_na(list(MDE = 1))

boxplot_plot <- ggplot(data = boxplot_data) +
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
boxplot_plot

ggsave('graphs/boxplot_arms_main.png', width = 6, height = 6)
  

# Facetted bar plot -------------------------------------------------------

pointplot_data <- sims %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>% 
  group_by(N, nArms, bayes, efftype, sigfun) %>% 
  summarize(pct_na = mean(is.na(MDE)),
            mean_MDE = mean(MDE, na.rm = TRUE),
            median_MDE = median(replace_na(MDE, Inf))) %>% 
  ungroup() %>% 
  filter(pct_na < 1/3)

pointplot_labeller <- labeller(bayes = c('TRUE' = 'Bayesian',
                                         'FALSE' = 'Frequentist'))

ggplot(data = pointplot_data,
       aes(x = N, 
           color = sigfun,
           shape = efftype,
           linetype = bayes,
           y = median_MDE)) +
  facet_grid(. ~ nArms, labeller = pointplot_labeller) +
  geom_point() +
  geom_line() +
  scale_color_discrete(limits = c('bayes_arms', 'freq_arms', 
                              'bayes_main', 'freq_main'),
                   breaks = c('bayes_arms', 'bayes_main',
                              'freq_arms', 'freq_main'),
                   labels = c('Bayes Arms', 'Bayes Main Eff.',
                              'Freq. Arms', 'Freq. Main Eff.'))


# Combined boxplot --------------------------------------------------------

boxplot_plot +
  geom_step(aes(x = sigfun, y = MDE))
  
