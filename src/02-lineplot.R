# source(here::here('PowerCalc/loadproject.R'))
# 
# sims <- readRDS(file.path("U:", "Projects", "PCI Paper", "PowerCalc", "save", "mde.1K.Rds"))

sims <- readRDS("../mde_full.Rds")

# Line plot median MDE ----------------------------------------------------

lineplot_data <- sims %>% 
  mutate(bayes = factor(bayes, levels = c(TRUE, FALSE),
                        labels = c('Bayesian', 'Classical'))) %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>%
  group_by(N, nArms, sigfun, bayes, efftype) %>% 
  summarise(MDE = if (mean(is.na(MDE)) > 1/3) NA else 
    median(replace_na(MDE, Inf)))

lineplot <- function(.x) ggplot(data = .x, 
                                aes(x = N, y = MDE)) +  
    geom_line(aes(col=factor(nArms, levels=c(16,72,180,250), 
                             labels=c('16 arms',
                                      '72 arms',
                                      '180 arms',
                                      '250 arms')),
                  linetype = bayes)) +
    geom_point(aes(col=factor(nArms, levels=c(16,72,180,250), 
                              labels=c('16 arms',
                                       '72 arms',
                                       '180 arms',
                                       '250 arms')),
                   shape = bayes))+
    guides(col=guide_legend(title=''), linetype = guide_legend(title = ''),
           shape = guide_legend(title = '')) +
    labs(x='Sample size') +
    theme_bw() + 
    scale_color_discrete(breaks = c('16 arms','72 arms','180 arms','250 arms'))+
    scale_y_continuous(limits = c(0, .7), name = 'Median MDE') + 
    scale_x_continuous(limits = c(0, 10000))

for (.eft in c('main', 'arms')) {
lineplot_data %>% 
  filter(efftype == .eft) %>% 
  lineplot() +
  ggsave(sprintf('graphs/lineplot_bf_%s.png', .eft), 
         width = 6.5, height = 5)
}
  
lineplot_data %>% 
  filter(efftype == 'arms', bayes == TRUE) %>% 
  lineplot() +
  ggsave('graphs/lineplot_arms.png', 
         width = 6.5, height = 5)
