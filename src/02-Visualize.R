# source(here::here('PowerCalc/loadproject.R'))
# 
# sims <- readRDS(file.path("U:", "Projects", "PCI Paper", "PowerCalc", "save", "mde.1K.Rds"))

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

# Line plot median MDE ----------------------------------------------------

lineplot_data <- sims %>% 
  mutate(bayes = factor(bayes, levels = c(TRUE, FALSE),
                        labels = c('Bayesian', 'Classical'))) %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>%
  group_by(N, nArms, sigfun, bayes, efftype) %>% 
  summarise(MDE = if (mean(is.na(MDE)) > 1/3) NA else 
    median(replace_na(MDE, Inf)))

lineplot <- function(.x) ggplot(data = lineplot_data %>% filter(efftype == .x), 
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
    scale_y_continuous(limits = c(0, .7)) + 
    scale_x_continuous(limits = c(0, 10000)) +
  ggsave(sprintf('graphs/lineplot_%s.png', .x), width = 6.5, height = 5)

lineplot1 <- function(.x) lineplot()

lineplot('arms')
lineplot('main')
