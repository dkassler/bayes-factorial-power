#unloadNamespace("plyr")
library(tidyverse)

raw_bayes <- read_csv("freq/BayesMDE.csv")
MDE_bayes <- raw_bayes %>% transmute(nArms = nA * nB * nC * nD, N, MDE)

raw_freq <- read_csv("data/power/classicalPowerRedux_output.csv")
MDE_freq <- raw_freq %>% 
  replace_na(list(MDE = 1)) %>% 
  group_by(nArms, N) %>% 
  summarise(MDE = if (mean(MDE >= 1) <= 1/3) median(MDE) else NA)

mde <- bind_rows(`Bayesian` = MDE_bayes, `Classical` = MDE_freq, .id = "type")
  
ggplot(data = mde, aes(x = N, y = MDE)) +  
  geom_line(aes(col=factor(nArms, levels=c(16,72,180,250), 
                           labels=c('16 arms',
                                    '72 arms',
                                    '180 arms',
                                    '250 arms')),
                linetype = type)) +
  geom_point(aes(col=factor(nArms, levels=c(16,72,180,250), 
                            labels=c('16 arms',
                                     '72 arms',
                                     '180 arms',
                                     '250 arms')),
                 shape = type))+
  guides(col=guide_legend(title=''), linetype = guide_legend(title = ''),
         shape = guide_legend(title = '')) +
  labs(x='Sample size') +
  theme_bw() + 
  scale_color_discrete(breaks = c('16 arms','72 arms','180 arms','250 arms'))+
  scale_y_continuous(limits = c(0, .7)) + 
  scale_x_continuous(limits = c(0, 10000))

ggsave(file='freq/BayesFreqPowerCompare.png', width=6.5, height=5)
ggsave(file='freq/BayesFreqPowerCompare.pdf', width=6.5, height=5)

