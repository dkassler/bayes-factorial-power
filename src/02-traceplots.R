sims_single <- readRDS('cache/sims_single.Rds')

library(gridExtra)
library(cowplot)

sims1 <- sims %>%
  mutate(extr = fit %>%
           map('result') %>%
           map(extract, pars = 'lp__', permute = FALSE) %>%
           map(. %>% drop %>% as.data.frame %>% mutate(n = seq(nrow(.))))) %>%
  unnest(extr)

sims2 <- sims1 %>%
  mutate(nArms = nA * nB * nC * nD) %>%
  gather('chain', 'lp__', starts_with('chain')) %>%
  mutate(chain = factor(parse_number(chain))) %>% 
  group_by(nArms, N) %>% 
  mutate(lp__ = (lp__ - mean(lp__)) / sd(lp__))

traceplot_labeller = labeller(
  N = function(x) paste('N =', x),
  nArms = c('16' = '2x2x2x2 = 16 arms',
            '72' = '3x4x3x2 = 72 arms',
            '180' = '4x3x5x3 = 180 arms',
            '250' = '5x5x5x2 = 250 arms')
)

ggplot(data = sims2, aes(x = n, y = lp__, color = chain)) +
  facet_grid(N ~ nArms, labeller = traceplot_labeller) +
  geom_path(alpha = 0.8) +
  scale_color_manual(values = c("#FF7F00", "#CDB5CD", "#473C8B")) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.grid = element_blank()) +
  labs(x = 'Sampler Iterations', y = '(Standardized) Log Posterior',
       color = 'Chain') +
  ggsave('graphs/traceplots.png', height = 5, width = 6.5)

sims1 <- sims %>% 
  mutate(tp = fit %>%
           map('result') %>% 
           map(traceplot, pars = 'lp__'))

plot_grid(plotlist = sims1$tp)
