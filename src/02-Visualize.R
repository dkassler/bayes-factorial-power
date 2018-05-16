# source(here::here('PowerCalc/loadproject.R'))
# 
# sims <- readRDS(file.path("U:", "Projects", "PCI Paper", "PowerCalc", "save", "mde.1K.Rds"))

sims <- readRDS("../mde_full.Rds")

# Boxplot -----------------------------------------------------------------

boxplot_data <- sims %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>% 
  replace_na(list(MDE = 1))

ggplot(data = boxplot_data) +
  facet_wrap(~sigfun) + 
  geom_boxplot(aes(x = factor(N), fill = nArms, y = MDE))

ggsave('graphs/boxplot.png')

# Line plot median MDE ----------------------------------------------------

lineplot_data <- sims %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>%
  group_by(N, nArms, sigfun, bayes, efftype) %>% 
  summarise(MDE = if (mean(is.na(MDE)) > 1/3) NA else 
    median(replace_na(MDE, Inf)))

lineplot <- function(.x) ggplot(data = lineplot_data %>% filter(efftype == .x), 
       aes(x = N, y = MDE, color = nArms, shape = bayes, linetype = bayes)) +
  geom_point() + 
  geom_line()

lineplot('arms')
ggsave('graphs/lineplot_arms.png')

lineplot('main')
ggsave('graphs/lineplot_main.png')
