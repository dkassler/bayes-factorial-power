source(here::here('PowerCalc/loadproject.R'))

sims <- readRDS(file.path("U:", "Projects", "PCI Paper", "PowerCalc", "save", "mde.1K.Rds"))

# Boxplot -----------------------------------------------------------------

boxplot_data <- sims %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>% 
  replace_na(list(MDE = 1))

ggplot(data = boxplot_data) +
  facet_wrap(~sigfun) + 
  geom_boxplot(aes(x = factor(N), fill = nArms, y = MDE))

ggsave('out/boxplot.png')

# Line plot median MDE ----------------------------------------------------

lineplot_data <- sims %>% 
  mutate(nArms = factor(nA * nB * nC * nD)) %>%
  group_by(N, nArms, sigfun) %>% 
  summarise(MDE = if (mean(is.na(MDE)) > 1/3) NA else 
    median(replace_na(MDE, Inf)))

ggplot(data = lineplot_data, aes(x = N, y = MDE, color = nArms)) +
  facet_wrap(~sigfun) + 
  geom_point() + 
  geom_line()

ggsave('out/lineplot.png')
