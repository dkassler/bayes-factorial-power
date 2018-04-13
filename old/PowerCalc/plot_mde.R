# Plot results of power calculations

library(tidyverse)
library(here)

outname <- "output_B"

if (!dir.exists(here("output", outname))) dir.create(here("output", outname))

outdir <- here("data/save", outname)
output <- dir(outdir) %>% 
  map_df(~readRDS(file.path(outdir, .)))

#####################################################################################
# for each design, summarize across simulations. for designs with more than 1/3 of  # 
# simulations returning NA, the summary will be equal to NA. otherwise, the summary #
# will be equal to the median value.                                                #
#####################################################################################
#output$MDE[is.na(output$MDE)] <- Inf
output$MDE[is.na(output$MDE)] <- 1

output <- filter(output, N > 400)

MDEs <- aggregate(MDE~N+nA+nB+nC+nD,           
                  data=output,
                  function(x)
                    return(ifelse(mean(is.na(x) | x == 1)>1/3, 
                                  NA, median(x, na.rm=T))))

########
# plot #
########
p <- ggplot(output, aes(factor(N), MDE))
p + #geom_abline(int=.25, sl=0, col='darkgrey') +
  geom_boxplot(aes(fill=factor(nA, levels=c(2,4,3,5),
                               labels=c('2x2x2x2 = 16 arms',
                                        '2x3x3x4 = 72 arms',
                                        '2x3x3x10 = 180 arms',
                                        '2x5x5x5 = 250 arms')))) +
  guides(fill=guide_legend(title='')) +
  labs(x='Sample Size') +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1), 
                     labels = c('0.00', '0.25', '0.50', '0.75', expression(phantom(x) >= '1.00'))) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(file=here('output', outname, 'distribution.pdf'), width=6.5, height=5)
ggsave(file=here('output', outname, 'distribution.png'), width=6.5, height=5)

ggplot(MDEs, aes(N, MDE)) +
  #geom_abline(int=.25, sl=0, col='darkgrey') +
  geom_line(aes(col=factor(nA, levels=c(2,4,3,5), 
                           labels=c('2x2x2x2 = 16 arms',
                                    '2x3x3x4 = 72 arms',
                                    '2x3x3x10 = 180 arms',
                                    '2x5x5x5 = 250 arms')))) +
  geom_point(aes(col=factor(nA, levels=c(2,4,3,5), 
                            labels=c('2x2x2x2 = 16 arms',
                                     '2x3x3x4 = 72 arms',
                                     '2x3x3x10 = 180 arms',
                                     '2x5x5x5 = 250 arms'))))+
  guides(col=guide_legend(title='')) +
  labs(x='Sample Size') +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1), 
                     labels = c('0.00', '0.25', '0.50', '0.75', expression(phantom(x) >= '1.00')),
                     limits = c(0, 0.75))+
  scale_x_continuous(limits = c(0, 10000))+
  theme_bw()
ggsave(file=here('output', outname, 'median.pdf'), width=6.5, height=5)
ggsave(file=here('output', outname, 'median.png'), width=6.5, height=5)

ggplot(data = MDEs, aes(x = nA * nB * nC * nD, y = MDE, group = N, color = factor(N))) +
  geom_point() + 
  geom_line()+
  labs(x = "Number of Arms", color = "Sample\nSize") +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1), labels = c('0.00', '0.25', '0.50', '0.75', expression(phantom(x) >= '1.00')))
ggsave(file=here('output', outname, 'mdenarm.pdf'), width=6.5, height=5)
ggsave(file=here('output', outname, 'mdenarm.png'), width=6.5, height=5)
