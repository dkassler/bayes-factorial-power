setwd("C:/Users/mfinucane/Dropbox/powerSim")

##############
# get set up #
##############
library(plyr)        # for data wrangling
library(ggplot2)     # for plotting
library(arm)         # for stable glm fitting

############################################################################
# how many times do you want to simulated an experiment under each design? #
############################################################################
nSims <- 500

###################################################
# which designs are you interested in looking at? #
###################################################
designs <- expand.grid(N=c(2000, 2500, 3000),
                       nA=c(2,72,180,250),
                       sd=c(.25))

designs$nB[designs$nA==2] <- 2   # 2 x 2 x 4
designs$nC[designs$nA==2] <- 4

designs$nB[designs$nA==72] <- 1   
designs$nC[designs$nA==72] <- 1

designs$nB[designs$nA==180] <- 1   
designs$nC[designs$nA==180] <- 1

designs$nB[designs$nA==250] <- 1   
designs$nC[designs$nA==250] <- 1

#########################################################
# create an object that will hold, for each simulation, # 
# information about the design and the MDE              #
#########################################################
output <- designs
for(i in 2:nSims)
  output <- rbind(output, designs)

############################################################################
# for each simulation, randomly choose a seed for random number generation #
############################################################################
set.seed(07012015)
output$seed <- sample.int(1e9, nrow(output))

#######################################################################
# the next block of code is just for testing. you shouldn't need this #
#######################################################################
{
#     i <- 6000
#     seed <- output$seed[i]
#     N <- output$N[i] 
#     nA <- output$nA[i]
#     nB <- output$nB[i]
#     nC <- output$nC[i]
#     sd <- output$sd[i]
}

##################################################################################
# this function simulates an experiment under a given design and returns the MDE #
##################################################################################
sim <- function(seed=1,
                N=2500,    # total sample size
                nA=2,      # the number of levels that each variable can take             
                nB=2, 
                nC=4, 
                sd=.25){ 
  
  #############################################
  # set the seed for random number generation #
  #############################################
  set.seed(seed)
  
  ###########################
  # randomly draw a 'truth' #
  ###########################
  # create a data frame for calculating the mean effectivenss (mu)
  coef <- data.frame(mu = rnorm(nA*nB*nC, 0, sd),
                     abc = 1:(nA*nB*nC))
  
  # center and scale the mean effectiveness values (mu) such the standard
  # deviation across the values is equal to 'sd', following the footnote on
  # page 19 of the impaq design report
  coef$mu <- coef$mu-mean(coef$mu)
  coef$mu <- coef$mu/sd(coef$mu)*sd 
  
  #######################################
  # generate random data from the model #
  #######################################
  armOrder <- sample(1:nrow(coef), nrow(coef))
  data <- coef[rep(armOrder, length=N),c('abc', 'mu')] 
  
  # generate outcomes
  data$y <- rnorm(N, data$mu, 1-.12) # [Note: R2=.12-->12% of variance explained
  #  by the covariates.]
  
  #################
  # fit the model #
  #################
  fit <-lm(y~factor(abc)-1, data)
  p <- summary(fit)$coef[,'Pr(>|t|)']
  pBH <- p.adjust(p, 'BH')  
  
## compare each mu to the mean of zero
  PPs <- data.frame(mu=coef$mu)
  PPs$diff <- abs(PPs$mu)
  PPs$SigBetter <- pBH < .05
  if(sum(PPs$SigBetter)>0){
    glmFit <- bayesglm(SigBetter~diff, data=PPs, family=binomial)
    esVals <- seq(0,max(PPs$diff),by=.0001)
    xBetaHat <- coef(glmFit)[1]+coef(glmFit)[2]*esVals
    pHat <- exp(xBetaHat)/(1+exp(xBetaHat))
    
    #     png('PPs.png', height=5, width=5, units='in', res=400)
#                 plot(PPs$diff, PPs$SigBetter, pch=19, cex=.5, col=grey(.1, alpha=.025),
#                      xlab='Difference between arms',
#                      ylab='Posterior probability of superiority of the better arm',
#                      main='How small a difference can we detect?',
#                      xlim=range(esVals))
#                 abline(h=0.975,col=2)
#                 abline(v=max(PPs$diff[PPs$PPbetter<.975]), col=2)
#                 lines(esVals, pHat)
#                 abline(h=.5, col=3)
#                 abline(v=esVals[min(which(pHat>.5))], col=3)
#                 abline(h=.8, col=4)
#                 abline(v=esVals[min(which(pHat>.8))], col=4)
#                 box()
    #     dev.off()
    
  # if the estiamted prob (pHat) for the largest difference (max(diff)) is <.8, 
  # then we're just extrapolating. might as well just set the MDE to 1.
#    if(pHat[esVals == round(max(PPs$diff),4)] > .8)
  if(pHat[length(pHat)] > .8)
    return(MDE=esVals[min(which(pHat>.8))]) else
      return(MDE=NA)} else 
      (return(NA))
}

####################
# do it in serial! #
####################
startTime <- Sys.time()
output$MDE <- NA
for(i in 1:nrow(output))
  output$MDE[i] <- sim(seed=output$seed[i],
      N=output$N[i],
      nA=output$nA[i],
      nB=output$nB[i],
      nC=output$nC[i],
      sd=output$sd[i])
Sys.time()-startTime

#####################################################################################
# for each design, summarize across simulations. for designs with more than 1/3 of  # 
# simulations returning NA, the summary will be equal to NA. otherwise, the summary #
# will be equal to the median value.                                                #
#####################################################################################
output$MDE[is.na(output$MDE)] <- Inf
MDEs <- aggregate(MDE~N+nA+nB+nC,           
                  data=output,
                  function(x)
                    return(ifelse(mean(x==Inf)>1/3, 
                                  NA, median(x, na.rm=T))))
########
# plot #
########
# p <- ggplot(output, aes(factor(N), MDE))
# p + geom_abline(int=.25, sl=0, col='darkgrey') +
#   geom_boxplot(aes(fill=factor(nA, levels=c(2,4,3,5), 
#                                labels=c('2x2x2x2 = 16 arms',
#                                         '2x3x3x4 = 72 arms',
#                                         '2x3x3x10 = 180 arms',
#                                         '2x5x5x5 = 250 arms')))) +
#   guides(fill=guide_legend(title='')) +
#   labs(x='Sample size')
# ggsave(file='powerSim4varNonAdaptive_distribution.pdf', width=6.5, height=5)

p2 <- ggplot(MDEs, aes(N, MDE))
p2 + geom_abline(int=.25, sl=0, col='darkgrey') +
  geom_line(aes(col=factor(nA, levels=c(2,72,180,250), 
                           labels=c('16 arms',
                                    '72 arms',
                                    '180 arms',
                                    '250 arms')))) +
  geom_point(aes(col=factor(nA, levels=c(2,72,180,250), 
                            labels=c('16 arms',
                                     '72 arms',
                                     '180 arms',
                                     '250 arms'))))+
  guides(col=guide_legend(title='')) +
  labs(x='Sample size') 
ggsave(file='VanillaAnovaPower.png', width=6.5, height=5)

########
# save #
########
write.csv(output, file='powerSim4varNonAdaptiveCompareToMean_output.csv', row.names=F)
write.csv(MDEs, file='powerSim4varNonAdaptiveCompareToMean_mde.csv', row.names=F)