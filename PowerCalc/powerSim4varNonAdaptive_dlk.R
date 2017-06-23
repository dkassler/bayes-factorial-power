##############
# get set up #
##############
library(rstan)              # for fitting Bayesian models
library(doParallel)         # for doing parallel computing
library(plyr)               # for data wrangling
library(ggplot2)            # for plotting
library(glm2)               # provides stability for models that may fail to converge using glm()

setwd(here::here("PowerCalc"))
options(error = kasslR::send_error())

###################################################
# which designs are you interested in looking at? #
###################################################
designs <- expand.grid(N=c(1000,2000,5000),
                       nA=c(2,3,4,5),
                       cauchyScale=c(5),
                       diff=c(.25))

designs$nB[designs$nA==2] <- 2   # 2 x 2 x 2 x 2
designs$nC[designs$nA==2] <- 2
designs$nD[designs$nA==2] <- 2

designs$nB[designs$nA==3] <- 3   # 3 x 3 x 10 x 2
designs$nC[designs$nA==3] <- 10
designs$nD[designs$nA==3] <- 2

designs$nB[designs$nA==4] <- 3   # 3 x 3 x 4 x 2
designs$nC[designs$nA==4] <- 3
designs$nD[designs$nA==4] <- 2

designs$nB[designs$nA==5] <- 5   # 5 x 5 x 5 x 2
designs$nC[designs$nA==5] <- 5
designs$nD[designs$nA==5] <- 2

############################################################################
# how many times do you want to simulated an experiment under each design? #
############################################################################
nSims <- 25

#########################################################
# create an object that will hold, for each simulation, # 
# information about the design and the MDE              #
#########################################################
output <- designs
for(i in 2:nSims)
  output <- rbind(output, designs)

#######################################################################
# the next block of code is just for testing. you shouldn't need this #
#######################################################################
{
#   i <- 4
#   N <- output$N[i] 
#   nA <- output$nA[i]
#   nB <- output$nB[i]
#   nC <- output$nC[i]
#   nD <- output$nD[i]
#   cauchyScale <- output$cauchyScale[i]
#   diff <- output$diff[i]
#   sA <- output$sA[i]
#   sB <- output$sB[i]
#   sC <- output$sC[i]
#   sD <- output$sD[i]
#   sAB <- output$sAB[i]
#   sAC <- output$sAC[i]
#   sBC <- output$sBC[i]
#   sAD <- output$sAD[i]
#   sBD <- output$sBD[i]
#   sCD <- output$sCD[i]
#   nIter <- 550 
#   nWarmup=50 
#   nChains=1
}

##################################################################################
# this function simulates an experiment under a given design and returns the MDE #
##################################################################################
sim <- function(N=5000,    # total sample size
                nA=5,      # the number of levels that each variable can take             
                nB=5, 
                nC=5, 
                nD=5, 
                cauchyScale=5,
                diff=.25,   # difference between mean and max outcome
                nWarmup=50, # mcmc settings
                nIter=550,
                nChains=1){             # difference in effectiveness between avg and most effective arms
  
  ########################################################
  # randomly draw variance components from a half cauchy #
  ########################################################
  sVals <- rcauchy(100,0,cauchyScale)
  sVals <- sVals[sVals>0][1:10]
  sA <- sVals[1]
  sB <- sVals[2]
  sC <- sVals[3]
  sD <- sVals[4]
  sAB <- sVals[5]
  sAC <- sVals[6]
  sBC <- sVals[7]
  sAD <- sVals[8]
  sBD <- sVals[9]
  sCD <- sVals[10]
  
  ###########################
  # randomly draw a 'truth' #
  ###########################
  # randomly draw values for each batch of coefficients
  A <- rnorm(nA, 0, sA)
  B <- rnorm(nA, 0, sB)
  C <- rnorm(nC, 0, sC)
  D <- rnorm(nD, 0, sD)  
  AB <- rnorm(nA*nB, 0, sAB)
  AC <- rnorm(nA*nC, 0, sAC)
  BC <- rnorm(nB*nC, 0, sBC)
  AD <- rnorm(nA*nD, 0, sAD)
  BD <- rnorm(nB*nD, 0, sBD)
  CD <- rnorm(nC*nD, 0, sCD)
  
  # create a data frame for calculating the mean effectivenss (mu) of
  # each arm as a sum of coefficient values
  coef <- expand.grid(a=1:nA, b=1:nB, c=1:nC, d=1:nD)
  coef$ab <- as.numeric(as.factor(paste(coef$a, coef$b)))
  coef$ac <- as.numeric(as.factor(paste(coef$a, coef$c)))
  coef$bc <- as.numeric(as.factor(paste(coef$b, coef$c)))
  coef$ad <- as.numeric(as.factor(paste(coef$a, coef$d)))
  coef$bd <- as.numeric(as.factor(paste(coef$b, coef$d)))
  coef$cd <- as.numeric(as.factor(paste(coef$c, coef$d)))
  coef$abcd<- 1:nrow(coef)
  coef$A <- A[coef$a]
  coef$B <- B[coef$b]
  coef$C <- C[coef$c]
  coef$D <- D[coef$d]
  coef$AB <- AB[coef$ab]
  coef$AC <- AC[coef$ac]
  coef$BC <- BC[coef$bc]
  coef$AD <- AD[coef$ad]
  coef$BD <- BD[coef$bd]
  coef$CD <- CD[coef$cd]
  
  coef$mu <- coef$A + coef$B + coef$C + coef$D +
    coef$AB + coef$AC + coef$BC +
    coef$AD + coef$BD + coef$CD
  
  # center and scale the mean effectiveness values (mu) such that the 
  # difference between the average effectineness and the maximum 
  # effectiveness is equal to the 'diff' values
  coef$mu <- coef$mu-mean(coef$mu)
  coef$mu <- coef$mu/max(coef$mu)*diff
  
  #######################################
  # generate random data from the model #
  #######################################
  armOrder <- sample(1:nrow(coef), nrow(coef))
  data <- coef[rep(armOrder, length=N),c('abcd', 'mu')] 
  
  # generate outcomes
  data$y <- rnorm(N, data$mu, 1-.12) # [Note: R2=.12-->12% of variance explained
  #  by the covariates.]
  
  coef$yBar <- aggregate(y~abcd, mean, data=data)$y[coef$abcd]
  coef$n    <- aggregate(y~abcd, length, data=data)$y[coef$abcd]
  
  #################
  # fit the model #
  #################
  # drop empty arms
  stanData <- list(N = nrow(coef),
                   nA = nA,
                   nB = nB,
                   nC = nC,
                   nD = nD,
                   nAB = nA*nB,
                   nAC = nA*nC,
                   nBC = nB*nC,
                   nAD = nA*nD,
                   nBD = nB*nD,
                   nCD = nC*nD,
                   nBatches = 10,
                   n = coef$n[!is.na(coef$yBar)],
                   yBar = coef$yBar[!is.na(coef$yBar)],
                   a = coef$a[!is.na(coef$yBar)],
                   b = coef$b[!is.na(coef$yBar)],
                   c = coef$c[!is.na(coef$yBar)],
                   d = coef$d[!is.na(coef$yBar)],
                   ab = coef$ab[!is.na(coef$yBar)],
                   ac = coef$ac[!is.na(coef$yBar)],
                   bc = coef$bc[!is.na(coef$yBar)],
                   ad = coef$ad[!is.na(coef$yBar)],
                   bd = coef$bd[!is.na(coef$yBar)],
                   cd = coef$cd[!is.na(coef$yBar)])
  
  #mod4var <- stan_model('powerSim4var.stan')
  #save(mod4var, file='mod4var.RData')
  load('mod4var.RData')
  fit <- sampling(mod4var, data = stanData, iter=nIter, warmup=nWarmup, 
                  chains=nChains, refresh=F, pars='yHat')
  
  sims <- extract(fit)    
  
  PPs <- data.frame(mu1=rep(NA,choose(length(coef$mu), 2)), 
                    mu2=rep(NA,choose(length(coef$mu), 2)), 
                    PP=rep(NA, choose(length(coef$mu), 2)))
  index <- 1
  for(i in 1:(length(coef$mu)-1))
    for(j in (i+1):length(coef$mu)){
      PPs$mu1[index] <- coef$mu[i]
      PPs$mu2[index] <- coef$mu[j]
      PPs$PP[index]  <- mean(sims$yHat[,i]>sims$yHat[,j])
      index <- index+1
    }
  
  PPs$diff <- abs(PPs$mu1-PPs$mu2)
  PPs$PPbetter <- PPs$PP
  PPs$PPbetter[PPs$mu1<PPs$mu2] <- 1-PPs$PP[PPs$mu1<PPs$mu2]
  PPs$SigBetter <- PPs$PPbetter>.975
  if(sum(PPs$SigBetter)>0){
    glmFit <- glm2(SigBetter~diff, data=PPs, family=binomial)
    esVals <- seq(0,1,length=10001)
    xBetaHat <- coef(glmFit)[1]+coef(glmFit)[2]*esVals
    pHat <- exp(xBetaHat)/(1+exp(xBetaHat))
    
    #     png('PPs.png', height=5, width=5, units='in', res=400)
#                 plot(PPs$diff, PPs$PPbetter, pch=19, cex=.5, col=grey(.1, alpha=.025),
#                      xlab='Difference between arms',
#                      ylab='Posterior probability of superiority of the better arm',
#                      main='How small a difference can we detect?')
#                 abline(h=0.975,col=2)
#                 abline(v=max(PPs$diff[PPs$PPbetter<.975]), col=2)
#                 lines(esVals, pHat)
#                 abline(h=.5, col=3)
#                 abline(v=esVals[min(which(pHat>.5))], col=3)
#                 abline(h=.8, col=4)
#                 abline(v=esVals[min(which(pHat>.8))], col=4)
#                 box()
    #     dev.off()
    
    return(MDE=esVals[min(which(pHat>.8))])} else 
      (return(NA))
}

##########
# do it! #
##########
cl <- makeCluster(64, outfile = "powerSim4varNonAdaptive_out.txt")
registerDoParallel(cl)
# registerDoSEQ()
startTime <- Sys.time()
output$MDE <- foreach(i = 1:nrow(output), .combine='c',
                     .packages=c('rstan', 'glm2', 'plyr')) %dopar%
# i <- 4
# s <- 
  sim(N=output$N[i],
      nA=output$nA[i],
      nB=output$nB[i],
      nC=output$nC[i],
      nD=output$nD[i],
      cauchyScale=output$cauchyScale[i],
      diff=.25,
      nWarmup=50,
      nIter=550,
      nChains=1)
Sys.time()-startTime

# kasslR::gandalf()

# save workspace in case of future problems
save.image(file = sprintf("powerSim_complete_%s.Rdata", Sys.Date()))
kasslR::send_mail(msg = sprintf("Stan code finished in %s", format(Sys.time()-startTime)))

#####################################################################################
# for each design, summarize across simulations. for designs with more than 1/3 of  # 
# simulations returning NA, the summary will be equal to NA. otherwise, the summary #
# will be equal to the median value.                                                #
#####################################################################################
output$MDE[is.na(output$MDE)] <- Inf
MDEs <- aggregate(MDE~N+nA+nB+nC+nD,           
                  data=output,
                  function(x)
                    return(ifelse(mean(is.na(x))>1/3, 
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
  labs(x='Sample size') 
ggsave(file='powerSim4varNonAdaptive_median.pdf', width=6.5, height=5)

########
# save #
########
write.csv(output, file='powerSim4varNonAdaptive_output.csv', row.names=F)
write.csv(MDEs, file='powerSim4varNonAdaptive_mde.csv', row.names=F)