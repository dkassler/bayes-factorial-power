#' My attempt at recreating Mariel's old power sims code with modern packages in
#' order to fix the bug. Borrowing as much from the old code as possible.
---
title: "Power Sims"
output: html_document
---

# Design setup ------------------------------------------------------------

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


# Simulate an experiment --------------------------------------------------

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
                nChains=1){
  
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
  
  
  
