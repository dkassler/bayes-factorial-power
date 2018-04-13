#' This script reproduces the frequentist power calculations in a format that
#' more closely matches the bayesian power calculations. It also moves
#' everything to a more modular functional form, so this should be a useful
#' place to refer back to if we ever need to put together a package that we want
#' to see the light of day.

nSims <- 1#000
design_params <- list(
  N = c(1E3, 2E3, 3.5E3, 5E3, 10E3),
  nArms = list(c(2, 2, 2, 2), 
               c(3, 3, 10, 2),
               c(3, 3, 4, 2),
               c(5, 5, 5, 2))
)

# Packages and setup ------------------------------------------------------
library(tidyverse)          # loading before plyr for backwards compatability
library(plyr)               # for data wrangling
library(glm2)               # provides stability for models that may fail to converge using glm()
library(doParallel)

# Define functions --------------------------------------------------------

gen_true_coef <- function(nA, nB, nC, nD, cauchyScale = 5, diff = 0.25) {
  # randomly draw variance components from a half cauchy #
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
  coef$abcd <- 1:nrow(coef)
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
  
  return(coef)
}

gen_study_data <- function(coef, N) {
  # generate random data from the model #
  armOrder <- sample(1:nrow(coef), nrow(coef))
  data <- coef[rep(armOrder, length=N),
               c('a', 'b', 'c', 'd', 'ab', 'ac', 'ad', 'bd', 'cd', 'mu')] 
  
  # generate outcomes
  data$y <- rnorm(N, data$mu, 1-.12) # [Note: R2=.12-->12% of variance explained
  #  by the covariates.]
  
  #coef$yBar <- aggregate(y~abcd, mean, data=data)$y[coef$abcd]
  #coef$n    <- aggregate(y~abcd, length, data=data)$y[coef$abcd]
  
  return(data)
}

vanilla_sig <- function(data, coef){
  data$a <- factor(data$a)
  data$b <- factor(data$b)
  data$c <- factor(data$c)
  data$d <- factor(data$d)
  
  # fit the model #
  fit <-lm(y ~ (a + b + c + d)^2 - 1, data)
  
  # function to get indicators for each arm
  get_arm_ind <- function(a, b, c, d) {
    nm1 <- paste0(c('a', 'b', 'c', 'd'), c(a, b, c, d))
    com <- combn(nm1, 2)
    nm2 <- paste(com[1,], com[2,], sep = ":")
    nm <- c(nm1, nm2)
    ind <- names(coef(fit)) %in% nm
    ind <- as.numeric(ind)
    ind <- setNames(ind, names(coef(fit)))
  }
  get_arm_ind.coef <- function(i) 
    get_arm_ind(coef[i, 'a'], coef[i, 'b'], coef[i, 'c'], coef[i, 'd'])
  
  coef$arm_ind <- pmap(coef[,c('a', 'b', 'c', 'd')], get_arm_ind)
  
  PPs <- data.frame(mu1=rep(NA,choose(length(coef$mu), 2)), 
                    mu2=rep(NA,choose(length(coef$mu), 2)), 
                    # k=NA,
                    pval=rep(NA, choose(length(coef$mu), 2)))
  
  # pt <- proc.time()
  pairwise <- cross_d(data_frame(i = 1:nrow(coef), j = i),
          .filter = function(i, j) {i >= j}) %>%
    mutate(
      mu1 = coef$mu[i],
      mu2 = coef$mu[j],
      diff = abs(mu1 - mu2)
    )
  k <- map2(pairwise$i, pairwise$j, ~ coef$arm_ind[[.x]] - coef$arm_ind[[.y]])
  k <- do.call(rbind, k)

  est <- (drop(k %*% coef(fit)))
  sd <- drop(sqrt(diag(k %*% vcov(fit) %*% t(k))))
  p <- 2 * pnorm(abs(est), sd = sd, lower.tail = FALSE)
  sig <- ifelse(sign(est) == with(pairwise, sign(mu1 - mu2)),
                p.adjust(p, 'BH') < .05, FALSE)
  pairwise[['SigBetter']] <- sig
  
  # print(proc.time() - pt)
  
  return(pairwise)
  
  #loop over pairs of arms
  # index <- 1
  # for(i in 1:(length(coef$mu)-1))
  #   for(j in (i+1):length(coef$mu)){
  #     PPs$mu1[index] <- coef$mu[i]
  #     PPs$mu2[index] <- coef$mu[j]
  #     
  #     k1 <- get_arm_ind.coef(i)
  #     k2 <- get_arm_ind.coef(j)
  #     k <- k1 - k2
  #     
  #     est <- k %*% coef(fit)
  #     sd <- sqrt(k %*% vcov(fit) %*% k)
  #     p <- 2 * pnorm(est, sd = sd, lower.tail = FALSE)
  #     PPs$pval[index] <- p
  #     
  #     index <- index+1
  #   }
  # #p <- summary(fit)$coef[,'Pr(>|t|)']
  # 
  # PPs$diff <- abs(PPs$mu1-PPs$mu2)
  # PPs$SigBetter <- p.adjust(PPs$pval, 'BH') < .05
  # 
  #print(proc.time() - pt)
  
  # return(PPs)
}

freq_sig_nointer <- function(data, coef){
  data$a <- factor(data$a)
  data$b <- factor(data$b)
  data$c <- factor(data$c)
  data$d <- factor(data$d)
  
  # fit the model #
  fit <-lm(y ~ (a + b + c + d) - 1, data)
  
  # function to get indicators for each arm
  get_arm_ind <- function(a, b, c, d) {
    nm1 <- paste0(c('a', 'b', 'c', 'd'), c(a, b, c, d))
    com <- combn(nm1, 2)
    nm2 <- paste(com[1,], com[2,], sep = ":")
    nm <- c(nm1, nm2)
    ind <- names(coef(fit)) %in% nm
    ind <- as.numeric(ind)
    ind <- setNames(ind, names(coef(fit)))
  }
  get_arm_ind.coef <- function(i) 
    get_arm_ind(coef[i, 'a'], coef[i, 'b'], coef[i, 'c'], coef[i, 'd'])
  
  coef$arm_ind <- pmap(coef[,c('a', 'b', 'c', 'd')], get_arm_ind)
  
  PPs <- data.frame(mu1=rep(NA,choose(length(coef$mu), 2)), 
                    mu2=rep(NA,choose(length(coef$mu), 2)), 
                    # k=NA,
                    pval=rep(NA, choose(length(coef$mu), 2)))
  
  # pt <- proc.time()
  pairwise <- cross_d(data_frame(i = 1:5, j = i),
                      .filter = function(i, j) {i >= j}) %>%
    mutate(
      mu1 = coef$mu[i],
      mu2 = coef$mu[j],
      diff = abs(mu1 - mu2)
    )
  k <- map2(pairwise$i, pairwise$j, ~ coef$arm_ind[[.x]] - coef$arm_ind[[.y]])
  k <- do.call(rbind, k)
  
  est <- drop(k %*% coef(fit))
  sd <- drop(sqrt(diag(k %*% vcov(fit) %*% t(k))))
  p <- 2 * pnorm(est, sd = sd, lower.tail = FALSE)
  pairwise[['SigBetter']] <- p.adjust(p, 'BH') < .05
  
  # print(proc.time() - pt)
  
  return(pairwise)
  
  #loop over pairs of arms
  # index <- 1
  # for(i in 1:(length(coef$mu)-1))
  #   for(j in (i+1):length(coef$mu)){
  #     PPs$mu1[index] <- coef$mu[i]
  #     PPs$mu2[index] <- coef$mu[j]
  #     
  #     k1 <- get_arm_ind.coef(i)
  #     k2 <- get_arm_ind.coef(j)
  #     k <- k1 - k2
  #     
  #     est <- k %*% coef(fit)
  #     sd <- sqrt(k %*% vcov(fit) %*% k)
  #     p <- 2 * pnorm(est, sd = sd, lower.tail = FALSE)
  #     PPs$pval[index] <- p
  #     
  #     index <- index+1
  #   }
  # #p <- summary(fit)$coef[,'Pr(>|t|)']
  # 
  # PPs$diff <- abs(PPs$mu1-PPs$mu2)
  # PPs$SigBetter <- p.adjust(PPs$pval, 'BH') < .05
  # 
  #print(proc.time() - pt)
  
  # return(PPs)
}

freq_sig_nofactor <- function(data, coef){
  data$a <- factor(data$a)
  data$b <- factor(data$b)
  data$c <- factor(data$c)
  data$d <- factor(data$d)
  data$arm <- with(data, interaction(a, b, c, d, lex.order = TRUE))
  
  # fit the model #
  fit <-lm(y ~ arm - 1, data)
  
  # function to get indicators for each arm
  get_arm_ind <- function(a, b, c, d) {
    #nm1 <- paste0(c('a', 'b', 'c', 'd'), c(a, b, c, d))
    nm1 <- paste0('arm')
    com <- combn(nm1, 2)
    nm2 <- paste(com[1,], com[2,], sep = ":")
    nm <- c(nm1, nm2)
    ind <- names(coef(fit)) %in% nm
    ind <- as.numeric(ind)
    ind <- setNames(ind, names(coef(fit)))
  }
  get_arm_ind.coef <- function(i) 
    get_arm_ind(coef[i, 'a'], coef[i, 'b'], coef[i, 'c'], coef[i, 'd'])
  
  coef$arm_ind <- pmap(coef[,c('a', 'b', 'c', 'd')], get_arm_ind)
  
  PPs <- data.frame(mu1=rep(NA,choose(length(coef$mu), 2)), 
                    mu2=rep(NA,choose(length(coef$mu), 2)), 
                    # k=NA,
                    pval=rep(NA, choose(length(coef$mu), 2)))
  
  # pt <- proc.time()
  pairwise <- cross_d(data_frame(i = 1:5, j = i),
                      .filter = function(i, j) {i >= j}) %>%
    mutate(
      mu1 = coef$mu[i],
      mu2 = coef$mu[j],
      diff = abs(mu1 - mu2)
    )
  k <- map2(pairwise$i, pairwise$j, ~ coef$arm_ind[[.x]] - coef$arm_ind[[.y]])
  k <- do.call(rbind, k)
  
  est <- drop(k %*% coef(fit))
  sd <- drop(sqrt(diag(k %*% vcov(fit) %*% t(k))))
  p <- 2 * pnorm(est, sd = sd, lower.tail = FALSE)
  pairwise[['SigBetter']] <- p.adjust(p, 'BH') < .05
  
  # print(proc.time() - pt)
  
  return(pairwise)
}

get_mde <- function(PPs) {
  if(sum(PPs$SigBetter)>0){
    glmFit <- glm2(SigBetter~diff, data=PPs, family=binomial)
    esVals <- seq(0,1,length=10001)
    xBetaHat <- coef(glmFit)[1]+coef(glmFit)[2]*esVals
    pHat <- exp(xBetaHat)/(1+exp(xBetaHat))
    
    return(MDE=esVals[min(which(pHat>.8))])} else 
      (return(NA))
}

sim <- function(nA, nB, nC, nD, N, cauchyScale = 5, diff = 0.25){
  coef <- gen_true_coef(nA, nB, nC, nD, cauchyScale = 5, diff = 0.25)
  data <- gen_study_data(coef, N)
  PPs <- vanilla_sig(data, coef)
  MDE <- get_mde(PPs)
  return(MDE)
}


# Run sims ----------------------------------------------------------------

run_sims <- function(design_params, sig_fun) {
designs <- expand.grid(design_params) %>% 
  dplyr::mutate(
    nA = map_dbl(nArms, 1),
    nB = map_dbl(nArms, 2),
    nC = map_dbl(nArms, 3),
    nD = map_dbl(nArms, 4),
    nArms = map_dbl(nArms, prod)
  )

  output <- dplyr::bind_rows(rep(list(designs), nSims))
  
  startTime <- Sys.time()
  pb <- winProgressBar(min = 0, max = nrow(output), label = "initializing...")

  output$MDE <- NA
  # for(i in 1:nrow(output)) {
  #cl <- makeCluster(64, outfile = "powerSim4varNonAdaptive_out.txt")
  #registerDoParallel(cl)
  registerDoSEQ()
  output$MDE <- foreach(i = 1:nrow(output), .combine = c) %do% {
    # output$MDE[i] <- sim(
    #     N=output$N[i],
    #     nA=output$nA[i],
    #     nB=output$nB[i],
    #     nC=output$nC[i],
    #     nD=output$nD[i]
    # )
    coef <- gen_true_coef(
      nA=output$nA[i],
      nB=output$nB[i],
      nC=output$nC[i],
      nD=output$nD[i]
    )
    data <- gen_study_data(coef, N=output$N[i])
    PPs <- sig_fun(data, coef)
    MDE <- get_mde(PPs); MDE
    setWinProgressBar(pb, i, label = format(Sys.time() - startTime))
    return(MDE)
    # output$MDE[i] <- MDE
  }
  
  print(Sys.time()-startTime)
  close(pb)
  kasslR::alert()
  
  return(output)
}

output <- run_sims(design_params, vanilla_sig)

kasslR::gandalf()

# validation --------------------------------------------------------------


# ggplot(output %>% mutate(MDE = plyr::mapvalues(MDE, Inf, 1))) + 
ggplot(output %>% replace_na(list(MDE = 1))) + 
  geom_boxplot(aes(x = factor(N), fill = factor(nArms), y = MDE))

# post-process sims -------------------------------------------------------

output$MDE[is.na(output$MDE)] <- Inf
MDEs <- aggregate(MDE~N+nA+nB+nC,           
                  data=output,
                  function(x)
                    return(ifelse(mean(x==Inf)>1/3, 
                                  NA, median(x, na.rm=T))))

write.csv(output, file='classicalPowerRedux_output.csv', row.names=F)
write.csv(MDEs, file='classicalPowerRedux_mde.csv', row.names=F)
