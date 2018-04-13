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
  B <- rnorm(nB, 0, sB)
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
  # coef$mu <- coef$mu-mean(coef$mu)
  # coef$mu <- coef$mu/max(coef$mu)*diff
  sel <- c('A', 'B', 'C', 'D', 'AB', 'AC', 'AD', 'BC', 'BD', 'CD', 'mu')
  coef[sel] <- sweep(coef[sel], 2, colMeans(coef[sel]), '-')
  s <- diff / max(coef$mu)
  coef[sel] <- sweep(coef[sel], 2, s, '*')
  
  return(coef)
}

gen_study_data <- function(coef, N) {
  stopifnot(N >= nrow(coef))
  # generate random data from the model #
  armOrder <- sample(1:nrow(coef), nrow(coef))
  data <- coef[rep(armOrder, length=N),
               c('a', 'b', 'c', 'd', 'ab', 'ac', 'ad', 'bd', 'cd', 'abcd', 'mu')] 
  
  # generate outcomes
  data$y <- rnorm(N, data$mu, 1-.12) # [Note: R2=.12-->12% of variance explained
  #  by the covariates.]
  
  # coef$yBar <- aggregate(y~abcd, mean, data=data)$y[coef$abcd]
  # coef$n    <- aggregate(y~abcd, length, data=data)$y[coef$abcd]
  
  return(data)
}

get_mde <- function(PPs) {
  if(sum(PPs$SigBetter)>0) {
    glmFit <- glm2::glm2(SigBetter~diff, data=PPs, family=binomial)
    esVals <- seq(0,1,length=10001)
    xBetaHat <- coef(glmFit)[1]+coef(glmFit)[2]*esVals
    pHat <- exp(xBetaHat)/(1+exp(xBetaHat))
    return(MDE=esVals[min(which(pHat>.8))])
  } else { 
      return(NA)
  }
}

simMDE <- function(nA, nB, nC, nD, N, fitfun, sigfun, cauchyScale = 5, diff = 0.25){
  coef <- gen_true_coef(nA, nB, nC, nD, cauchyScale, diff)
  data <- gen_study_data(coef, N)
  fit <- fitfun()
  PPs <- sigfun(data, coef)
  MDE <- get_mde(PPs)
  return(MDE)
}