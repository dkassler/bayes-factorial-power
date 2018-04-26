# A fit function takes arguments data and coef. Other arguments can be partially
# applied away.

# Bayesian Fit Functions --------------------------------------------------

fit_bayes <- function(data, coef, model, ...) {
  coef$yBar <- aggregate(y~abcd, mean, data=data)$y[coef$abcd]
  coef$n    <- aggregate(y~abcd, length, data=data)$y[coef$abcd]

  stanData <- list(N = nrow(coef),
                   nA = max(coef$a),
                   nB = max(coef$b),
                   nC = max(coef$c),
                   nD = max(coef$d))
  stanData <- c(stanData, with(stanData, list(
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
  ))
  
  fit <- sampling(model, data = stanData, ...)
  return(fit)
}

model <- stan_model_cache('powerSim4var_update')


# Frequentist fit functions -----------------------------------------------

fit_freq_inter <- function(data, coef) {
  data$a <- factor(data$a)
  data$b <- factor(data$b)
  data$c <- factor(data$c)
  data$d <- factor(data$d)
  
  fit <-lm(y ~ (a + b + c + d)^2 - 1, data)
}

fit_freq_main <- function(data, coef) {
  data$a <- factor(data$a)
  data$b <- factor(data$b)
  data$c <- factor(data$c)
  data$d <- factor(data$d)
  
  fit <-lm(y ~ a + b + c + d - 1, data)
}
