# Specification for significance functions:
#' @param data a data frame
#' @param coef data frame of coefs
#' 
#' @return a data frame with one row per pair of arms, columns `SigBetter` and `diff`
#' 


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

  return(pairwise)
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

  return(pairwise)
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
  
  return(pairwise)
}



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

sigfun.bayes_arms <- function(data, coef) {
  fit <- fit_bayes(data, coef, model, chains = 1,
                   control = list(adapt_delta = 0.99, max_treedepth = 20))
  sims <- extract(fit, pars = 'yHat')    
  
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
  
  return(PPs)
}
sig.bayes_arms <- function(coef, fit) {
  sims <- extract(fit, pars = 'yHat')    
  
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
  
  return(PPs)
}

sigfun.bayes_maineff <- function(data, coef) {
  fit <- fit_bayes(data, coef, model, chains = 1, 
                   control = list(adapt_delta = 0.99, max_treedepth = 20))
  levs <- coef[1:10]
  sims <- extract(fit, toupper(names(levs)))
  
  f <- function(a, x, i) {
    z1 <- levs[levs[[a]] == x, grep(a, names(levs))]
    z2 <- map2(toupper(names(z1)), z1, ~ sims[[.x]][i, .y])
    z3 <- do.call(cbind, z2)
    z4 <- colMeans(z3)
    sum(z4)
  }
  g <- function(a, x) {
    z <- coef[coef[[tolower(a)]] == x, 
                   grep(toupper(a), names(coef))]
    sum(colMeans(z))
  }
  
  levpairs <- map_dfr(names(levs)[1:4], function(x) {
    n <- n_distinct(levs[[x]])
    com <- combn(1:n, 2)
    data_frame(a = x, x1 = com[1,], x2 = com[2,])
  })
  
  ni = nrow(as.data.frame(fit))
  PPs <- levpairs %>% group_by_all() %>% do(
    PP = {
      mean(sapply(1:ni, function(i) {
        f(.[[1, 'a']], .[[1, 'x1']], i) > f(.[[1, 'a']], .[[1, 'x2']], i)
      }))
    },
    diff = {g(.[[1, 'a']], .[[1, 'x1']]) - g(.[[1, 'a']], .[[1, 'x2']])}
  ) %>% 
    mutate_at(4:5, unlist) %>%
    mutate(PP = ifelse(diff < 0, 1 - PP, PP),
           diff = abs(diff),
           SigBetter = PP > .975) 
  return(PPs)
}
sig.bayes_main <- function(coef, fit) {
  levs <- coef[1:10]
  sims <- extract(fit, toupper(names(levs)))
  
  f <- function(a, x, i) {
    z1 <- levs[levs[[a]] == x, grep(a, names(levs))]
    z2 <- map2(toupper(names(z1)), z1, ~ sims[[.x]][i, .y])
    z3 <- do.call(cbind, z2)
    z4 <- colMeans(z3)
    sum(z4)
  }
  g <- function(a, x) {
    z <- coef[coef[[tolower(a)]] == x, 
                   grep(toupper(a), names(coef))]
    sum(colMeans(z))
  }
  
  levpairs <- map_dfr(names(levs)[1:4], function(x) {
    n <- n_distinct(levs[[x]])
    com <- combn(1:n, 2)
    data_frame(a = x, x1 = com[1,], x2 = com[2,])
  })
  
  ni = nrow(as.data.frame(fit))
  PPs <- levpairs %>% group_by_all() %>% do(
    PP = {
      mean(sapply(1:ni, function(i) {
        f(.[[1, 'a']], .[[1, 'x1']], i) > f(.[[1, 'a']], .[[1, 'x2']], i)
      }))
    },
    diff = {g(.[[1, 'a']], .[[1, 'x1']]) - g(.[[1, 'a']], .[[1, 'x2']])}
  ) %>% 
    mutate_at(4:5, unlist) %>%
    mutate(PP = ifelse(diff < 0, 1 - PP, PP),
           diff = abs(diff),
           SigBetter = PP > .975) 
  return(PPs)
}
