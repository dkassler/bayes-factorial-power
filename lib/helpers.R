design_frame <- function(nA, nB, nC, nD) {
  data_frame(nA = nA, nB = nB, nC = nC, nD = nD)
}

fit_sig_type <- function(bayes, efftype, fitfun, sigfun) {
  data_frame(bayes = bayes, efftype = efftype, fitfun = fitfun, sigfun = sigfun)
}
