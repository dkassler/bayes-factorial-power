functions {
  int idx_main(int i1, int i2) {
    
  }
}

data {
  int nObs;              // number of participants
  int nFact;             // number of factors
  int nLev[nFact];       // number of levels per factor
  int lev[nObs, nFact];  // assignment of factor levels    
}

transformed data {
  int nInter = nFact * (nFact - 1) / 2;
  int maxLev = max(nLev);
  int nInterLev[nInter];
  
  for (i in 1:(nFact-1)) {
    for (j in (i+1):nFact) {
      nInterLev[i*(i-1)/2 + j] = nLev[i] * nLev[j];
    }
  }
}

parameters {
   real main[sum(nLev)];
   real inter[sum(nInterLev)];
}

transformed parameters {
  
}

model {
  
}