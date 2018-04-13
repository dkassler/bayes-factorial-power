data { 
  int<lower=1> N;
  int<lower=1> nA;
  int<lower=1> nB;
  int<lower=1> nC;
  int<lower=1> nD;
  int<lower=1> nAB;
  int<lower=1> nAC;
  int<lower=1> nBC;
  int<lower=1> nAD;
  int<lower=1> nBD;
  int<lower=1> nCD;
  int<lower=1> nBatches;
  int<lower=1> n[N];
  real yBar[N];
  int<lower = 1, upper = nA> a[N];
  int<lower = 1, upper = nB> b[N];
  int<lower = 1, upper = nC> c[N];
  int<lower = 1, upper = nD> d[N];
  int<lower = 1, upper = nAB> ab[N];
  int<lower = 1, upper = nAC> ac[N];
  int<lower = 1, upper = nBC> bc[N];
  int<lower = 1, upper = nAD> ad[N];
  int<lower = 1, upper = nBD> bd[N];
  int<lower = 1, upper = nCD> cd[N];
}

parameters {
  real mu;
  // real <lower = 0.01, upper = 10> lambda;
  real <lower = 0.01, upper = 10> tau;
  // real <lower = 0.01, upper=pi()/2> sigma_unif[nBatches];
  vector[nA] A_untransformed;
  vector[nB] B_untransformed;
  vector[nC] C_untransformed;
  vector[nD] D_untransformed;
  vector[nAB] AB_untransformed;
  vector[nAC] AC_untransformed;
  vector[nBC] BC_untransformed;
  vector[nAD] AD_untransformed;
  vector[nBD] BD_untransformed;
  vector[nCD] CD_untransformed;
  vector<lower = 0>[nBatches] sigma;
}

transformed parameters {
  vector[nA] A;
  vector[nB] B;
  vector[nC] C;
  vector[nD] D;
  vector[nAB] AB;
  vector[nAC] AC;
  vector[nBC] BC;
  vector[nAD] AD;
  vector[nBD] BD;
  vector[nCD] CD;
  real yHat[N];
  // for(i in 1:nBatches)
  //   sigma[i] <- lambda * tan(sigma_unif[i]); // sigma~cauchy(0,lambda)
  A  <- A_untransformed * sigma[1];
  B  <- B_untransformed * sigma[2];
  C  <- C_untransformed * sigma[3];
  D  <- D_untransformed * sigma[4];
  AB <- AB_untransformed * sigma[5];
  AC <- AC_untransformed * sigma[6];
  BC <- BC_untransformed * sigma[7];
  AD <- AD_untransformed * sigma[8];
  BD <- BD_untransformed * sigma[9];
  CD <- CD_untransformed * sigma[10];
  for(i in 1:N)
    yHat[i] <- mu +
               A[a[i]] +
               B[b[i]] +
               C[c[i]] +
               D[d[i]] +
               AB[ab[i]] +
               AC[ac[i]] +
               BC[bc[i]] +
               AD[ad[i]] +
               BD[bd[i]] +
               CD[cd[i]];
}

model {
  real s[N];
  for(i in 1:N)
    s[i] <- tau/sqrt(n[i]);
  yBar ~ normal(yHat, s);
  A_untransformed ~ normal(0,1);
  B_untransformed ~ normal(0,1);
  C_untransformed ~ normal(0,1);
  D_untransformed ~ normal(0,1);
  AB_untransformed ~ normal(0,1);
  AC_untransformed ~ normal(0,1);
  BC_untransformed ~ normal(0,1);
  AD_untransformed ~ normal(0,1);
  BD_untransformed ~ normal(0,1);
  CD_untransformed ~ normal(0,1);
  sigma ~ normal(0, 1);
}