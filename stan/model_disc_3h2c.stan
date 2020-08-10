data {
  int<lower = 1> N; //Number of patient
  int<lower = 1> nX;
  int<lower = 0, upper = 1> Y_Smear[N];
  int<lower = 0, upper = 1> Y_Mgit[N];
  int<lower = 0, upper = 1> Y_Xpert[N];
  int<lower = 0, upper = 1> Y_Img[N];
  int<lower = 0, upper = 1> Y_CSF[N];
  matrix<lower=0>[N, nX] X; //Covariates
}

parameters {
  // For probit regression
  real a0; //intercept
  vector[nX] a; //slope
  real b_HIV[3]; //adjustment of RE with HIV;
  vector[3] b_CSF;
  vector[N] RE; //random effects without HIV; 
  real<lower=0> b_RE; //coef of random effect
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  ordered[2] z_Img;
  ordered[2] z_CSF;
}

model {
  //Prevalence of TBM
  real theta[N] = to_array_1d(Phi(a0 + X*a));
  
  // z of each test positivitiy with random effect
  matrix[N,2] z_Smear_RE; 
  matrix[N,2] z_Mgit_RE; 
  matrix[N,2] z_Xpert_RE;
  matrix[N,2] z_CSF_pred;
  
  //Probs of each test become positive
  vector[2] p_Img;
  matrix[N,2] p_CSF;
  matrix[N,2] p_Smear;
  matrix[N,2] p_Mgit;
  matrix[N,2] p_Xpert;
  
  // Priors of covariates
  a0 ~ normal(0,1);
  a ~ normal(0,1);
  
  //Random effects covariates
  RE ~ normal(0,1);
  b_HIV ~ normal(0,1);
  b_CSF ~ normal(0,1);
  
  //1-Specificity of each test
  z_Xpert[1] ~ normal(inv_Phi(.005), .7);
  z_Mgit[1] ~ normal(-3.023, .89);//2.378
  z_Smear[1] ~ normal(-3.023, .89);
  
  //Sensitivity of each test
  z_Xpert[2] ~ normal(inv_Phi(.593), .117);
  z_Mgit[2] ~ normal(inv_Phi(.665), .217);
  z_Smear[2] ~ normal(inv_Phi(.786), .405);
  
  //Same prior of sens and spcs for two classes
  z_CSF ~ normal(0,1);
  z_Img ~ normal(0,1);
  
  //Add random effects
  z_Smear_RE = rep_matrix(z_Smear', N);
  z_Mgit_RE = rep_matrix(z_Mgit', N);
  z_Xpert_RE = rep_matrix(z_Xpert', N);
  
  z_Smear_RE[,2] += b_RE*RE + b_HIV[1]*X[,1];
  z_Mgit_RE[,2] += b_RE*RE + b_HIV[2]*X[,1];
  z_Xpert_RE[,2] += b_RE*RE + b_HIV[3]*X[,1];
  
  z_CSF_pred = rep_matrix(z_CSF', N) + rep_matrix(X[,8:10]*b_CSF, 2);
  
  // Phi to convert to the probabilities of tests returning +
  p_Smear = Phi(z_Smear_RE);
  p_Mgit = Phi(z_Mgit_RE);
  p_Xpert = Phi(z_Xpert_RE);
  p_Img = Phi(z_Img);
  p_CSF = Phi(z_CSF_pred);
  
  
  for (n in 1:N){
    real x = 0;
    target += log_mix(theta[n],
    bernoulli_lpmf(Y_Xpert[n] | p_Xpert[n,2]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[n,2]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[n,2]) + bernoulli_lpmf(Y_CSF[n] | p_CSF[n,2]) + bernoulli_lpmf(Y_Img[n] | p_Img[2]),
    bernoulli_lpmf(Y_Xpert[n] | p_Xpert[n,1]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[n,1]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[n,1]) + bernoulli_lpmf(Y_CSF[n] | p_CSF[n,1]) + bernoulli_lpmf(Y_Img[n] | p_Img[1]));
  }
}
