data {
  int<lower = 1> N; //Number of patient
  int<lower = 1> nX;
  int<lower = 0, upper = 1> Y_Smear[N];
  int<lower = 0, upper = 1> Y_Mgit[N];
  int<lower = 0, upper = 1> Y_Xpert[N];
  matrix[N, nX] X; //Covariates
}

parameters {
  // For logit regression
  real a0; //intercept
  vector[nX] a; //slope
  real b_HIV; //adjustment of RE with HIV X[,1]
  real b_age; //adjustment of RE with age X[,2];
  // real b_mil; //adjustment of RE with age X[,nX];
  real<lower=0> b[3];
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  
  vector[N] RE; //base random effect;
}

transformed parameters{
  //Prevalence of TBM
  real theta[N] = to_array_1d(inv_logit(a0 + X*a));
  
  // z of each test positivitiy with random effect
  matrix[N,2] z_Smear_RE; 
  matrix[N,2] z_Mgit_RE; 
  matrix[N,2] z_Xpert_RE;
  
  //Bacterial load
  vector[N] bac_load;
  bac_load = RE + b_HIV*X[,1] + b_age*X[,2]; //+ b_mil*X[,nX];
  
  //Add random effects
  z_Smear_RE = rep_matrix(z_Smear' , N);
  z_Mgit_RE  = rep_matrix(z_Mgit' , N);
  z_Xpert_RE = rep_matrix(z_Xpert', N);
  
  z_Smear_RE[,2] += b[1]*bac_load;
  z_Mgit_RE[,2]  += b[2]*bac_load;
  z_Xpert_RE[,2] += b[3]*bac_load;
}

model {
  
  //Probs of each test become positive
  // Priors of covariates
  a0 ~ student_t(5,0,10);
  a ~ student_t(5,0,2.5);
  
  //Random effects covariates
  RE ~ normal(0,1);
  b_age ~ student_t(5,0,2.5);
  b_HIV ~ student_t(5,0,2.5);
  b ~ student_t(5,0,1);
  
  //1-Specificity of each test
  z_Xpert[1] ~ normal(inv_Phi(.005), .7);
  z_Mgit[1] ~ normal(-3.023, .89);//2.378
  z_Smear[1] ~ normal(-3.023, .89);
  
  //Sensitivity of each test
  z_Xpert[2] ~ normal(inv_Phi(.593), .117);
  z_Mgit[2]  ~ normal(inv_Phi(.665), .217);
  z_Smear[2] ~ normal(inv_Phi(.786), .405);
  
  for (n in 1:N){
    target += log_mix(theta[n],
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,2]),// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[2]),
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,1]));// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[1]));
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = log_mix(theta[n],
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,2]),// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[2]),
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,1]));// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[1]));
  }
}
