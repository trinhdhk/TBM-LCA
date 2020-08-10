data {
 int<lower = 1> N; //Number of patient
 int<lower = 1> nX;
 int<lower = 0, upper = 1> Y_Smear[N];
 int<lower = 0, upper = 1> Y_Mgit[N];
 int<lower = 0, upper = 1> Y_Xpert[N];
 int<lower = 0, upper = 1> Y_Img[N];
 int<lower = 0, upper = 1> Y_CSF[N];
 matrix[N, nX] X; //Covariates
}

parameters {
  // For probit regression
  real a0; //intercept
  vector[nX] a; //slope
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  ordered[2] z_Img;
  ordered[2] z_CSF;
}

model {
  vector[2] p_Smear = Phi(z_Smear); 
  vector[2] p_Mgit = Phi(z_Mgit);
  vector[2] p_Xpert = Phi(z_Xpert);
  vector[2] p_Img = Phi(z_Img);
  vector[2] p_CSF = Phi(z_CSF);
  
  real theta[N] = to_array_1d(Phi(a0 + X*a));
  
  //Specificity of each test
  z_Xpert[1] ~ normal(inv_Phi(.005), .7);
  z_Mgit[1] ~ normal(-3.023, .75);
  z_Smear[1] ~ normal(-3.023, .89);

  // Other class
  z_Xpert[2] ~ normal(inv_Phi(.593), .117);
  z_Mgit[2] ~ normal(inv_Phi(.665), .217);
  z_Smear[2] ~ normal(inv_Phi(.786), .405);
  
 //Same prior for two classes
 z_Img ~ normal(0, 1);
 z_CSF ~ normal(0, 1);
 
 //Prior for covariates
 a0 ~ normal(0, 1);
 a ~ normal(0, 1);
 
 for (n in 1:N){
   real x = 0;
   target += log_mix(theta[n],
      bernoulli_lpmf(Y_Xpert[n] | p_Xpert[2]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[2]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[2]) + bernoulli_lpmf(Y_Img[n] | p_Img[2]) + bernoulli_lpmf(Y_CSF[n] | p_CSF[2]),
      bernoulli_lpmf(Y_Xpert[n] | p_Xpert[1]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[1]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[1]) + bernoulli_lpmf(Y_Img[n] | p_Img[1]) + bernoulli_lpmf(Y_CSF[n] | p_CSF[1]));
 }
}
