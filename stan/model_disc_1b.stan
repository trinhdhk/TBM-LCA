data {
 int<lower = 1> N; //Number of patient
 int<lower = 1> nX;
 int<lower = 0, upper = 1> Y_Smear[N];
 int<lower = 0, upper = 1> Y_Mgit[N];
 int<lower = 0, upper = 1> Y_Xpert[N];
 int<lower = 0, upper = 1> Y_Img[N];
 int<lower = 0, upper = 1> Y_LymGluRatio[N];
 matrix[N, nX] X; //Covariates
}

parameters {
  // For probit regression
  real a0; //intercept
  vector[nX] a; //slope
  
  //Probability of each vars
  ordered[2] p_Smear; 
  ordered[2] p_Mgit;
  ordered[2] p_Xpert;
  ordered[2] p_Img;
  ordered[2] p_LymGluRatio;
}

model {
  real theta[N] = to_array_1d(Phi(a0 + X*a));
  
  //Specificity of each test
  p_Xpert[1] ~ beta_proportion(.001, 1);
  p_Mgit[1] ~ beta_proportion(.001, 1);
  p_Smear[1] ~ beta_proportion(.001, 1);
  
  // Other class
  p_Xpert[2] ~ beta_proportion(.593, 379);
  p_Mgit[2] ~ beta_proportion(.665, 379);
  p_Smear[2] ~ beta_proportion(.786, 379);
    
  //Same for two classes
  p_Img ~ beta(0.5, 1);
  p_LymGluRatio ~ beta(0.5, 1);
 
  //Prior for covariates
  a0 ~ normal(0, 1);
  a ~ normal(0, 1);
  
  for (n in 1:N){
    real x = 0;
    target += log_mix(theta[n],
    bernoulli_lpmf(Y_Xpert[n] | p_Xpert[2]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[2]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[2]) + bernoulli_lpmf(Y_Img[n] | p_Img[2]) + bernoulli_lpmf(Y_LymGluRatio[n] | p_LymGluRatio[2]),
    bernoulli_lpmf(Y_Xpert[n] | p_Xpert[1]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[1]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[1]) + bernoulli_lpmf(Y_Img[n] | p_Img[1]) + bernoulli_lpmf(Y_LymGluRatio[n] | p_LymGluRatio[1]));
  }
}
