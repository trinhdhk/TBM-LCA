data {
  int<lower=1> N; //Number of patient
  int<lower=1> nX;
  // For probit regression
  real a0; //intercept
  vector[nX] a; //slope
  // real b_age; //adjustment of RE with age;
  real b[3];
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  ordered[2] z_Img;
  positive_ordered[2] lambda_brainImg;
  matrix<lower=0>[N, nX] X; //Covariates
}

parameters {
  vector<lower=0>[N] RE; //base random effect; 
}

model {
  //Random effects covariates
  RE ~ normal(0, .8);
}

generated quantities{
  int Y_Xpert[N]; 
  int Y_Mgit[N];
  int Y_Smear[N];
  // int Y_Img[N];
  int Y_brainImg[N];
  int C[N];
  
  //Prevalence of TBM
  real theta[N] = to_array_1d(Phi(a0 + X*a));
  
  // z of each test positivitiy with random effect
  matrix[N,2] z_Smear_RE; 
  matrix[N,2] z_Mgit_RE; 
  matrix[N,2] z_Xpert_RE;
  
  //Probs of each test become positive
  // vector[2] p_Img;
  // matrix[N,2] p_Smear;
  // matrix[N,2] p_Mgit;
  // matrix[N,2] p_Xpert;
  // 
    //Add random effects
  z_Smear_RE = rep_matrix(z_Smear', N);
  z_Mgit_RE = rep_matrix(z_Mgit', N);
  z_Xpert_RE = rep_matrix(z_Xpert', N);
  
  // z_Smear_RE[,2] += b[1]*(RE + b_age*X[,1]);
  // z_Mgit_RE[,2] += b[2]*(RE + b_age*X[,1]);
  // z_Xpert_RE[,2] += b[3]*(RE + b_age*X[,1]);
  
  z_Smear_RE[,2] += b[1]*RE;
  z_Mgit_RE[,2] += b[2]*RE;
  z_Xpert_RE[,2] += b[3]*RE;
  
  // Phi to convert to the probabilities of tests returning +
  // p_Smear = Phi(z_Smear_RE);
  // p_Mgit = Phi(z_Mgit_RE);
  // p_Xpert = Phi(z_Xpert_RE);
  // p_Img = Phi(z_Img);
  
  // p_Smear = inv_logit(z_Smear_RE);
  // p_Mgit = inv_logit(z_Mgit_RE);
  // p_Xpert = inv_logit(z_Xpert_RE);
  // p_Img = inv_logit(z_Img);
  
  // for (n in 1:N){
  //   Y_Xpert[n] = bernoulli_rng(p_Xpert[n,1]*(1-theta[n]) + p_Xpert[n,2]*theta[n]); 
  //   Y_Mgit[n] = bernoulli_rng(p_Mgit[n,1]*(1-theta[n]) + p_Mgit[n,2]*theta[n]); 
  //   Y_Smear[n] = bernoulli_rng(p_Smear[n,1]*(1-theta[n]) + p_Smear[n,2]*theta[n]); 
  //   Y_Img[n] = bernoulli_rng(p_Img[1]*(1-theta[n]) + p_Img[2]*theta[n]); 
  //   Y_brainImg[n] = poisson_rng(lambda_brainImg[1]*(1-theta[n]) + lambda_brainImg[2]*theta[n]);
  // }
  
  C = bernoulli_rng(theta);
  
  for (n in 1:N){
    if (C[n] == 1){
      Y_Xpert[n] = bernoulli_logit_rng(z_Xpert_RE[n,2]);
      Y_Mgit[n] = bernoulli_logit_rng(z_Mgit_RE[n,2]);
      Y_Smear[n] = bernoulli_logit_rng(z_Smear_RE[n,2]);
      Y_brainImg[n] = poisson_rng(lambda_brainImg[2]);
    } else {
      Y_Xpert[n] = bernoulli_logit_rng(z_Xpert_RE[n,1]);
      Y_Mgit[n] = bernoulli_logit_rng(z_Mgit_RE[n,1]);
      Y_Smear[n] = bernoulli_logit_rng(z_Smear_RE[n,1]);
      Y_brainImg[n] = poisson_rng(lambda_brainImg[1]);
    }
  }
}
