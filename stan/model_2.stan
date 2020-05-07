data {
 int<lower = 1> N; //Number of patient
 int<lower = 1> nX;
 int<lower = 0, upper = 1> Y_Smear[N];
 int<lower = 0, upper = 1> Y_Mgit[N];
 int<lower = 0, upper = 1> Y_Xpert[N];
 int<lower = 0, upper = 1> Y_Img[N];
 // real Y_lnLymGlu[N];
 real<lower=0> Y_lympho[N];
 matrix[N, nX] X; //Covariates
}

parameters {
  // For probit regression
  real a0; //intercept
  vector[nX] a; //slope
  // random effect = sd_RE * z_RE ~ normal(0, sd_RE) where z_RE ~ normal(0,1) - diverged so set sigma as 1
  // real<lower=0> sd_RE; //sd of random effects
  // vector[N] z_RE; //scale of random effects
  vector[N] RE; //random effects
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  ordered[2] z_Img;
  // ordered[2] mu_lnLymGlu;
  // real<lower=0> sigma_lnLymGlu[2];
  ordered[2] mu_lympho;
  real<lower=0> sigma_lympho[2];
}

// transformed parameters{
//   
// }

model {
  //Prevalence of TBM
  real theta[N] = to_array_1d(Phi(a0 + X*a));
  
  //Random effects
  // vector[N] RE;

  // z of each test positivitiy with random effect
  matrix[N,2] z_Smear_RE; 
  matrix[N,2] z_Mgit_RE; 
  matrix[N,2] z_Xpert_RE;
  
  //Probs of each test become positive
  vector[2] p_Img;
  matrix[N,2] p_Smear;
  matrix[N,2] p_Mgit;
  matrix[N,2] p_Xpert;
  
  // Priors of covariates
  a0 ~ normal(0,1);
  a ~ normal(0,1);
  
  //Random effects covariates
  // sd_RE ~ uniform(0,1);
  // z_RE ~ normal(0,1);
  // RE = sd_RE * z_RE;
  RE ~ normal(0,1);
  
  //1-Specificity of each test
  z_Xpert[1] ~ normal(inv_Phi(.005), .7^2);
  z_Mgit[1] ~ normal(-3.023, 2.378);//2.378
  z_Smear[1] ~ normal(-3.023, 2.378);
  
  //Sensitivity of each test
  z_Xpert[2] ~ normal(inv_Phi(.593), .25^2);
  z_Mgit[2] ~ normal(inv_Phi(.665), .1^2);
  z_Smear[2] ~ normal(inv_Phi(.786), .1^2);
  
  z_Smear_RE = rep_matrix(z_Smear', N);
  z_Mgit_RE = rep_matrix(z_Mgit', N);
  z_Xpert_RE = rep_matrix(z_Xpert', N);
  
  z_Smear_RE[,2] += RE;
  z_Mgit_RE[,2] += RE;
  z_Xpert_RE[,2] += RE;
  
  // Phi to convert to the probabilities of tests returning +
  p_Smear = Phi(z_Smear_RE);
  p_Mgit = Phi(z_Mgit_RE);
  p_Xpert = Phi(z_Xpert_RE);
  p_Img = Phi(z_Img);
  
  //Same prior of sens and spcs for two classes
  // mu_lnLymGlu ~ normal(0,1);
  // sigma_lnLymGlu ~ normal(0,1);  
  mu_lympho ~ normal(0,1);
  sigma_lympho ~ normal(0,1);
  z_Img ~ normal(0, 1);
  
  
  for (n in 1:N){
    real x = 0;
    target += log_mix(theta[n],
    bernoulli_lpmf(Y_Xpert[n] | p_Xpert[n,1]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[n,1]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[n,1]) + bernoulli_lpmf(Y_Img[n] | p_Img[1]) + lognormal_lpdf(Y_lympho[n] | mu_lympho[1], sigma_lympho[1]),// normal_lpdf(Y_lnLymGlu[n] | mu_lnLymGlu[1], sigma_lnLymGlu[1]), //+ normal_lpdf(Y_protein[n] | mu_protein[1], sigma_protein[1]),// - normal_lccdf(0 | mu_protein[1], sigma_protein[1]), //normal_lpdf(lnY_lympho[n] | mu_lympho[1], 1), //+ normal_lpdf(Y_lnLymGlu[n]|mu_lnLymGlu[1], 5),
    bernoulli_lpmf(Y_Xpert[n] | p_Xpert[n,2]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[n,2]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[n,2]) + bernoulli_lpmf(Y_Img[n] | p_Img[2]) + lognormal_lpdf(Y_lympho[n] | mu_lympho[1], sigma_lympho[1]));// normal_lpdf(Y_lnLymGlu[n] | mu_lnLymGlu[2], sigma_lnLymGlu[2])); //+ normal_lpdf(Y_protein[n] | mu_protein[2], sigma_protein[2]));// - normal_lccdf(0 | mu_protein[2], sigma_protein[2]));//normal_lpdf(lnY_lympho[n] | mu_lympho[2], 1)); //+normal_lpdf(Y_lnLymGlu[n]|mu_lnLymGlu[2], 5));
  }
}
