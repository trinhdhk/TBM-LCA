data {
  int<lower=1> N_pt; // number of observations
  // For logistic regression
  // logit(C) = a0 + X*a
  int<lower=0> N_cov;
  matrix[N_pt, N_cov] X; //covariate design matrix
  
  // For latent class model
  // Currently have no covariate for manifest variables
  // For binary manifest var (mv): P(U[i] == 1|c) = b[c,i] + Z[i]e[i] 
  // where Z[i] is cov design matrix for indicator i and e[i] is param vector for indicator i.
  // For continuous mv: mean[U[i]] = b[c,i] + Z[i]e[i] and U[i] ~ normal(mean[U[i]], s) 
  int<lower=0, upper=1> U[N_pt, 3]; //indicator design matrix
}

parameters {
  
  
  // For logistic regression
  real a0; //intercept
  vector[N_cov] a; //slope
  
  // For latent model
  matrix[2, 3] b;
  real b_re; //random effect param
  vector[N_pt] re; //random effect value
}

transformed parameters{
  // vector<lower=0, upper=1>[N_pt] Pc = exp(logitPc)./(1 + exp(logitPc)); //logit of prob latent classes
  vector<lower=0, upper=1>[N_pt] Pc = rep_vector(1, N_pt);
  vector[N_pt] logitPc;
  matrix[N_pt, 3] P; //whole probability of each test
  matrix[N_pt, 3] spc; //marginal probablity of each test=0 given C=0
  matrix[N_pt, 3] sen; //marginal probablity of each test=1 given C=1
  
  for (i in 1:3){
    Pc = Pc .* (to_vector(U[:, i]) .* P[:, i] + ((1 - to_vector(U[:,i]) .* (1 - P[:,i]))));
  }
  
  logitPc = logit(Pc); //logit of prob latent classes
  
  for (i in 1:3){
    spc[:,i] = Phi(rep_vector(b[1,i], N_pt));
    sen[:,i] = Phi(rep_vector(b[2,i], N_pt) + b_re * re);
  }
  
  for (i in 1:3){
    P[:, i] = (1 - spc[:, i]).*(1 - Pc) + (sen[:, i]).*(Pc);
  }
}

model {
  // Latent class model sampling
  // C = 1
  b[2,:] ~ normal(0,1); 
  b_re ~ uniform(0,5); //random effect param
  re ~ normal(0,1); //random effect value
  
  // C = 0
  b[1,1] ~ normal(2.886, 3.923);//to be change. for xpert.
  b[1,2] ~ normal(3.023, 3.923);//to be change. for culture.
  b[1,3] ~ normal(2, 3.923);
  
  
  //Expected prob of each indicator 
  for (i in 1:3){
    U[:, i] ~ bernoulli(P[:, i]);
  }
  
  // Logistic regression sampling
  a0 ~ normal(0,1);
  a[1:N_cov-1] ~ normal(0,2);
  a[N_cov] ~ skew_normal(2.886, 3.923, -2); //for xray
  
  logitPc ~ normal(a0 + X*a, 1);
}
