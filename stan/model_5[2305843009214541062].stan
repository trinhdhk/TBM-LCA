data {
 int<lower = 1> N; //Number of patient
 int<lower = 0, upper = 1> Y_Cl[N];
 int<lower = 0, upper = 1> Y_Xp[N];
 int<lower = 0, upper = 1> Y_Sm[N];
 matrix[N, 6] X; //Covariates
 int<lower=1> C; //Number of latent class
}

parameters {
  // vector[2] mu;
  // real<lower=0> sigma[2]
  // For logistic regression
  real a0; //intercept
  vector[6] a; //slope
  
  //Probability of each vars
  ordered[C] z_Cl; 
  ordered[C] z_Xp;
  ordered[C] z_Sm;
  
  real z_theta;
}

transformed parameters{
  // simplex[C] theta; //probabily of each class
  real<lower=0, upper=1> theta = Phi(z_theta);
  
  ordered[C] p_Cl = Phi(z_Cl); 
  ordered[C] p_Xp = Phi(z_Xp);
  ordered[C] p_Sm = Phi(z_Sm);
}

model {
 //Specificity of each test
 z_Xp[1] ~ normal(2.886, 7.923);
 z_Cl[1] ~ normal(3.023, 7.923);
 z_Sm[1] ~ normal(1, 7.923);

// Other class
 z_Xp[2] ~ normal(0, 1);
 z_Cl[2] ~ normal(0, 1);
 z_Sm[2] ~ normal(0, 1);
 
 a0 ~ normal(0, 1);
 a ~ normal(0, 1);
 
 z_theta ~ normal(a0 + X*a, 1);
 
 for (n in 1:N)
   target += log_mix(theta,
                     bernoulli_lpmf(Y_Xp[n] | p_Xp[1]) + bernoulli_lpmf(Y_Cl[n] | p_Cl[1]) + bernoulli_lpmf(Y_Sm[n] | p_Sm[1]),
                     bernoulli_lpmf(Y_Xp[n] | p_Xp[2]) + bernoulli_lpmf(Y_Cl[n] | p_Cl[2]) + bernoulli_lpmf(Y_Sm[n] | p_Sm[2]));
}
