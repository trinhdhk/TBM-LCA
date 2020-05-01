data {
 int<lower = 1> N; //Number of patient
 int<lower = 0, upper = 1> Y_Cl[N];
 int<lower = 0, upper = 1> Y_Xp[N];
 int<lower = 0, upper = 1> Y_Sm[N];
 
 matrix[N, 6] X; //Covariates
}

parameters {
  // vector[2] mu;
  // real<lower=0> sigma[2]
  // For logistic regression
  real a0; //intercept
  vector[6] a; //slope
  
  //Probability of each vars
  ordered[2] z_Cl; 
  ordered[2] z_Xp;
  ordered[2] z_Sm;
  
  real z_theta;
  row_vector[N] r; //random effect
  real b_r; //random effect coef
  real mu_r;
  
  int C[N];
}

transformed parameters{
  // simplex[C] theta; //probabily of each class
  real<lower=0, upper=1> theta = inv_logit(z_theta);
  
}

model {
  matrix[N, 2] Z_Cl = rep_matrix(z_Cl', N) + [rep_vector(0, N)', b_r * r]';
  matrix[N, 2] Z_Xp = rep_matrix(z_Xp', N) + [rep_vector(0, N)', b_r * r]';
  matrix[N, 2] Z_Sm = rep_matrix(z_Sm', N) + [rep_vector(0, N)', b_r * r]';
  
  matrix[N, 2] P_Cl;
  matrix[N, 2] P_Xp;
  matrix[N, 2] P_Sm;
  
  real P_Cl_true[N];
  real P_Xp_true[N];
  real P_Sm_true[N];
  
  
  
  z_theta ~ normal(a0 + X*a, 5);
  C ~ bernoulli(theta);
  
  for (n in 1:N){
    for (c in 1:2){
      P_Cl[n,c] = Phi(Z_Cl[n,c]);
      P_Xp[n,c] = Phi(Z_Xp[n,c]);
      P_Sm[n,c] = Phi(Z_Sm[n,c]);
    }
  }
  
 //Specificity of each test
 z_Xp[1] ~ normal(2.886, 5);
 z_Cl[1] ~ normal(3.023, 5);
 z_Sm[1] ~ normal(1.5, 2);

// Other class
 z_Xp[2] ~ normal(0, 1);
 z_Cl[2] ~ normal(0, 1);
 z_Sm[2] ~ normal(0, 1);
 
 //random effect
 mu_r ~ normal(0, 5);
 r ~ normal(mu_r, 5);
 b_r ~ uniform(0, 5);
 
 // theta ~ dirichlet(rep_vector(2.0, C));
 a0 ~ normal(0, 1);
 a ~ normal(0, .1);
 
 for (n in 1:N){
   P_Xp_true[n] = ((1-P_Xp[n,1])^(1-C[n])) * (P_Xp[n,2]^C[n]); 
   P_Cl_true[n] = ((1-P_Cl[n,1])^(1-C[n])) * (P_Cl[n,2]^C[n]); 
   P_Sm_true[n] = ((1-P_Sm[n,1])^(1-C[n])) * (P_Sm[n,2]^C[n]); 
 }
 
 Y_Xp ~ bernoulli(P_Xp_true);
 Y_Cl ~ bernoulli(P_Cl_true);
 Y_Sm ~ bernoulli(P_Sm_true);
 
}
