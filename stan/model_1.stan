functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}

data {
 int<lower = 1> N; //Number of patient
 int<lower = 0, upper = 1> Y_Smear[N];
 int<lower = 0, upper = 1> Y_Mgit[N];
 int<lower = 0, upper = 1> Y_Xpert[N];
 int<lower = 0> Y_Img[N];
 real Y_lnLymGlu[N]; 
 matrix[N, 13] X; //Covariates
}

parameters {
  // For probit regression
  real a0; //intercept
  vector[13] a; //slope
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  ordered[2] eta_Img;
  ordered[6] c_Img1;
  ordered[6] c_Img2;
  ordered[2] mu_lnLymGlu;
}

transformed parameters{
  vector<lower=0, upper=1>[N] theta = Phi(a0 + X*a);
  
  ordered[2] p_Smear = Phi(z_Smear); 
  ordered[2] p_Mgit = Phi(z_Mgit);
  ordered[2] p_Xpert = Phi(z_Xpert);
}

model {
  
 //Specificity of each test
 z_Xpert[1] ~ normal(2.886, 7.923);
 z_Mgit[1] ~ normal(3.023, 7.923);
 z_Smear[1] ~ normal(1, 7.923);

// Other class
 z_Xpert[2] ~ normal(0, 1);
 z_Mgit[2] ~ normal(0, 1);
 z_Smear[2] ~ normal(0, 1);
 
 a0 ~ normal(0, 1);
 a ~ normal(0, 1);

  mu_lnLymGlu[1] ~ normal(0, 1);
  mu_lnLymGlu[2] ~ normal(5, 2);
  c_Img1 ~ induced_dirichlet(rep_vector(1, 7), 0);
  c_Img2 ~ induced_dirichlet(rep_vector(1, 7), 0);
 
 // z_theta ~ normal(a0 + X*a, 1);
 
 for (n in 1:N)
   target += log_mix(theta[n],
                     bernoulli_lpmf(Y_Xpert[n] | p_Xpert[1]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[1]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[1]) + ordered_probit_lpmf(Y_Img[n] | eta_Img[1], c_Img1) + normal_lpdf(Y_lnLymGlu[n]|mu_lnLymGlu[1], 5),
                     bernoulli_lpmf(Y_Xpert[n] | p_Xpert[2]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[2]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[2]) + ordered_probit_lpmf(Y_Img[n] | eta_Img[2], c_Img2) + normal_lpdf(Y_lnLymGlu[n]|mu_lnLymGlu[2], 5));
}
