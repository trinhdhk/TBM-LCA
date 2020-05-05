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
 int<lower = 1> nX;
 int<lower = 0, upper = 1> Y_Smear[N];
 int<lower = 0, upper = 1> Y_Mgit[N];
 int<lower = 0, upper = 1> Y_Xpert[N];
 int<lower = 1, upper = 7> Y_Img[N];
 // real Y_lnLymGlu[N]; 
 // real Y_lympho[N];
 // real<lower=0> Y_protein[N];
 matrix[N, nX] X; //Covariates
}

transformed data{
  // real lnY_lympho[N] = log(Y_lympho);
}

parameters {
  // For probit regression
  real a0; //intercept
  vector[nX] a; //slope
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  positive_ordered[2] lambda_Img;
  // ordered[2] mu_lnLymGlu;
  // ordered[2] mu_lympho;
  // positive_ordered[2] mu_protein;
  // real<lower=0> sigma_protein[2];
}

// transformed parameters{
  // vector<lower=0, upper=1>[N] theta = Phi(a0 + X*a);
  // ordered[2] p_Smear = Phi(z_Smear); 
  // ordered[2] p_Mgit = Phi(z_Mgit);
  // ordered[2] p_Xpert = Phi(z_Xpert);
// }

model {
  vector[2] p_Smear = Phi(z_Smear); 
  vector[2] p_Mgit = Phi(z_Mgit);
  vector[2] p_Xpert = Phi(z_Xpert);
  
  real theta[N] = to_array_1d(Phi(a0 + X*a));
  
 //Specificity of each test
 z_Xpert[1] ~ normal(inv_Phi(.005), .7^2);
 z_Mgit[1] ~ normal(-3.023, 7);
 z_Smear[1] ~ normal(-3.023, 7);

// Other class
 z_Xpert[2] ~ normal(inv_Phi(.593), .25^2);
 z_Mgit[2] ~ normal(inv_Phi(.665), .25^2);
 z_Smear[2] ~ normal(inv_Phi(.786), .25^2);
 
 a0 ~ normal(0, 1);
 a ~ normal(0, 1);
 
 lambda_Img[1] ~ normal(3,3);
 lambda_Img[2] ~ normal(0,3);
 
 // mu_protein[1] ~ normal(0.25, 1);
 // mu_protein[2] ~ normal(2.5, 1);
 // sigma_protein ~ normal(0, 1);
 
 // mu_lympho[1] ~ normal(log(20), 10);
 // mu_lympho[2] ~ normal(log(200), 10);

 // mu_lnLymGlu[1] ~ normal(0, 1);
 // mu_lnLymGlu[2] ~ normal(0, 1);
 // c_Img1 ~ induced_dirichlet(rep_vector(1, 7), 0);
 // c_Img2 ~ induced_dirichlet(rep_vector(1, 7), 0);
 
 // z_theta ~ normal(a0 + X*a, 1);
 
 for (n in 1:N){
   real x = 0;
   // if (Y_protein[n]<0) target += negative_infinity();
   target += log_mix(theta[n],
      bernoulli_lpmf(Y_Xpert[n] | p_Xpert[1]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[1]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[1]) + poisson_lpmf(Y_Img[n] | lambda_Img[1]), //+ normal_lpdf(Y_protein[n] | mu_protein[1], sigma_protein[1]),// - normal_lccdf(0 | mu_protein[1], sigma_protein[1]), //normal_lpdf(lnY_lympho[n] | mu_lympho[1], 1), //+ normal_lpdf(Y_lnLymGlu[n]|mu_lnLymGlu[1], 5),
      bernoulli_lpmf(Y_Xpert[n] | p_Xpert[2]) + bernoulli_lpmf(Y_Mgit[n] | p_Mgit[2]) + bernoulli_lpmf(Y_Smear[n] | p_Smear[2]) + poisson_lpmf(Y_Img[n] | lambda_Img[2])); //+ normal_lpdf(Y_protein[n] | mu_protein[2], sigma_protein[2]));// - normal_lccdf(0 | mu_protein[2], sigma_protein[2]));//normal_lpdf(lnY_lympho[n] | mu_lympho[2], 1)); //+normal_lpdf(Y_lnLymGlu[n]|mu_lnLymGlu[2], 5));
  // if (Y_protein[n] < 0){
  //   target+= negative_infinity();
  // } else {
  //   target+= log_mix(theta[n],-normal_lccdf(0 | mu_protein[1], .2),-normal_lccdf(0 | mu_protein[2], 1.5));
  // }
 }
  
  
}
