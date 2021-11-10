  // CSF PCA
  vector[m_csf] L_t0_csf;   // lower diagonal elements of L
  vector<lower=0>[D_csf] L_d_csf;   // lower diagonal elements of L
  real<lower=0> psi0_csf;         // vector of variances
  real<lower=0>   mu_psi_csf;
  real<lower=0>  sigma_psi_csf;
  real mu_lt_csf;
  real<lower=0>  sigma_lt_csf; 
  vector[D_csf] U0_csf[N];
