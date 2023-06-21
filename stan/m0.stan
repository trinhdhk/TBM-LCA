functions{
#include includes/functions.stan
}

data {
  int<lower = 1> N_all; //Number of patient
  int<lower=0, upper=1> unsure_spc;

#include includes/data/Y.stan
#include includes/cross_validation/data.stan
}

transformed data{
  // * Global variables -------------------------------------------------------
#include includes/cross_validation/transform_data_Y.stan
}
 
parameters {
  // Parameters of the logistics regression -----------------------------------
  real a0; //intercept
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
}

model {
  int nu = 6;
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
  
  for (n in 1:N){
    real theta = inv_logit(a0);
    
    target += log_mix(theta,
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
  }
}


generated quantities {
  vector[N_all] log_lik;
  real p_Smear;
  real p_Mgit;
  real p_Xpert;
  real<lower=0, upper=1> theta;

  theta = inv_logit(a0);

  {
    p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta * inv_logit(z_Smear[2]);
    p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta * inv_logit(z_Mgit[2]);
    p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta * inv_logit(z_Xpert[2]);
    
    for (n in 1:N_all) {
      log_lik[n] = log_mix(theta,
      bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[2]),
      bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1]));
    }
  }
}
