functions{
#include ./includes/functions.stan
}

data {
  int<lower = 1> N_all; //Number of patient
  int<lower=0, upper=1> Y_Smear_all[N_all];
  int<lower=0, upper=1> Y_Mgit_all[N_all];
  int<lower=0, upper=1> Y_Xpert_all[N_all];
  
  // Hold-out for cross-validation
  int<lower=0, upper=1> keptin[N_all];
}

transformed data{
  int timestamp = 41461735; 
  //this is just a time stamp to force Stan to recompile the code and not used.  
  
  // * Global variables -------------------------------------------------------
#include includes/cross_validation/transform_data_Y.stan
}
 
parameters {
  // Parameters of the logistics regression -----------------------------------
  
  real a0; //intercept
  real<lower=0> b_RE;
  vector[N] RE; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
}

model {
  int nu = 5;
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/m_RE.stan
  
  for (n in 1:N){
    real theta = inv_logit(a0);
    
    real z_Smear_RE = z_Smear[2] + b_RE*RE[n];
    real z_Mgit_RE  = z_Mgit [2] + b_RE*RE[n];
    real z_Xpert_RE = z_Xpert[2] + b_RE*RE[n];
    
    target += log_mix(theta,
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
  }
}


generated quantities {
  vector[N_all] log_lik;
  vector[N_all] p_Smear;
  vector[N_all] p_Mgit;
  vector[N_all] p_Xpert;
  real<lower=0, upper=1> theta;

  theta = inv_logit(a0);

  {
    vector[N_all] RE_all;
    RE_all[which(keptin)] = RE;
    for (n in which_not(keptin)) RE_all[n] = normal_rng(0, 1);
    p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta * inv_logit(z_Smear[2] + b_RE*RE_all);
    p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta * inv_logit(z_Mgit[2]  + b_RE*RE_all);
    p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta * inv_logit(z_Xpert[2] + b_RE*RE_all);
    
    for (n in 1:N_all) {
      log_lik[n] = log_mix(theta,
      bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[2] + b_RE*RE_all[n]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[2]  + b_RE*RE_all[n]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[2] + b_RE*RE_all[n]),
      bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1]));
    }
  }
}
