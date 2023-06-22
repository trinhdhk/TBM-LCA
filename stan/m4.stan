functions{
#include includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=0> nB; //Number of added RE  
  int<lower=1> B[nB];
  int<lower=0, upper=1> unsure_spc;
  int<lower=0, upper=1> quad_RE;
  int<lower=0, upper=1> obs_Smear_all[N_all];
  int<lower=0, upper=1> obs_Mgit_all [N_all];
  int<lower=0, upper=1> obs_Xpert_all[N_all];
  
#include includes/data/a.stan 
#include includes/data/X.stan
#include includes/data/Y.stan
#include includes/cross_validation/data.stan
#include includes/data/penalty.stan  
#include includes/impute_model/data.stan
}

transformed data{
  // Penalty term adaptation
  int adapt_penalty[2];
  
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariates 
  int nA = nX + nQ; // Number of coef
#include includes/cross_validation/transform_data_Y.stan
#include includes/cross_validation/transform_data_X.stan
#include includes/impute_model/transform_data.stan

  for (i in 1:2) adapt_penalty[i] = penalty_term[i] == 0 ? 1 : 0;
}
 
parameters {
    // Parameters of the imputation model ---------------------------------------
#include includes/impute_model/parameters.stan
  
  // Parameters of the logistics regression -----------------------------------
#include includes/parameters/a.stan
  real<lower=0> b_HIV_raw; //adjustment of RE with HIV Xd[,1]
  vector[nB] b_raw;
  vector<lower=0>[3] b_RE_raw;
  vector<lower=0>[3] b_FE_raw;
  vector[N] RE[2]; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
 
#include includes/parameters/penalty.stan
}


transformed parameters {
#include includes/transform_parameters/a_variable_declaration.stan
  real b_HIV;
  vector[nB] b;
  vector<lower=0>[3] b_RE;
  vector<lower=0>[3] b_FE;
#include includes/impute_model/transform_parameters.stan
  
  {
#include includes/transform_parameters/penalty.stan
#include includes/transform_parameters/a_transform.stan
    
    b_FE = b_FE_raw * SP[2];
    b_RE = b_RE_raw * SP[2];
    b_HIV = b_HIV_raw * SP[2] / mean(b_FE);
    b = b_raw * SP[2] / mean(b_FE);
  }
}


model {
  int nu = 6;
   // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/a.stan
#include includes/main_prior/penalty.stan
  for (j in 1:2) RE[j] ~ normal(0, 1);
  if (nB > 0) {
     int j = 1;
     for (i in B){
       b_raw[j]  ~ normal(0, inv(sd_X[i]));
       j += 1;
     }
  }
  if (penalty_family == 0){
    b_HIV_raw ~ student_t(nu, 0, 2);
    b_RE_raw  ~ student_t(nu, 0, 1);
    b_FE_raw  ~ student_t(nu, 0, 1);
  }
  if (penalty_family == 1){
    b_HIV_raw  ~ double_exponential(0, 2);
    b_RE_raw  ~ double_exponential(0, 1);
    b_FE_raw  ~ double_exponential(0, 1);
  }
  if (penalty_family == 2){
    b_HIV_raw  ~ normal(0, 2);
    b_RE_raw  ~ normal(0, 1);
    b_FE_raw  ~ normal(0, 1);
  }
  for (n in 1:N){
    // If HIV is observed
    if (obs_Xd[n,1]==1){
#include includes/impute_model/impute_priors.observedHIV.stan
      int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
      real bac_load   = b_HIV*Xd_imp[n,1] + dot_product(b, Xc_imp[n,B]) + RE[1,n];
      real z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE[2,n] + square(RE[2,n])*quad_RE);
      real z_Mgit_RE  = z_Mgit[2]  + b_FE[2]*bac_load + b_RE[2]*(RE[2,n] + square(RE[2,n])*quad_RE);
      real z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE[2,n] + square(RE[2,n])*quad_RE);
      
      // Symptoms or motor palsy is missing
      if (N_Xd_miss > 0){
        int N_pattern = int_power(2, N_Xd_miss);
        vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,2:3], obs_Xd[n, 2:3], a[2:3]);
        vector[N_pattern] log_liks;
        pat_thetas[2] += a0 + a[1]*Xd_imp[n,1] + dot_product(a[4:], X_compl[n]);
        
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
          
          log_liks[i] = logprob_theta + log_mix(theta, 
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
        target += log_sum_exp(log_liks);
        // The normal way
      } else {
        row_vector[nA] X = append_col(Xd_imp[n,:], X_compl[n,:]);
        real theta = inv_logit(a0 + dot_product(a, X));
        
        target += log_mix(theta, 
        bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
        bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
      }
    } else { // HIV is missing, the situation got worse
#include includes/impute_model/impute_priors.unobservedHIV.stan
#include includes/impute_model/impute_priors.unobsHIVpos.stan
#include includes/impute_model/impute_priors.unobsHIVneg.stan
  
      // if HIV is positive
      {
        int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
        
        real bac_load   = b_HIV + dot_product(b, Xc_imp[n,B]) + RE[1,n];
        real z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE[2,n] + square(RE[2,n])*quad_RE);
        real z_Mgit_RE  = z_Mgit[2]  + b_FE[2]*bac_load + b_RE[2]*(RE[2,n] + square(RE[2,n])*quad_RE);
        real z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE[2,n] + square(RE[2,n])*quad_RE);
        // Symptoms or motor palsy is missing
        if (N_Xd_miss > 0){
          int N_pattern = int_power(2, N_Xd_miss);
          row_vector[2] csmp_imp = [obs_Xd[n,2] ? Xd[n,2] : theta_cs_hiv2[1], obs_Xd[n,3] ? Xd[n,3] : theta_mp_hiv2[1]];
          vector[N_pattern] pat_thetas[2] = get_patterns(csmp_imp, obs_Xd[n,2:3], a[2:3]);
          vector[N_pattern] log_liks;
          pat_thetas[2] += a0 + a[1] + dot_product(a[4:], X_compl[n]);
          
          for (i in 1:N_pattern){
            real logprob_theta = pat_thetas[1,i];
            real theta = inv_logit(pat_thetas[2,i]);
            
            log_liks[i] = logprob_theta + log_mix(theta, 
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          }
          ll_HIV[1] = log_sum_exp(log_liks);
        } else {
          row_vector[nA-1] X = append_col(Xd_imp[n,2:3], X_compl[n,:]);
          real theta = inv_logit(a0 + a[1] + dot_product(a[2:nA], X));
          
          ll_HIV[1] = log_mix(theta, 
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      }
      
      // If HIV is negative
      {
        int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
        real bac_load   = dot_product(b, Xc_imp[n,B]) + RE[1,n];
        real z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE[2,n] + square(RE[2,n])*quad_RE);
        real z_Mgit_RE  = z_Mgit[2]  + b_FE[2]*bac_load + b_RE[2]*(RE[2,n] + square(RE[2,n])*quad_RE);
        real z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE[2,n] + square(RE[2,n])*quad_RE);
            
        // Symptoms or motor palsy is missing
        if (N_Xd_miss > 0){
          int N_pattern = int_power(2, N_Xd_miss);
          row_vector[2] csmp_imp = [obs_Xd[n,2] ? Xd[n,2] : theta_cs_hiv2[2], obs_Xd[n,3] ? Xd[n,3] : theta_mp_hiv2[2]];
          vector[N_pattern] pat_thetas[2] = get_patterns(csmp_imp, obs_Xd[n, 2:3], a[2:3]);
          vector[N_pattern] log_liks;
          pat_thetas[2] += a0 + dot_product(a[4:], X_compl[n]);
          
          for (i in 1:N_pattern){
            real logprob_theta = pat_thetas[1,i];
            real theta = inv_logit(pat_thetas[2,i]);
            
            log_liks[i] = logprob_theta + log_mix(theta, 
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          }
          ll_HIV[2] = log_sum_exp(log_liks);
        } else {
          row_vector[nA-1] X = append_col(Xd_imp[n,2:3], X_compl[n,:]);
          real theta = inv_logit(a0 + dot_product(a[2:nA], X));
          
          ll_HIV[2] = log_mix(theta, 
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      }
#include includes/impute_model/target.stan
    }
  }
}

generated quantities {
  vector[N_all] log_lik;
  vector[N_all] p_Smear;
  vector[N_all] p_Mgit;
  vector[N_all] p_Xpert;
  vector[3] pairwise_corr;
  vector[N_all] theta;
  vector[N_all] z_theta;
  matrix[N_all, nA] X;

  {
#include includes/impute_model/generate_X_CV.stan
    z_theta = a0 + X*a;
    theta = inv_logit(z_theta);
    {
      vector[N_all] RE_all[2];
      vector[N_all] bac_load;
      int B2[nB];
      for (i in 1:nB) B2[i] = B[i]+nXd;
      for (j in 1:2){
        RE_all[j,which(keptin)] = RE[j];
        for (n in which_not(keptin)) RE_all[j,n] = normal_rng(0, 1);
      }
      
      bac_load = b_HIV*X[:,1] + X[:,B2]*b + RE_all[1];
      vector[N_all] z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE_all[2] + square(RE_all[2])*quad_RE);
      vector[N_all] z_Mgit_RE  = z_Mgit[2]  + b_FE[2]*bac_load + b_RE[2]*(RE_all[2] + square(RE_all[2])*quad_RE);
      vector[N_all] z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE_all[2] + square(RE_all[2])*quad_RE);
      p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta .* inv_logit(z_Smear_RE);
      p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta .* inv_logit(z_Mgit_RE);
      p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta .* inv_logit(z_Xpert_RE);
      
#include includes/generated_quantities/pairwise_corr.stan
      
      for (n in 1:N_all){
        log_lik[n] = log_mix(theta[n],
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert_RE[n]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit_RE[n]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear_RE[n]),
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1])
        );
      }
    }
  }
}
