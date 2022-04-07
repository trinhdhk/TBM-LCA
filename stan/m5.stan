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
  vector<lower=0>[3] b_HIV_raw; //adjustment of RE with HIV Xd[,1]
  // vector[3] b_cs_raw;
  matrix[nB,3] b_raw;
  vector<lower=0>[3] b_RE_raw;
  vector[N] RE; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  
#include includes/parameters/penalty.stan
}


transformed parameters {
#include includes/transform_parameters/a_variable_declaration.stan
  vector<lower=0>[3] b_HIV;
  // vector[3] b_cs;
  matrix[nB,3] b;
  vector<lower=0>[3] b_RE;
#include includes/impute_model/transform_parameters.stan
   {
#include includes/transform_parameters/penalty.stan
#include includes/transform_parameters/a_transform.stan
    b = b_raw * SP[2];
    b_HIV = b_HIV_raw * SP[2];
    // b_cs = b_cs_raw * SP[2];
    b_RE = b_RE_raw * SP[2];
  } 
}


model {
  int nu = 4;
  
  // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
    
    // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/a.stan
#include includes/main_prior/penalty.stan

    RE ~ normal(0,1);
    if (nB > 0) {
     int j = 1;
     for (i in B){
       b_raw[j]  ~ normal(0, inv(sd_X[i]));
       j += 1;
     }
    }
    
    if (penalty_family == 0){
      // to_vector(b_raw)  ~ student_t(nu, 0, 1);
      // b_cs_raw  ~ student_t(nu, 0,2);
      b_RE_raw ~ student_t(nu, 0, 1);
      b_HIV_raw ~ student_t(nu, 0,2);
    }
    if (penalty_family == 1){
      // to_vector(b_raw)  ~ double_exponential(0, 1);
      // b_cs_raw  ~ double_exponential(0, 2);
      b_RE_raw ~ double_exponential(0, 1);
      b_HIV_raw ~ double_exponential(0, 2);
    }
    if (penalty_family == 2){
      // to_vector(b_raw)  ~ normal(0, 1);
      // b_cs_raw   ~ normal(0, 2);
      b_RE_raw  ~ normal(0, 1);
      b_HIV_raw  ~ normal(0, 2);
    }
    
    for (n in 1:N){
       // If HIV is observed
    if (obs_Xd[n,1]==1){
#include includes/impute_model/impute_priors.observedHIV.stan
      int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
      vector[3] bac_load   = b_HIV*Xd_imp[n, 1] + to_vector(Xc_imp[n,B]' * b);
      real z_Smear_RE = z_Smear[2] + bac_load[1] + b_RE[1]*(RE[n] + square(RE[n])*quad_RE);
      real z_Mgit_RE  = z_Mgit [2] + bac_load[2] + b_RE[2]*(RE[n] + square(RE[n])*quad_RE);
      real z_Xpert_RE = z_Xpert[2] + bac_load[3] + b_RE[3]*(RE[n] + square(RE[n])*quad_RE);
      
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
      if (is_nan(
        multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,1] + cs_a[:,2]*Xc_imp[n,1], L_Omega_cs) +
        multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,2]*Xc_imp[n,1], L_Omega_cs) +
        multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,1] + mp_a[:,2]*Xc_imp[n,1], L_Omega_mp) +
        multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + cs_a[:,2]*Xc_imp[n,1], L_Omega_mp)))
        // This is to suppress the program from complaining at the start
        target += not_a_number();

#include includes/impute_model/impute_priors.unobsHIVpos.stan
#include includes/impute_model/impute_priors.unobsHIVneg.stan
  
      // if HIV is positive
      {
        int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
        
        vector[3] bac_load   = b_HIV + to_vector(Xc_imp[n,B]' * b);
        real z_Smear_RE = z_Smear[2] + bac_load[1] + b_RE[1]*(RE[n] + square(RE[n])*quad_RE);
        real z_Mgit_RE  = z_Mgit [2] + bac_load[2] + b_RE[2]*(RE[n] + square(RE[n])*quad_RE);
        real z_Xpert_RE = z_Xpert[2] + bac_load[3] + b_RE[3]*(RE[n] + square(RE[n])*quad_RE);
        // Symptoms or motor palsy is missing
        if (N_Xd_miss > 0){
          int N_pattern = int_power(2, N_Xd_miss);
          vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,2:3], obs_Xd[n,2:3], a[2:3]);
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
        vector[3] bac_load   = to_vector(Xc_imp[n,B]' * b);
        real z_Smear_RE = z_Smear[2] + bac_load[1] + b_RE[1]*(RE[n] + square(RE[n])*quad_RE);
        real z_Mgit_RE  = z_Mgit [2] + bac_load[2] + b_RE[2]*(RE[n] + square(RE[n])*quad_RE);
        real z_Xpert_RE = z_Xpert[2] + bac_load[3] + b_RE[3]*(RE[n] + square(RE[n])*quad_RE);
            
        // Symptoms or motor palsy is missing
        if (N_Xd_miss > 0){
          int N_pattern = int_power(2, N_Xd_miss);
          vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,2:3], obs_Xd[n, 2:3], a[2:3]);
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
      target += log_mix(p_HIV[n], 
      ll_HIV[1] + ll_z_mp[1] + ll_z_cs[1] + ll_Xc_imp_2[1],
      ll_HIV[2] + ll_z_mp[2] + ll_z_cs[2] + ll_Xc_imp_2[2]);
    }
  }
}

generated quantities {
  vector[N_all] log_lik;
  vector[N_all] p_Smear;
  vector[N_all] p_Mgit;
  vector[N_all] p_Xpert;
  vector[3] pairwise_corr;
  vector[N_all] z_theta;
  vector[N_all] theta;
  matrix[N_all, nA] X;

  {
#include includes/impute_model/generate_X_CV.stan
    
    z_theta = a0 + X*a;
    theta = inv_logit(z_theta);
    {
      vector[N_all] RE_all;
      vector[N_all] bac_load[3];
      int B2[nB];
      for (i in 1:nB) B2[i] = B[i]+nXd;
      RE_all[which(keptin)] = RE;
      for (n in which_not(keptin)) RE_all[n] = normal_rng(0, 1);
      bac_load[1] = b_HIV[1]*X[:,1] + X[:,B2]*b[:,1];
      bac_load[2] = b_HIV[2]*X[:,1] + X[:,B2]*b[:,2];
      bac_load[3] = b_HIV[3]*X[:,1] + X[:,B2]*b[:,3];
      vector[N_all] z_Smear_RE = z_Smear[2] + b_RE[1]*(RE_all + square(RE_all)*quad_RE) + bac_load[1];
      vector[N_all] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(RE_all + square(RE_all)*quad_RE) + bac_load[2];
      vector[N_all] z_Xpert_RE = z_Xpert[2] + b_RE[3]*(RE_all + square(RE_all)*quad_RE) + bac_load[3];
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
