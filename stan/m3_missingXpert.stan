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
  int obs_Smear2_all[N_all] = rep_array(1, N_all);
  int obs_Mgit2_all[N_all] = rep_array(1, N_all);
  int<lower=0, upper=1> obs_Smear[N] = obs_Smear2_all[which(keptin)];
  int<lower=0, upper=1> obs_Mgit [N] = obs_Mgit2_all [which(keptin)];
  int<lower=0, upper=1> obs_Xpert[N] = obs_Xpert_all[which(keptin)];
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
  vector[N] RE; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  ordered[2] z_obs;
 
#include includes/parameters/penalty.stan
}


transformed parameters {
#include includes/transform_parameters/a_variable_declaration.stan
  real<lower=0> b_HIV;
  vector[nB] b;
  vector<lower=0>[3] b_RE;
#include includes/impute_model/transform_parameters.stan
  
  {
#include includes/transform_parameters/penalty.stan
#include includes/transform_parameters/a_transform.stan
#include includes/transform_parameters/b_transform.stan   
  }
}


model {
  int nu = 5;
   // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/a.stan
#include includes/main_prior/m_RE.stan
#include includes/main_prior/penalty.stan
  {
    int j = 1;
    for (i in B){
      b_raw[j]  ~ normal(0, inv(sd_X[i]));
      j += 1;
    }
  }
  z_obs[1] ~ logistic(0, .5);
  z_obs[2] ~ logistic(logit(.99), 1);
  
  for (n in 1:N){
    // If HIV is observed
    if (obs_Xd[n,1]==1){
#include includes/impute_model/impute_priors.observedHIV.stan
      int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
      real bac_load   = b_HIV*Xd_imp[n,1] + dot_product(b, Xc_imp[n, B]);
      real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
      real z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
      real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
      
      // Symptoms or motor palsy is missing
      if (N_Xd_miss > 0){
        int N_pattern = int_power(2, N_Xd_miss);
        vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,2:3], obs_Xd[n, 2:3], a[2:3]);
        vector[N_pattern] log_liks;
        pat_thetas[2] += a0 + a[1]*Xd_imp[n,1] + dot_product(a[4:], X_compl[n]);
        
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
          
          real ll_Smear[2] = 
          obs_Smear[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE) } :
          { 0, 0 };
          
          real ll_Mgit[2] = 
          obs_Mgit[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) } :
          { 0, 0 };
          
          real ll_Xpert[2] = 
          obs_Xpert[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) } :
          { 0, 0 };
          
          int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] == 3 ? 1 : 0;  
          log_liks[i] = logprob_theta + log_mix(theta, 
          bernoulli_logit_lpmf(obs | z_obs[2]) + ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
          bernoulli_logit_lpmf(obs | z_obs[1]) + ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
        }
        target += log_sum_exp(log_liks);
        
        // The normal way
      } else {
        row_vector[nA] X = append_col(Xd_imp[n,:], X_compl[n,:]);
        real theta = inv_logit(a0 + dot_product(a, X));
        
        real ll_Smear[2] = 
        obs_Smear[n] == 1 ?
        { bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE) } :
        { 0, 0 };
        
        real ll_Mgit[2] = 
        obs_Mgit[n] == 1 ?
        { bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) } :
        { 0, 0 };
        
        real ll_Xpert[2] = 
        obs_Xpert[n] == 1 ?
        { bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) } :
        { 0, 0 };
        
        int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] == 3 ? 1 : 0;  
        target += log_mix(theta, 
        bernoulli_logit_lpmf(obs | z_obs[2]) + ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
        bernoulli_logit_lpmf(obs | z_obs[1]) + ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
        
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
        
        real bac_load   = b_HIV + dot_product(b, Xc_imp[n, B]);
        real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
        real z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
        real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
        // Symptoms or motor palsy is missing
        if (N_Xd_miss > 0){
          int N_pattern = int_power(2, N_Xd_miss);
          vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,2:3], obs_Xd[n,2:3], a[2:3]);
          vector[N_pattern] log_liks;
          pat_thetas[2] += a0 + a[1] + dot_product(a[4:], X_compl[n]);
          
          for (i in 1:N_pattern){
            real logprob_theta = pat_thetas[1,i];
            real theta = inv_logit(pat_thetas[2,i]);
            
            real ll_Smear[2] = 
            obs_Smear[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE) } :
            { 0, 0 };
            
            real ll_Mgit[2] = 
            obs_Mgit[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) } :
            { 0, 0 };
            
            real ll_Xpert[2] = 
            obs_Xpert[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) } :
            { 0, 0 };
            
            int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] == 3 ? 1 : 0;  
            log_liks[i] = logprob_theta + log_mix(theta, 
            bernoulli_logit_lpmf(obs | z_obs[2]) + ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
            bernoulli_logit_lpmf(obs | z_obs[1]) + ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
            
          }
          ll_HIV[1] = log_sum_exp(log_liks);
        } else {
          row_vector[nA-1] X = append_col(Xd_imp[n,2:3], X_compl[n,:]);
          real theta = inv_logit(a0 + a[1] + dot_product(a[2:nA], X));
          
          real ll_Smear[2] = 
          obs_Smear[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE) } :
          { 0, 0 };
          
          real ll_Mgit[2] = 
          obs_Mgit[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) } :
          { 0, 0 };
          
          real ll_Xpert[2] = 
          obs_Xpert[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) } :
          { 0, 0 };
          
          int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] == 3 ? 1 : 0;  
          ll_HIV[1] = log_mix(theta, 
          bernoulli_logit_lpmf(obs | z_obs[2]) + ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
          bernoulli_logit_lpmf(obs | z_obs[1]) + ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
        }
      }
      
      // If HIV is negative
      {
        int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
        real bac_load   = dot_product(b, Xc_imp[n, B]);
        real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
        real z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
        real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            
        // Symptoms or motor palsy is missing
        if (N_Xd_miss > 0){
          int N_pattern = int_power(2, N_Xd_miss);
          vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,2:3], obs_Xd[n, 2:3], a[2:3]);
          vector[N_pattern] log_liks;
          pat_thetas[2] += a0 + dot_product(a[4:], X_compl[n]);
          
          for (i in 1:N_pattern){
            real logprob_theta = pat_thetas[1,i];
            real theta = inv_logit(pat_thetas[2,i]);
            
             real ll_Smear[2] = 
            obs_Smear[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE) } :
            { 0, 0 };
            
            real ll_Mgit[2] = 
            obs_Mgit[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) } :
            { 0, 0 };
            
            real ll_Xpert[2] = 
            obs_Xpert[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) } :
            { 0, 0 };
            
            int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] == 3 ? 1 : 0;  
            log_liks[i] = logprob_theta + log_mix(theta, 
            bernoulli_logit_lpmf(obs | z_obs[2]) + ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2], 
            bernoulli_logit_lpmf(obs | z_obs[1]) + ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
            
          }
          ll_HIV[2] = log_sum_exp(log_liks);
        } else {
          row_vector[nA-1] X = append_col(Xd_imp[n,2:3], X_compl[n,:]);
          real theta = inv_logit(a0 + dot_product(a[2:nA], X));
          
          real ll_Smear[2] = 
          obs_Smear[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE) } :
          { 0, 0 };
          
          real ll_Mgit[2] = 
          obs_Mgit[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) } :
          { 0, 0 };
          
          real ll_Xpert[2] = 
          obs_Xpert[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) } :
          { 0, 0 };
          
          int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] == 3 ? 1 : 0;  
          ll_HIV[2] = log_mix(theta, 
            bernoulli_logit_lpmf(obs | z_obs[2]) + ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
            bernoulli_logit_lpmf(obs | z_obs[1]) + ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
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
  vector[N_all] theta;
  vector[N_all] z_theta;
  vector[N_all] Y0_theta[3];
  vector[N_all] Y00_theta[3];
  vector[N_all] Y000_theta;
  matrix[N_all, nA] X;
  vector[N_all] z_Smear_RE;
  vector[N_all] z_Mgit_RE;
  vector[N_all] z_Xpert_RE;

  {
#include includes/impute_model/generate_X_CV.stan
    
    z_theta = a0 + X*a;
    theta = inv_logit(z_theta);
    {
      vector[N_all] RE_all;
      vector[N_all] bac_load;
      int B2[nB];
      for (i in 1:nB) B2[i] = B[i]+nXd;
      RE_all[which(keptin)] = RE;
      for (n in which_not(keptin)) RE_all[n] = normal_rng(0, 1);
      bac_load = b_HIV*X[:,1] + X[:, B2]*b + RE_all + square(RE_all)*quad_RE;
      z_Smear_RE = z_Smear[2] + b_RE[1]*bac_load;
      z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*bac_load;
      z_Xpert_RE = z_Xpert[2] + b_RE[3]*bac_load;
      p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta .* inv_logit(z_Smear_RE);
      p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta .* inv_logit(z_Mgit_RE);
      p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta .* inv_logit(z_Xpert_RE);
      
#include includes/generated_quantities/pairwise_corr.stan
#include includes/generated_quantities/post_test_prob.stan
      
      for (n in 1:N_all){
        real ll_Smear[2] = 
          obs_Smear2_all[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear_RE[n]) } :
          { 0, 0 };
          
        real ll_Mgit[2] = 
          obs_Mgit2_all[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit_RE[n]) } :
          { 0, 0 };
        
        real ll_Xpert[2] = 
          obs_Xpert_all[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert_RE[n]) } :
          { 0, 0 };
        int obs = obs_Smear2_all[n] + obs_Mgit2_all[n] + obs_Xpert_all[n] == 3 ? 1 : 0;  
        log_lik[n] = log_mix(theta[n],
        bernoulli_logit_lpmf(obs | z_obs[2]) + 
        ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
        bernoulli_logit_lpmf(obs | z_obs[1]) + 
        ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]
        
        );
      }
    }
  }
}
