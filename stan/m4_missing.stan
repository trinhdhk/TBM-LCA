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
  // int obs_either_all[N_all];
  // int obs_either[N];
  
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariates 
  int nA = nX + nQ; // Number of coef
#include includes/cross_validation/transform_data_Y.stan
  int<lower=0, upper=1> obs_Smear[N] = obs_Smear_all[which(keptin)];
  int<lower=0, upper=1> obs_Mgit [N] = obs_Mgit_all [which(keptin)];
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
  real b_HIV_raw; //adjustment of RE with HIV Xd[,1]
  vector[nB] b_raw;
  vector<lower=0>[3] b_RE_raw;
  vector<lower=0>[3] b_FE_raw;
  vector[N] RE; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  ordered[2] z_obs;
 
  //Penalty terms
  vector<lower=0>[adapt_penalty[1]+adapt_penalty[2]] sp;
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
  int nu = 4;
   // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/a.stan
#include includes/main_prior/penalty.stan
  z_obs[1] ~ logistic(0, .5);
  z_obs[2] ~ logistic(logit(.99), 1);
  RE    ~ normal(0, 1);
  if (penalty_family == 0){
    b_raw     ~ student_t(nu, 0, 1);
    b_HIV_raw ~ student_t(nu, 0, 1);
    b_RE_raw  ~ student_t(nu, 0, 1);
    b_FE_raw  ~ student_t(nu, 0, 1);
  }
  if (penalty_family == 1){
    b_raw     ~ double_exponential(0, 1);
    b_HIV_raw ~ double_exponential(0, 1);
    b_RE_raw  ~ double_exponential(0, 1);
    b_FE_raw  ~ double_exponential(0, 1);
  }
  if (penalty_family == 2){
    b_raw     ~ normal(0, 1);
    b_HIV_raw ~ normal(0, 1);
    b_RE_raw  ~ normal(0, 1);
    b_FE_raw  ~ normal(0, 1);
  }
  for (n in 1:N){
    int N_Xd_miss = 3 - sum(obs_Xd[n, 1:3]);
#include includes/impute_model/impute_priors.loop_part.stan

    if (N_Xd_miss > 0){ //if there is some discrete variables missing
      int N_pattern = int_power(2, N_Xd_miss);
      vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,:], obs_Xd[n, 1:3], a[1:3]);
      vector[N_pattern] log_liks;
      pat_thetas[2] += a0 + dot_product(a[4:], X_compl[n]);
    
      //check if HIV is missing
      if (obs_Xd[n,1] == 1){
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
            
          real bac_load   = b_HIV*Xd_imp[n,1] + dot_product(b, Xc_imp[n,B]);
          real z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE[n] + square(RE[n])*quad_RE);
          real z_Mgit_RE  = z_Mgit[2]  + b_FE[2]*bac_load + b_RE[2]*(RE[n] + square(RE[n])*quad_RE);
          real z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE[n] + square(RE[n])*quad_RE);
          
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
            
          int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] > 0 ? 1 : 0;  
          log_liks[i] = logprob_theta + log_mix(theta, 
            // (bernoulli_logit_lpmf(obs_Smear[n] | z_obs[1,2]) + bernoulli_logit_lpmf(obs_Mgit[n] | z_obs[2,2]) + bernoulli_logit_lpmf(obs_Xpert[n] | z_obs[3,2]))/3 + 
            bernoulli_logit_lpmf(obs | z_obs[2]) + 
            ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
            // (bernoulli_logit_lpmf(obs_Smear[n] | z_obs[1,1]) + bernoulli_logit_lpmf(obs_Mgit[n] | z_obs[2,1]) + bernoulli_logit_lpmf(obs_Xpert[n] | z_obs[3,1]))/3 + 
            bernoulli_logit_lpmf(obs | z_obs[1]) + 
            ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
        }
      } else {
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
          
          vector[2] pat_bac_load[2] = get_patterns([Xd_imp[n,1]], {0}, [b_HIV]');
          vector[2] logprob_Y = pat_bac_load[1];
          vector[2] bac_load = pat_bac_load[2] + dot_product(b, Xc_imp[n,B]);
          
          vector[2] z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE[n] + square(RE[n])*quad_RE);
          vector[2] z_Mgit_RE  = z_Mgit[2]  + b_FE[2]*bac_load + b_RE[2]*(RE[n] + square(RE[n])*quad_RE);
          vector[2] z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE[n] + square(RE[n])*quad_RE);
          
          real ll_Smear[3] = 
            obs_Smear[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1]), bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2])} :
            { 0, 0, 0};
            
          real ll_Mgit[3] = 
            obs_Mgit[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]), bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2])} :
            { 0, 0, 0};
            
          real ll_Xpert[3] = 
            obs_Xpert[n] == 1 ?
            { bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[1]), bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[2])} :
            { 0, 0, 0};
            
          int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] > 0 ? 1 : 0;
          log_liks[i] = logprob_theta + log_mix(theta, 
            // (bernoulli_logit_lpmf(obs_Smear[n] | z_obs[1,2]) + bernoulli_logit_lpmf(obs_Mgit[n] | z_obs[2,2]) + bernoulli_logit_lpmf(obs_Xpert[n] | z_obs[3,2]))/3 +
            bernoulli_logit_lpmf(obs | z_obs[2]) +
            log_sum_exp(
              logprob_Y[1] + 
              ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
              logprob_Y[2] + 
              ll_Smear[3] + ll_Mgit[3] + ll_Xpert[3]
            ),
            // (bernoulli_logit_lpmf(obs_Smear[n] | z_obs[1,1]) + bernoulli_logit_lpmf(obs_Mgit[n] | z_obs[2,1]) + bernoulli_logit_lpmf(obs_Xpert[n] | z_obs[3,1]))/3 + 
            bernoulli_logit_lpmf(obs | z_obs[1]) +
            ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
        }
      }
      // Sum everything up
      target += log_sum_exp(log_liks);
      
    } else {
      // The normal way
      row_vector[nA] X = append_col(Xd_imp[n,:], X_compl[n,:]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load   = b_HIV*Xd_imp[n, 1] + dot_product(b, Xc_imp[n,B]);
      real z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE[n] + square(RE[n])*quad_RE);
      real z_Mgit_RE  = z_Mgit [2] + b_FE[2]*bac_load + b_RE[2]*(RE[n] + square(RE[n])*quad_RE);
      real z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE[n] + square(RE[n])*quad_RE);
      
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
    
      int obs = obs_Smear[n] + obs_Mgit[n] + obs_Xpert[n] > 0 ? 1 : 0;
      target += log_mix(theta, 
            // (bernoulli_logit_lpmf(obs_Smear[n] | z_obs[1,2]) + bernoulli_logit_lpmf(obs_Mgit[n] | z_obs[2,2]) + bernoulli_logit_lpmf(obs_Xpert[n] | z_obs[3,2]))/3 + 
            bernoulli_logit_lpmf(obs | z_obs[2]) +
            ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
            // (bernoulli_logit_lpmf(obs_Smear[n] | z_obs[1,1]) + bernoulli_logit_lpmf(obs_Mgit[n] | z_obs[2,1]) + bernoulli_logit_lpmf(obs_Xpert[n] | z_obs[3,1]))/3 + 
            bernoulli_logit_lpmf(obs | z_obs[1]) + 
            ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]);
    }
  }
}

generated quantities {
  vector[N_all] log_lik;
  vector[N_all] p_Smear;
  vector[N_all] p_Mgit;
  vector[N_all] p_Xpert;
  vector[N_all] theta;
  matrix[N_all, nX] X;

  {
#include includes/impute_model/generate_X_CV.stan
    
    theta = inv_logit(a0 + X*a);
    {
      vector[N_all] RE_all;
      vector[N_all] bac_load;
      int B2[nB];
      for (i in 1:nB) B2[i] = B[i]+nXd;
      RE_all[which(keptin)] = RE;
      for (n in which_not(keptin)) RE_all[n] = normal_rng(0, 1);
      bac_load = b_HIV*X[:,1] + X[:,B2]*b;
      vector[N_all] z_Smear_RE = z_Smear[2] + b_FE[1]*bac_load + b_RE[1]*(RE_all + square(RE_all)*quad_RE);
      vector[N_all] z_Mgit_RE  = z_Mgit[2]  + b_FE[2]*bac_load + b_RE[2]*(RE_all + square(RE_all)*quad_RE);
      vector[N_all] z_Xpert_RE = z_Xpert[2] + b_FE[3]*bac_load + b_RE[3]*(RE_all + square(RE_all)*quad_RE);
      p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta .* inv_logit(z_Smear_RE);
      p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta .* inv_logit(z_Mgit_RE);
      p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta .* inv_logit(z_Xpert_RE);
      
      for (n in 1:N_all){
        real ll_Smear[2] = 
          obs_Smear_all[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1]), bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear_RE[n]) } :
          { 0, 0 };
          
        real ll_Mgit[2] = 
          obs_Mgit_all[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]), bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit_RE[n]) } :
          { 0, 0 };
        
        real ll_Xpert[2] = 
          obs_Xpert_all[n] == 1 ?
          { bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]), bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert_RE[n]) } :
          { 0, 0 };
        int obs = obs_Smear_all[n] + obs_Mgit_all[n] + obs_Xpert_all[n] > 0 ? 1 : 0;  
        log_lik[n] = log_mix(theta[n],
        // (bernoulli_logit_lpmf(obs_Smear_all[n] | z_obs[1,2]) + bernoulli_logit_lpmf(obs_Mgit_all[n] | z_obs[2,2]) + bernoulli_logit_lpmf(obs_Xpert_all[n] | z_obs[3,2]))/3 + 
        bernoulli_logit_lpmf(obs | z_obs[2]) + 
        ll_Smear[2] + ll_Mgit[2] + ll_Xpert[2],
        // (bernoulli_logit_lpmf(obs_Smear_all[n] | z_obs[1,1]) + bernoulli_logit_lpmf(obs_Mgit_all[n] | z_obs[2,1]) + bernoulli_logit_lpmf(obs_Xpert_all[n] | z_obs[3,1]))/3 + 
        bernoulli_logit_lpmf(obs | z_obs[1]) + 
        ll_Smear[1] + ll_Mgit[1] + ll_Xpert[1]
        
        );
      }
    }
  }
}
