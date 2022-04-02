functions{
#include includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=0, upper=1> unsure_spc;

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
  int nA = nX + nQ;

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
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  
#include includes/parameters/penalty.stan
}


transformed parameters {
#include includes/transform_parameters/a_variable_declaration.stan
  {
#include includes/transform_parameters/penalty.stan
#include includes/transform_parameters/a_transform.stan
  }
#include includes/impute_model/transform_parameters.stan
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

  for (n in 1:N){
    int N_Xd_miss = 3 - sum(obs_Xd[n, 1:3]);
#include includes/impute_model/impute_priors.loop_part.stan

    if (N_Xd_miss > 0){ //if there is some discrete variables missing
      int N_pattern = int_power(2, N_Xd_miss);
      vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,:], obs_Xd[n, 1:3], a[1:3]);
      vector[N_pattern] log_liks;
      pat_thetas[2] += a0 + dot_product(a[4:], X_compl[n]);
    
      for (i in 1:N_pattern){
        real logprob_theta = pat_thetas[1,i];
        real theta = inv_logit(pat_thetas[2,i]);
        
        log_liks[i] = logprob_theta + log_mix(theta, 
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
      }
      
      // Sum everything up
      target += log_sum_exp(log_liks);
      
    } else {
      // The normal way
      row_vector[nA] X = append_col(Xd_imp[n,:], X_compl[n,:]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      target += log_mix(theta, 
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
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
      p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta * inv_logit(z_Smear[2]);
      p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta * inv_logit(z_Mgit[2]);
      p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta * inv_logit(z_Xpert[2]);
      
      {
        vector[3] mu_i;
        vector[3] mu_ij;
        vector[N_all] p_Smear_Mgit;
        vector[N_all] p_Smear_Xpert;
        vector[N_all] p_Mgit_Xpert;
        {
          matrix[N_all,3] y;
          for (n in 1:N_all){
            y[n,1] = bernoulli_rng(p_Smear[n]);
            y[n,2] = bernoulli_rng(p_Mgit[n]);
            y[n,3] = bernoulli_rng(p_Xpert[n]);
          }
          for (i in 1:3) mu_i[i] = mean(y[:,i]);
        }
        
        p_Smear_Mgit = (1 - theta)*inv_logit(z_Smear[1])*inv_logit(z_Mgit[1]) +
          theta * inv_logit(z_Smear[2]) * inv_logit(z_Mgit[2]); 
        p_Smear_Xpert= (1 - theta)*inv_logit(z_Smear[1])*inv_logit(z_Xpert[1]) +
          theta * inv_logit(z_Smear[2]) * inv_logit(z_Xpert[2]); 
        p_Mgit_Xpert = (1 - theta)*inv_logit(z_Mgit[1])*inv_logit(z_Xpert[1]) +
          theta * inv_logit(z_Mgit[2]) * inv_logit(z_Xpert[2]); 
          {
            matrix[N_all,3] y;
            for (n in 1:N_all){
              y[n,1] = bernoulli_rng(p_Smear_Mgit[n]);
              y[n,2] = bernoulli_rng(p_Smear_Xpert[n]);
              y[n,3] = bernoulli_rng(p_Mgit_Xpert[n]);
            }
            for (i in 1:3) mu_ij[i] = mean(y[:,i]);
          }
    
        {
          int k = 1;
          for (i in 1:2){
            for (j in (i+1):3){
              pairwise_corr[k] = (mu_ij[k] - mu_i[i]*mu_i[j])/sqrt(mu_i[i]*(1-mu_i[i])*mu_i[j]*(1-mu_i[j]));
              k += 1;
            }
          }
        }
      }

      
      for (n in 1:N_all){
        log_lik[n] = log_mix(theta[n],
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[2]),
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1])
        );
      }
    }
  }
}
