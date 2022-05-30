functions{
#include includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=0, upper=1> unsure_spc;
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
#include includes/impute_model/transform_parameters.stan
  {
#include includes/transform_parameters/penalty.stan
#include includes/transform_parameters/a_transform.stan
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

  for (n in 1:N){
     if (obs_Xd[n,1]==1){
#include includes/impute_model/impute_priors.observedHIV.stan
      int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
      
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
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
        target += log_sum_exp(log_liks);
        // The normal way
      } else {
        row_vector[nA] X = append_col(Xd_imp[n,:], X_compl[n,:]);
        real theta = inv_logit(a0 + dot_product(a, X));
        
        target += log_mix(theta, 
        bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
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
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          }
          ll_HIV[1] = log_sum_exp(log_liks);
        } else {
          row_vector[nA-1] X = append_col(Xd_imp[n,2:3], X_compl[n,:]);
          real theta = inv_logit(a0 + a[1] + dot_product(a[2:nA], X));
          
          ll_HIV[1] = log_mix(theta, 
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      }
      
      // If HIV is negative
      {
        int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
            
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
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          }
          ll_HIV[2] = log_sum_exp(log_liks);
        } else {
          row_vector[nA-1] X = append_col(Xd_imp[n,2:3], X_compl[n,:]);
          real theta = inv_logit(a0 + dot_product(a[2:nA], X));
          
          ll_HIV[2] = log_mix(theta, 
          bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2]),
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
  vector[N_all] Y0_theta[3];
  vector[N_all] Y00_theta[3];
  vector[N_all] Y000_theta;

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
      
      {
        real pY000_theta = inv_logit(-z_Smear[2]) * inv_logit(-z_Mgit[2]) * inv_logit(-z_Xpert[2]);
        real pY000_1mtheta = inv_logit(-z_Smear[1]) * inv_logit(-z_Mgit[1]) * inv_logit(-z_Xpert[1]);
        real pY00_theta[3];
        real pY00_1mtheta[3];
        
        pY00_theta[1] = inv_logit(-z_Smear[2]) * inv_logit(-z_Mgit[2]);
        pY00_theta[2] = inv_logit(-z_Smear[2]) * inv_logit(-z_Xpert[2]);
        pY00_theta[3] = inv_logit(-z_Mgit[2])  * inv_logit(-z_Xpert[2]);
        
        pY00_1mtheta[1] = inv_logit(-z_Smear[1]) * inv_logit(-z_Mgit[1]) ;
        pY00_1mtheta[2] = inv_logit(-z_Smear[1]) * inv_logit(-z_Xpert[1]) ;
        pY00_1mtheta[3] = inv_logit(-z_Mgit[1]) * inv_logit(-z_Mgit[1]) ;
        
        Y0_theta[1] = (theta * inv_logit(-z_Smear[2])) ./ (1 - p_Smear);
        Y0_theta[2] = (theta * inv_logit(-z_Mgit[2])) ./ (1 - p_Mgit);
        Y0_theta[3] = (theta * inv_logit(-z_Xpert[2])) ./ (1 - p_Xpert);
        
        for (i in 1:3){
          Y00_theta[i] = (theta * pY00_theta[i]) ./ (theta * pY00_theta[i] + (1 - theta) * pY00_1mtheta[i]);
        }
        Y000_theta = (theta * pY000_theta) ./ (theta * pY000_theta + (1 - theta) * pY000_1mtheta);
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
