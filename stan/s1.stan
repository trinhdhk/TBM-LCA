functions{
#include includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
#include includes/data/X.stan
#include includes/data/Y.stan
  //Customisable penalty
  /*
  Different penalty values will have different prior families:
    0 = student_t(4);
    1 = double_exponential();
    2 = normal();
  
  penalty: if == 0 will perform auto adaptation for penalty terms,
  currently, RE model and prevalence will have different penalty if adpatation
  is turned on.
  if != 0 will create manual penalty.
  */
  int<lower=0, upper=2> penalty_family;
  real<lower=0> penalty_term;

  vector[N_all] mu_ztheta;
  vector<lower=0>[N_all] sd_ztheta;
  int<lower=0> nX_extra;
  matrix[N_all, nX_extra] X_extra_all;
#include includes/cross_validation/data.stan
}

transformed data{
  // int N = N_all;
  int nQ = 0;
  int quad_Xc_idx[nQ];
  int idxS[10] = {1,2,3,4,5,6,7,10,17};
  // Penalty term adaptation
  int adapt_penalty = penalty_term == 0 ? 1 : 0;
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariates 
  int nA = nX + nQ;
#include includes/cross_validation/transform_data_Y.stan
#include includes/cross_validation/transform_data_X.stan
  matrix[N, nX_extra] X_extra = X_extra_all[which(keptin)];
#include includes/impute_model/transform_data.stan
}
 
parameters {
    // Parameters of the imputation model ---------------------------------------
#include includes/impute_model/parameters.stan
  
  //Penalty terms
  vector<lower=0>[adapt_penalty] sp;
  
  
  //Full coef
  // vector[N] ztheta;
  //Simplified coefs;
  real a0; //intercept
  vector<lower=0>[1] a_pos;
  vector[8+nX_extra] a_; //coefs
  // real<lower=0> sigma_ztheta;
}

transformed parameters{
  vector[9+nX_extra] a_raw;
  vector[9+nX_extra] a;
  
  {
    real SP = (adapt_penalty == 1) ? sp[1] : penalty_term;
    a_raw = append_row(a_pos, a_);
    a = a_raw * SP;
  }
#include includes/impute_model/transform_parameters.stan
}

model {
  int nu = 4;
   // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
#include includes/main_prior/penalty.stan
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/penalty.stan
  a0 ~ student_t(nu, 0, 5);
  
 if (penalty_family == 0){
    int j = 1;
    a_raw[idxS[1:7]] ~ student_t(nu, 0, 2);
    for (i in {10-nXd, 17-nXd}){
      a_raw[j] ~ student_t(nu, 0, inv(sd_X[i]));
      j += 1;
    }
    for (i in 10:(9+nX_extra)) a_raw[i] ~ student_t(nu, 0, 2);
  }
  if (penalty_family == 1){
    int j = 1;
    a_raw[idxS[1:7]] ~ double_exponential(0, 2);
    for (i in {10-nXd, 17-nXd}){
      a_raw[j] ~ double_exponential(0, inv(sd_X[i]));
      j += 1;
    }
    for (i in 10:(9+nX_extra)) a_raw[i] ~ double_exponential(0, 2);
  }
  if (penalty_family == 2){
    int j = 1;
    a_raw[idxS[1:7]] ~ normal(0, 2);
    for (i in {10-nXd, 17-nXd}){
      a_raw[j] ~ normal(0, inv(sd_X[i]));
      j += 1;
    }
    for (i in 10:(9+nX_extra)) a_raw[i] ~ normal(0, 2);
  }
  // ztheta ~ logistic(mu_ztheta[which(keptin)], sd_ztheta[which(keptin)] * sqrt(3) / pi());
  // target += 100*logistic_lpdf(ztheta | mu_ztheta[which(keptin)], sd_ztheta[which(keptin)] * sqrt(3) / pi());
  // sigma_ztheta ~ normal(0,1);
  
  for (n in 1:N){
    // If HIV is observed
    if (obs_Xd[n,1]==1){
#include includes/impute_model/impute_priors.observedHIV.stan
      int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
      
      // Symptoms or motor palsy is missing
      if (N_Xd_miss > 0){
        int N_pattern = int_power(2, N_Xd_miss);
        vector[N_pattern] pat_lambda[2] = get_patterns(Xd_imp[n,:], obs_Xd[n, 2:3], a[2:3]);
        vector[N_pattern] log_liks;
        int idxS_2[7] = {1,2,3,4,7,14};
        pat_lambda[2] += a0 + a[1]*Xd_imp[n,1] + dot_product(a[4:9], X_compl[n, idxS_2]) + dot_product(a[10:], X_extra[n,:]);
    
        for (i in 1:N_pattern){
          // log_liks[i] = pat_lambda[1,i] + logistic_lpdf(ztheta[n] | pat_lambda[2,i], 1);
          log_liks[i] = pat_lambda[1,i] + logistic_lpdf(pat_lambda[2,i] | mu_ztheta[which(keptin)][n], sd_ztheta[which(keptin)][n] * sqrt(3) / pi());
        }
        // Sum everything up
        target += log_sum_exp(log_liks);
        
        // The normal way
      } else {
        row_vector[nX] X = append_col(Xd_imp[n,:], X_compl[n,:]);
        real lambda = a0 + dot_product(a[1:9], X[idxS]) + dot_product(a[10:],X_extra[n,:]); 
        // print(logistic_lpdf(ztheta[n] | lambda, sigma_ztheta));
        target += logistic_lpdf(lambda | mu_ztheta[which(keptin)][n], sd_ztheta[which(keptin)][n] * sqrt(3) / pi());
        // target += logistic_lpdf(ztheta[n] | lambda, 1);
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
          vector[N_pattern] pat_lambda[2] = get_patterns(Xd_imp[n,:], obs_Xd[n, 2:3], a[2:3]);
          vector[N_pattern] log_liks;
          int idxS_2[7] = {1,2,3,4,7,14};
          pat_lambda[2] += a0 + a[1] + dot_product(a[4:9], X_compl[n, idxS_2]) + dot_product(a[10:], X_extra[n,:]);
          
          for (i in 1:N_pattern){
            // log_liks[i] = pat_lambda[1,i] + logistic_lpdf(ztheta[n] | pat_lambda[2,i], 1);
            log_liks[i] = pat_lambda[1,i] + logistic_lpdf(pat_lambda[2,i] | mu_ztheta[which(keptin)][n], sd_ztheta[which(keptin)][n] * sqrt(3) / pi());
          }
          // Sum everything up
          ll_HIV[1] = log_sum_exp(log_liks);
        } else {
          row_vector[nX] X = append_col(Xd_imp[n,:], X_compl[n,:]);
          real lambda = a0 + a[1] + dot_product(a[2:9], X[idxS[2:]]) + dot_product(a[10:],X_extra[n,:]); 
          // print(logistic_lpdf(ztheta[n] | lambda, sigma_ztheta));
          // ll_HIV[1] = logistic_lpdf(ztheta[n] | lambda, 1);
          ll_HIV[1] = logistic_lpdf(lambda | mu_ztheta[which(keptin)][n], sd_ztheta[which(keptin)][n] * sqrt(3) / pi());
        }
      }
      
      // If HIV is negative
      {
        
        int N_Xd_miss = 2 - sum(obs_Xd[n, 2:3]);
        
        // Symptoms or motor palsy is missing
        if (N_Xd_miss > 0){
          int N_pattern = int_power(2, N_Xd_miss);
          vector[N_pattern] pat_lambda[2] = get_patterns(Xd_imp[n,:], obs_Xd[n, 2:3], a[2:3]);
          vector[N_pattern] log_liks;
          int idxS_2[7] = {1,2,3,4,7,14};
          pat_lambda[2] += a0 + dot_product(a[4:9], X_compl[n, idxS_2]) + dot_product(a[10:],X_extra[n,:]);
          
          for (i in 1:N_pattern){
            // log_liks[i] = pat_lambda[1,i] + logistic_lpdf(ztheta[n] | pat_lambda[2,i], 1);
            log_liks[i] = pat_lambda[1,i] + logistic_lpdf(pat_lambda[2,i] | mu_ztheta[which(keptin)][n], sd_ztheta[which(keptin)][n] * sqrt(3) / pi());
          }
          // Sum everything up
          ll_HIV[2] = log_sum_exp(log_liks);
        } else {
          row_vector[nX] X = append_col(Xd_imp[n,:], X_compl[n,:]);
          real lambda = a0 + dot_product(a[2:9], X[idxS[2:]]) + dot_product(a[10:],X_extra[n,:]); 
          // print(logistic_lpdf(ztheta[n] | lambda, sigma_ztheta));
          // ll_HIV[2] = logistic_lpdf(ztheta[n] | lambda, 1);
          ll_HIV[2] = logistic_lpdf(lambda | mu_ztheta[which(keptin)][n], sd_ztheta[which(keptin)][n] * sqrt(3) / pi());
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
  // vector[N_all] ztheta_all;
  vector[N_all] lambda;
  matrix[N_all, nX] X;
  
  {
    
#include includes/impute_model/generate_X_CV.stan
    
    // ztheta_all[which(keptin)] = ztheta;
    // ztheta_all[which_not(keptin)] = to_vector(logistic_rng(mu_ztheta[which_not(keptin)], sd_ztheta[which_not(keptin)] * sqrt(3) / pi()));
    lambda = a0 + X[:, idxS]*a[1:9] + X_extra_all*a[10:];
    
    for (n in 1:N_all){
      // log_lik[n] = logistic_lpdf(ztheta_all[n] | lambda[n], 1);
      log_lik[n] = logistic_lpdf(lambda[n] | mu_ztheta[n], sd_ztheta[n] * sqrt(3) / pi());
    }
  }
}
