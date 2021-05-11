functions{
#include ./includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=1> nXc;     //Number of continuous X
  int<lower=1> nXd;     //Number of discrete X
  int<lower=1> nTd;     //Number of disc. aux covariates for imputation model
  int<lower=1> nTc;     //Number of cont. aux covariates for imputation model
  
  int<lower=0, upper=1> Y_Smear_all[N_all];
  int<lower=0, upper=1> Y_Mgit_all[N_all];
  int<lower=0, upper=1> Y_Xpert_all[N_all];
  
  real Xc_all[N_all, nXc]; //Continuous covariates
  int  Xd_all[N_all, nXd]; //Discrete covariates
  int  Td_all[N_all, nTd]; //Auxillary covariates - discrete
  real Tc_all[N_all, nTc]; //Auxillary covariates - continuous
  
  //Matrices of observation
  int<lower=0, upper=1> obs_Xc_all[N_all, nXc]; 
  int<lower=0, upper=1> obs_Xd_all[N_all, nXd];
  int<lower=0, upper=1> obs_Td_all[N_all, nTd]; 
  int<lower=0, upper=1> obs_Tc_all[N_all, nTc]; 
  
  // Hold-out for cross-validation
  int<lower=0, upper=1> keptin[N_all];
}

transformed data{
  int timestamp = 41961735; 
  //this is just a time stamp to force Stan to recompile the code and not used.  
  
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariates 

#include includes/cross_validation/transform_data_Y.stan
#include includes/cross_validation/transform_data_X.stan
#include includes/impute_model/transform_data.stan
}
 
parameters {
    // Parameters of the imputation model ---------------------------------------
#include /includes/impute_model/parameters.stan
  
  // Parameters of the logistics regression -----------------------------------
  real a0; //intercept
  real<lower=0> a_pos; // assuming HIV must have positive effect. 
  vector[nX - 1] a_; // Extra 1 is for sqrt(glu_ratio) = sqrt(csf_glu)/sqrt(bld_glu)
  real b_HIV; //adjustment of RE with HIV Xd[,1]
  real<lower=0> b_RE;
  vector[N] RE; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
}


transformed parameters {
  vector[nX] a = append_row(a_pos, a_); //Add HIV coef to the vector of coef
#include ./includes/impute_model/transform_parameters.stan
}


model {
  int nu = 5;
   // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
  
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/mN.stan
#include includes/main_prior/m_RE.stan

  for (n in 1:N){
    int N_Xd_miss = 3 - sum(obs_Xd[n, 1:3]);
#include /includes/impute_model/impute_priors.loop_part.stan

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
          
          real bac_load   = b_HIV*Xd_imp[n,1];
          real z_Smear_RE = z_Smear[2] + b_RE*(bac_load + RE[n]);
          real z_Mgit_RE  = z_Mgit[2]  + b_RE*(bac_load + RE[n]);
          real z_Xpert_RE = z_Xpert[2] + b_RE*(bac_load + RE[n]);
          
          log_liks[i] = logprob_theta + log_mix(theta, 
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      } else {
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
          
          vector[2] pat_bac_load[2] = get_patterns([Xd_imp[n,1]], {0}, [b_HIV]');
          vector[2] logprob_Y = pat_bac_load[1];
          vector[2] bac_load = pat_bac_load[2];
          
          vector[2] z_Smear_RE = z_Smear[2] + b_RE*(bac_load + RE[n]);
          vector[2] z_Mgit_RE  = z_Mgit[2]  + b_RE*(bac_load + RE[n]);
          vector[2] z_Xpert_RE = z_Xpert[2] + b_RE*(bac_load + RE[n]);
          
          log_liks[i] = logprob_theta + log_mix(theta, log_sum_exp(
            logprob_Y[1] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1])),
            logprob_Y[2] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2]))
            ),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      }
      // Sum everything up
      target += log_sum_exp(log_liks);
      
    } else {
      // The normal way
      row_vector[nX + 1] X = append_col(Xd_imp[n,:], X_compl[n,:]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load   = b_HIV*Xd_imp[n, 1];
      real z_Smear_RE = z_Smear[2] + b_RE*(bac_load + RE[n]);
      real z_Mgit_RE  = z_Mgit [2] + b_RE*(bac_load + RE[n]);
      real z_Xpert_RE = z_Xpert[2] + b_RE*(bac_load + RE[n]);
      
      target += log_mix(theta, 
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
    }
  }
}

generated quantities {
  vector[N_all] log_lik;
  vector[N_all] p_Smear;
  vector[N_all] p_Mgit;
  vector[N_all] p_Xpert;
  vector[N_all] theta;

  {
    matrix[N_all, nX + 1] X;
#include /includes/impute_model/generate_X_CV.stan
    
    theta = inv_logit(a0 + X*a);
    
    {
      vector[N_all] RE_all;
      vector[N_all] bac_load;
      RE_all[which(keptin)] = RE;
      for (n in which_not(keptin)) RE_all[n] = normal_rng(0, 1);
      bac_load = b_HIV*X[:,1] + RE_all;
      p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta .* inv_logit(z_Smear[2] + b_RE*bac_load);
      p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta .* inv_logit(z_Mgit[2]  + b_RE*bac_load);
      p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta .* inv_logit(z_Xpert[2] + b_RE*bac_load);
      
      for (n in 1:N_all){
        log_lik[n] = log_mix(theta[n],
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[2] + b_RE*bac_load[n]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[2] + b_RE*bac_load[n]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[2] + b_RE*bac_load[n]),
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1])
        );
      }
    }
  }
}
