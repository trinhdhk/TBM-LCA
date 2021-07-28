functions{
#include includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=0, upper=1> unsure_spc;
  
#include includes/data/X.stan
#include includes/data/Y.stan
#include includes/cross_validation/data.stan
#include includes/data/penalty.stan  
}

transformed data{
  // Penalty term adaptation
  int adapt_penalty[2];
  
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariates 
  int<lower=1, upper=4> Y_Xpert_C_all[N_all];
  int<lower=1, upper=4> Y_Xpert_C[sum(keptin)];

#include includes/cross_validation/transform_data_Y.stan
#include includes/cross_validation/transform_data_X.stan
#include includes/impute_model/transform_data.stan

  for (i in 1:2) adapt_penalty[i] = penalty_term[i] == 0 ? 1 : 0;
  for (n in 1:N_all) {
    if (Y_Xpert_all[n]==0) Y_Xpert_C_all[n] = 1;
    else {
      if (Tc_all[n,6]==1) Y_Xpert_C_all[n] = 2;
      if (Tc_all[n,7]==1) Y_Xpert_C_all[n] = 3;
      if (Tc_all[n,8]==1) Y_Xpert_C_all[n] = 4;
    }
  }
  Y_Xpert_C = Y_Xpert_C_all[which(keptin)];
}
 
parameters {
    // Parameters of the imputation model ---------------------------------------
#include includes/impute_model/parameters.stan
  
  // Parameters of the logistics regression -----------------------------------
  real a0; //intercept
  real<lower=0> a_pos; // assuming HIV must have positive effect. 
  vector[nX - 1] a_; // Extra 1 is for sqrt(glu_ratio) = sqrt(csf_glu)/sqrt(bld_glu)
  real b_HIV; //adjustment of RE with HIV Xd[,1]
  vector<lower=0>[3] b_RE;
  vector[N] RE; //base random effect;
  ordered[3] c_Xpert; //cut-off of xpert
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
 
 //Penalty terms
  real<lower=0> sp1[adapt_penalty[1]]; //sigma for the penalty
  real<lower=0> sp2[adapt_penalty[2]]; //sigma for the penalty
}


transformed parameters {
  vector[nX] a = append_row(a_pos, a_); //Add HIV coef to the vector of coef
#include includes/impute_model/transform_parameters.stan
}


model {
  int nu = 4;
  real SP[2] = {adapt_penalty[1] == 1 ? sp1[1] : penalty_term[1], adapt_penalty[2] == 1 ? sp2[1] : penalty_term[2]};
   // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
  
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/m.stan
#include includes/main_prior/m_RE.stan
#include includes/main_prior/penalty.stan
// Priors for model != m0
/**
if (penalty_family == 0) b_XpertLevel ~ student_t(nu, 0, SP[2]);
if (penalty_family == 1) b_XpertLevel ~ double_exponential(0, SP[2]);
if (penalty_family == 2) b_XpertLevel ~ normal(0, SP[2]);
**/

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
          
          real bac_load   = b_HIV*Xd_imp[n,1];
          real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n]);
          real z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n]);
          real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n]);
          
          log_liks[i] = logprob_theta + log_mix(theta, 
            ordered_logistic_lpmf(Y_Xpert_C[n] | z_Xpert_RE, c_Xpert) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
            ordered_logistic_lpmf(Y_Xpert_C[n] | z_Xpert[1], c_Xpert) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      } else {
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
          
          vector[2] pat_bac_load[2] = get_patterns([Xd_imp[n,1]], {0}, [b_HIV]');
          vector[2] logprob_Y = pat_bac_load[1];
          vector[2] bac_load = pat_bac_load[2];
          
          vector[2] z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n]);
          vector[2] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n]);
          vector[2] z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n]);
          
          log_liks[i] = logprob_theta + log_mix(theta, log_sum_exp(
            logprob_Y[1] + (ordered_logistic_lpmf(Y_Xpert_C[n] | z_Xpert_RE[1], c_Xpert) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1])),
            logprob_Y[2] + (ordered_logistic_lpmf(Y_Xpert_C[n] | z_Xpert_RE[2], c_Xpert) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2]))
            ),
            ordered_logistic_lpmf(Y_Xpert_C[n] | z_Xpert[1], c_Xpert) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      }
      // Sum everything up
      target += log_sum_exp(log_liks);
      
    } else {
      // The normal way
      row_vector[nX + 1] X = append_col(Xd_imp[n,:], X_compl[n,:]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load   = b_HIV*Xd_imp[n, 1];
      real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n]);
      real z_Mgit_RE  = z_Mgit [2] + b_RE[2]*(bac_load + RE[n]);
      real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n]);
      
      target += log_mix(theta, 
      ordered_logistic_lpmf(Y_Xpert_C[n] | z_Xpert_RE, c_Xpert) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
      ordered_logistic_lpmf(Y_Xpert_C[n] | z_Xpert[1], c_Xpert) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
    }
  }
}

 generated quantities {
   vector[N_all] log_lik;
   vector[N_all] p_Smear;
   vector[N_all] p_Mgit;
   matrix[N_all, 4] p_Xpert;
   vector[N_all] theta;
   matrix[N_all, nX] X;
 
   {
 #include includes/impute_model/generate_X_CV.stan
     
     theta = inv_logit(a0 + X*a);
     
     {
       vector[N_all] RE_all;
       vector[N_all] bac_load;
       RE_all[which(keptin)] = RE;
       for (n in which_not(keptin)) RE_all[n] = normal_rng(0, 1);
       bac_load = b_HIV*X[:,1] + RE_all;
       vector[N_all] z_Smear_RE = z_Smear[2] + b_RE[1]*bac_load;
       vector[N_all] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*bac_load;
       vector[N_all] z_Xpert_RE = z_Xpert[2] + b_RE[3]*bac_load;
       p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta .* inv_logit(z_Smear_RE);
       p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta .* inv_logit(z_Mgit_RE);
       p_Xpert[:,1] = (1 - theta) * (1 - inv_logit(z_Xpert[1] - c_Xpert[1])) + theta .* (1 - inv_logit(z_Xpert_RE - c_Xpert[1]));
       for (k in 2:3) {
         p_Xpert[:,k] = (1 - theta) * (inv_logit(z_Xpert[1] - c_Xpert[k-1]) - inv_logit(z_Xpert[1] - c_Xpert[k])) + theta .* (inv_logit(z_Xpert_RE - c_Xpert[k-1]) - inv_logit(z_Xpert_RE - c_Xpert[k]));
       }
       p_Xpert[:,4] = (1 - theta) * inv_logit(z_Xpert[1] - c_Xpert[3]) + theta .* inv_logit(z_Xpert_RE - c_Xpert[3]);
       for (n in 1:N_all){
         log_lik[n] = log_mix(theta[n],
         ordered_logistic_lpmf(Y_Xpert_C_all[n] | z_Xpert_RE[n],c_Xpert) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit_RE[n]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear_RE[n]),
         ordered_logistic_lpmf(Y_Xpert_C_all[n] | z_Xpert[1],c_Xpert) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1])
         );
       }
     }
   }
}
