functions{
#include includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=0> nB; //Number of added RE  
  int<lower=1> B[nB];
  int<lower=0, upper=1> unsure_spc;
  int<lower=0, upper=1> quad_RE;
  int<lower=0> nA_neg;
  int<lower=0> nA_pos;
  int<lower=0> A_neg[nA_neg]; //position of negative a(s)
  int<lower=0> A_pos[nA_pos]; //position of positive a(s)
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

#include includes/cross_validation/transform_data_Y.stan
#include includes/cross_validation/transform_data_X.stan
#include includes/impute_model/transform_data.stan

for (i in 1:2) adapt_penalty[i] = penalty_term[i] == 0 ? 1 : 0;
}
 
parameters {
    // Parameters of the imputation model ---------------------------------------
#include includes/impute_model/parameters.stan
  
  // Parameters of the logistics regression -----------------------------------
  real a0; //intercept
  vector<lower=0>[nA_pos] a_pos; // assuming HIV must have positive effect. 
  vector<upper=0>[nA_neg] a_neg;
  vector[nX - nA_pos - nA_neg] a_;
  real b_HIV; //adjustment of RE with HIV Xd[,1]
  real b_cs;
  vector[nB] b;
  vector<lower=0>[3] b_RE_raw;
  vector[N] RE; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
 
  //Penalty terms
  vector<lower=0>[adapt_penalty[1]+adapt_penalty[2]] sp;
}


transformed parameters {
  vector[nX] a; //Add HIV coef to the vector of coef
  vector<lower=0>[3] b_RE;
#include includes/impute_model/transform_parameters.stan
  {
    int nA_posneg = nA_pos + nA_neg;
    int A_[nX - nA_posneg];
    int A_posneg[nA_posneg] = append_array(A_pos, A_neg);
    int A_posneg_asc[nA_posneg] = sort_asc(A_posneg);
    int J = 1;
    int k = 1;
    for (i in 1:nX) {
      int match = 0;
      if (J <= (nA_posneg))
        for (j in J:(nA_posneg)){
          if (A_posneg_asc[j] == i) {
            match = 1;
            J += 1;
            break;
          }
        }
      if (match == 0) {
        A_[k] = i;
        k += 1;
      }
    }
    a[A_posneg] = append_row(a_pos, a_neg);
    a[A_] = a_;
  }
  
  {
    real SP[2];
    SP[1] = (adapt_penalty[1] == 1) ? sp[1] : penalty_term[1];
    SP[2] = (penalty_term[2] == 0) ? sp[num_elements(sp)] : (penalty_term[2] == -1) ? SP[1] : penalty_term[2];
    b_RE  = b_RE_raw/mean(fabs(append_row(b, b_HIV)))*SP[2];
  }
}


model {
  int nu = 4;
  real SP[2];
   // Imputation model ---------------------------------------------------------
#include includes/impute_model/variables_declaration.stan 
#include includes/impute_model/impute_priors.main_part.stan 
  
  SP[1] = (adapt_penalty[1] == 1) ? sp[1] : penalty_term[1];
  SP[2] = (penalty_term[2] == 0) ? sp[num_elements(sp)] : (penalty_term[2] == -1) ? SP[1] : penalty_term[2];
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/m.stan
#include includes/main_prior/m_RE.stan
#include includes/main_prior/penalty.stan
  if (penalty_family == 0){
    b  ~ student_t(nu, 0, SP[2]);
    b_cs  ~ student_t(nu, 0, SP[2]);
  }
  if (penalty_family == 1){
    b  ~ double_exponential(0, SP[2]);
    b_cs  ~ double_exponential(0, SP[2]);
  }
  if (penalty_family == 2){
    b  ~ normal(0, SP[2]);
    b_cs  ~ normal(0, SP[2]);
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
          
          if (obs_Xd[n,2] == 1){
            
            real bac_load   = b_HIV*Xd_imp[n,1] + b_cs*Xd_imp[n,2] + dot_product(b, Xc_imp[n,B]);
            real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            real z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
          
            log_liks[i] = logprob_theta + log_mix(theta, 
              bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
              bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          } else {
            vector[2] pat_bac_load[2] = get_patterns([Xd_imp[n,2]], {0}, [b_cs]');
            vector[2] logprob_Y = pat_bac_load[1];
            vector[2] bac_load = pat_bac_load[2] + b_HIV*Xd_imp[n,1] + dot_product(b, Xc_imp[n,B]);
          
            vector[2] z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            vector[2] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            vector[2] z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
          
            log_liks[i] = logprob_theta + log_mix(theta, log_sum_exp(
            logprob_Y[1] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1])),
            logprob_Y[2] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2]))
            ),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          }
        }
      } else {
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
          
          if (obs_Xd[n,2] == 1){
          
            vector[2] pat_bac_load[2] = get_patterns([Xd_imp[n,1]], {0}, [b_HIV]');
            vector[2] logprob_Y = pat_bac_load[1];
            vector[2] bac_load = pat_bac_load[2] + b_cs*Xd_imp[n,2] + dot_product(b, Xc_imp[n,B]);
            
            vector[2] z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            vector[2] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            vector[2] z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            
            log_liks[i] = logprob_theta + log_mix(theta, log_sum_exp(
              logprob_Y[1] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1])),
              logprob_Y[2] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2]))
              ),
              bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          } else {
            
            vector[4] pat_bac_load[2] = get_patterns(Xd_imp[n,1:2], {0, 0}, [b_HIV, b_cs]');
            vector[4] logprob_Y = pat_bac_load[1];
            vector[4] bac_load = pat_bac_load[2] + dot_product(b, Xc_imp[n,B]);
            
            vector[4] z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            vector[4] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            vector[4] z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
            
            log_liks[i] = logprob_theta + log_mix(theta, log_sum_exp([
              logprob_Y[1] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1])),
              logprob_Y[2] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2])),
              logprob_Y[3] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[3]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[3]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[3])),
              logprob_Y[4] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[4]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[4]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[4]))
              ]),
              bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
          }
        }
      }
      // Sum everything up
      target += log_sum_exp(log_liks);
      
    } else {
      // The normal way
      row_vector[nX] X = append_col(Xd_imp[n,:], X_compl[n,:]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load   = b_HIV*Xd_imp[n, 1] + b_cs*Xd_imp[n,2] + dot_product(b, Xc_imp[n,B]);
      real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n] + square(RE[n])*quad_RE);
      real z_Mgit_RE  = z_Mgit [2] + b_RE[2]*(bac_load + RE[n] + square(RE[n])*quad_RE);
      real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n] + square(RE[n])*quad_RE);
      
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
      bac_load = b_HIV*X[:,1] + b_cs*X[:,2] + RE_all + X[:,B2]*b;
      for (n in 1:N_all) bac_load += square(RE_all[n])*quad_RE;
      vector[N_all] z_Smear_RE = z_Smear[2] + b_RE[1]*bac_load;
      vector[N_all] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*bac_load;
      vector[N_all] z_Xpert_RE = z_Xpert[2] + b_RE[3]*bac_load;
      p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta .* inv_logit(z_Smear_RE);
      p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta .* inv_logit(z_Mgit_RE);
      p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta .* inv_logit(z_Xpert_RE);
      
      for (n in 1:N_all){
        log_lik[n] = log_mix(theta[n],
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert_RE[n]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit_RE[n]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear_RE[n]),
        bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1])
        );
      }
    }
  }
}
