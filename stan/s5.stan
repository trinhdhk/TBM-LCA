functions{
#include includes/functions.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=1> nXc;     //Number of continuous X
  int<lower=1> nXd;     //Number of discrete X
  real Xc_all[N_all, nXc]; //Continuous covariates
  int  Xd_all[N_all, nXd]; //Discrete covariates
  int  Y_all[N_all]; //class
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
#include includes/cross_validation/data.stan
}

transformed data{
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariates 
  int<lower=1> N = sum(keptin);
  int<lower=0> N_valid = N_all - N;
  real Xc[N, nXc] = Xc_all[which(keptin),:];
  int  Xd[N, nXd] = Xd_all[which(keptin),:];
  int  Y[N]       = Y_all[which(keptin)];
  real sd_X[nXc];
  int adapt_penalty = penalty_term == 0 ? 1 : 0;
  for (i in 1:nXc) sd_X[i] = sd(Xc[:,i]);
}
 
parameters {
    // Parameters of the imputation model ---------------------------------------
  //Penalty terms
  vector<lower=0>[adapt_penalty] sp;
  real a0; //intercept
  vector<lower=0>[1] a_pos;
  vector[nX-1] a_; //coefs
}

transformed parameters{
  vector[nX] a_raw;
  vector[nX] a;
  
  {
    real SP = (adapt_penalty == 1) ? sp[1] : penalty_term;
    a_raw = append_row(a_pos, a_);
    a = a_raw * SP;
  }
}

model {
  int nu = 4;
  matrix[N, nX] X = append_col(to_matrix(Xd), to_matrix(Xc));
   // Imputation model ---------------------------------------------------------
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/penalty.stan
  a0 ~ student_t(nu, 0, 5);
  
 if (penalty_family == 0){
    a_raw[1:nXd] ~ student_t(nu, 0, 2);
    for (i in 1:nXc){
      a_raw[nXd+i] ~ student_t(nu, 0, inv(sd_X[i]));
    }
  }
  if (penalty_family == 1){
    a_raw[1:nXd] ~ double_exponential(0, 2);
     for (i in 1:nXc){
      a_raw[nXd+i] ~ double_exponential(0, inv(sd_X[i]));
    }
  }
  if (penalty_family == 2){
    a_raw[1:nXd] ~ normal(0, 2);
     for (i in 1:nXc){
      a_raw[nXd+i] ~ normal(0, inv(sd_X[i]));
    }
  }
  
  Y ~ bernoulli_logit_glm(X, a0, a);
  
}

generated quantities {
  vector[N_all] log_lik;
  vector[N_all] lambda;
  vector[N_all] theta;
  
  {
    matrix[N_all, nX] X = append_col(to_matrix(Xd_all), to_matrix(Xc_all));
    lambda = a0 + X*a; 
    theta = inv_logit(lambda);
    for (n in 1:N_all){
      log_lik[n] = bernoulli_logit_glm_lpmf(Y_all | X, a0, a);
    }
  }
}
