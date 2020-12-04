functions{
#include /includes/functions.stan
}

data {
  int<lower=1> N;   //Number of patient
  int<lower=1> nXc; //Number of continuous X
  int<lower=1> nXd; //Number of discrete X
  int<lower=1> nTd; //Number of disc. aux covariates for imputation model
  
  int<lower=0, upper=1> Y_Smear[N];
  int<lower=0, upper=1> Y_Mgit[N];
  int<lower=0, upper=1> Y_Xpert[N];
  
  real Xc[N, nXc]; //Continuous covariates
  int  Xd[N, nXd]; //Discrete covariates
  int  Td[N, nTd]; //Auxillary covariates - discrete
  
  //Matrices of observation
  int<lower=0, upper=1> obs_Xc[N, nXc]; 
  int<lower=0, upper=1> obs_Xd[N, nXd];
  int<lower=0, upper=1> obs_Td[N, nTd]; 
}

transformed data {
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariates
  
  // * Clinical symptoms Td[1:3] ----------------------------------------------
  int Td_cs[N,3] = Td[:,1:3];
  int obs_cs[N,3] = obs_Td[:,1:3];
  int<lower=0> N_miss_cs;
  int<lower=1,upper=N> n_miss_cs[(N*3) - sum2d(obs_cs)];
  int<lower=1,upper=3> d_miss_cs[size(n_miss_cs)];
  int<lower=0> N_pos_cs;
  int<lower=1,upper=N> n_pos_cs[sum2d_with_missing(Td_cs, obs_cs)];
  int<lower=1,upper=3> d_pos_cs[size(n_pos_cs)];
  int<lower=0> N_neg_cs;
  int<lower=1,upper=N> n_neg_cs[sum2d(obs_cs) - size(n_pos_cs)];
  int<lower=1,upper=3> d_neg_cs[size(n_neg_cs)];
  
   // * Motor palsy Td[4:6] ---------------------------------------------------
  int Td_mp[N,3] = Td[,4:6];
  int obs_mp[N,3] = obs_Td[,4:6];
  int<lower=0> N_miss_mp;
  int<lower=1,upper=N> n_miss_mp[(N*3) - sum2d(obs_mp)];
  int<lower=1,upper=3> d_miss_mp[size(n_miss_mp)];
  int<lower=0> N_pos_mp;
  int<lower=1,upper=N> n_pos_mp[sum2d_with_missing(Td_mp, obs_mp)];
  int<lower=1,upper=3> d_pos_mp[size(n_pos_mp)];
  int<lower=0> N_neg_mp;
  int<lower=1,upper=N> n_neg_mp[sum2d(obs_mp) - size(n_pos_mp)];
  int<lower=1,upper=3> d_neg_mp[size(n_neg_mp)];
  
  // * CSF laboratory test Xc[3:7] --------------------------------------------
  int obs_bld_glu = sum(obs_Xc[:,3]);
  int obs_csf_glu = sum(obs_Xc[:,4]);
  int obs_csf_other = sum2d(obs_Xc[:,5:7]);
  
  // * Var assignments --------------------------------------------------------
  N_pos_cs  = size(n_pos_cs);
  N_neg_cs  = size(n_neg_cs);
  N_miss_cs = size(n_miss_cs);
  {
    int i = 1;
    int j = 1;
    int k = 1;
    for (n in 1:N) {
      for (d in 1:3) {
        if (!obs_cs[n,d]){
          n_miss_cs[k] = n;
          d_miss_cs[k] = d;
          k += 1;
        } else if (Td_cs[n,d]) {
          n_pos_cs[i] = n;
          d_pos_cs[i] = d;
          i += 1;
        } else {
          n_neg_cs[j] = n;
          d_neg_cs[j] = d;
          j += 1;
        }
      }
    }
  }
  
  N_pos_mp = size(n_pos_mp);
  N_neg_mp = size(n_neg_mp);
  N_miss_mp = size(n_miss_mp);
  {
    int i = 1;
    int j = 1;
    int k = 1;
    for (n in 1:N) {
      for (d in 1:3) {
        if (!obs_mp[n,d]){
          n_miss_mp[k] = n;
          d_miss_mp[k] = d;
          k += 1;
        } else if (Td_mp[n,d]) {
          n_pos_mp[i] = n;
          d_pos_mp[i] = d;
          i += 1;
        } else {
          n_neg_mp[j] = n;
          d_neg_mp[j] = d;
          j += 1;
        }
      }
    }
  }
}

parameters{
  // * Impute HIV Xd[1] -------------------------------------------------------
  real HIV_a0;
  
  // Impute symptoms Td[1:3] for components variables -------------------------
  // Xd[2] is the combined one aka clin_symptoms
  vector[3] cs_a0;
  matrix[3,2] cs_a; //For HIV & illness days
  cholesky_factor_corr[3] L_Omega_cs;
  vector<lower=0>[N_pos_cs] z_pos_cs;
  vector<upper=0>[N_neg_cs] z_neg_cs;
  vector[N_miss_cs] z_miss_cs;
  
  // Impute motor palsy Td[4:6] for components variables ----------------------
  // Xd[3] is the combined one aka clin_motor_palsy
  vector[3] mp_a0;
  matrix[3, 2] mp_a; //For HIV & illness days
  cholesky_factor_corr[3] L_Omega_mp;
  vector<lower=0>[N_pos_mp] z_pos_mp;
  vector<upper=0>[N_neg_mp] z_neg_mp;
  vector[N_miss_mp] z_miss_mp;
  
  // Impute age Xc[1] ---------------------------------------------------------
  real age_a0;
  real age_a; //for HIV
  real<lower=0> age_sigma;
  real age_imp[N - sum(obs_Xc[:,1])];
  
  // Impute illness day Xc[2]
  real id_a0;
  real id_a; //for HIV
  real<lower=0> id_sigma;
  real id_imp[N - sum(obs_Xc[:,2])];
  
  // Impute csf lab tests as mvNormal -----------------------------------------
  vector[5] csf_a0;
  vector[2] glu_a; //For diabetes
  cholesky_factor_corr[5] L_Omega_csf;
  vector<lower=0>[5] L_sigma_csf;
  real<lower=1> bld_glu_imp[N - obs_bld_glu]; //lower = 1 for numerical stability
  real<lower=0> csf_glu_imp[N - obs_csf_glu]; 
  real csf_other_imp[(N*3) - obs_csf_other];
  
  // Parameters of the logistics regression -----------------------------------
  real<lower=3> nu;
  real a0; //intercept
  // vector[nX + 1] a; //slope. 1 is for sqrt(glu_ratio) = sqrt(csf_glu)/sqrt(bld_glu)
  real<lower=0> a_pos; // assuming HIV must have positive effect. 
  vector[nX + 1 - 1] a_; // Extra 1 is for sqrt(glu_ratio) = sqrt(csf_glu)/sqrt(bld_glu)
  real b_HIV; //adjustment of RE with HIV X[,1]
  real<lower=0> b_RE;
  vector[N] RE; //base random effect;
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
}

transformed parameters {
  vector[nX + 1] a = append_row(a_pos, a_); //Add HIV coef to the vector of coef
  
  vector[3] z_cs[N]; // - Clinical symptoms
  vector[3] z_mp[N]; // - Motor palsy
  
  for (n in 1:N_miss_cs) z_cs[n_miss_cs[n], d_miss_cs[n]] = z_miss_cs[n];
  for (n in 1:N_pos_cs)  z_cs[n_pos_cs[n] , d_pos_cs[n]]  = z_pos_cs[n];
  for (n in 1:N_neg_cs)  z_cs[n_neg_cs[n] , d_neg_cs[n]]  = z_neg_cs[n];
  
  for (n in 1:N_miss_mp) z_mp[n_miss_mp[n], d_miss_mp[n]] = z_miss_mp[n];
  for (n in 1:N_pos_mp)  z_mp[n_pos_mp[n] , d_pos_mp[n]]  = z_pos_mp[n];
  for (n in 1:N_neg_mp)  z_mp[n_neg_mp[n] , d_neg_mp[n]]  = z_neg_mp[n];
}

model {
  real p_HIV = Phi(HIV_a0); //Populational Probability of having HIV
  matrix[N, 3] Xd_imp; //fully imputed discrete X
  vector[nXc + 1] Xc_imp[N]; //fully imputed cont X
  matrix[N, nX + 1 - 3] X_compl;
  Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], rep_array(Phi(HIV_a0), N)); //HIV
  Xd_imp[:,2] = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
  Xd_imp[:,3] = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
  Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
  Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
  Xc_imp[:,3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp)); //CSF Lab tests
  Xc_imp[:,8] = to_array_1d(to_vector(Xc_imp[:,4]) ./ to_vector(Xc_imp[:,3])); //Glucose ratio
  if (nXc > 7) // Other if exists
    for (j in 7:nXc) Xc_imp[:, j + 1] = Xc[:, j];
  
  X_compl = append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp)); //complement part of X
  
  
  nu ~ gamma(2,0.1);
  // Imputation ---------------------------------------------------------------
  // - HIV
  HIV_a0 ~ student_t(nu, 0, 1);
  // - Clinical symptoms
  L_Omega_cs ~ lkj_corr_cholesky(4);
  cs_a0 ~ student_t(nu, 0, 1);
  to_vector(cs_a) ~ student_t(5, 0, 1);
  // - Motor palsy
  L_Omega_mp ~ lkj_corr_cholesky(4);
  mp_a0 ~ student_t(nu, 0, 1);
  to_vector(mp_a) ~ student_t(nu, 0, 1);
  // - Age
  age_a0 ~ student_t(nu, 0, 1);
  age_a  ~ student_t(nu, 0, 1);
  age_sigma ~ normal(0, 2.5);
  
  // - Illness day
  id_a0 ~ student_t(nu, 0, 1);
  id_a  ~ student_t(nu, 0, 1);
  id_sigma ~ normal(0, 2.5);
  
  // - CSF lab tests
  to_vector(glu_a) ~ student_t(nu, 0, 1);
  L_Omega_csf ~ lkj_corr_cholesky(4);
  csf_a0 ~  student_t(nu, 0, 1);
  
  {
    vector[5] csf_mu[N];
    matrix[5, 5] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
    for (n in 1:N){
      csf_mu[n, 1:2] = csf_a0[1:2] + glu_a*Td[n, 7];
      csf_mu[n, 3:5] = csf_a0[3:5];
    }
    
    Xc_imp[:,3:7] ~ multi_normal_cholesky(csf_mu, L_Sigma_csf);
  }
  
  // Main model ---------------------------------------------------------------
  
  //Probs of each test become positive
  // Priors of covariates
  a0       ~ student_t(nu, 0, 1);
  a        ~ student_t(nu, 0, 1);
  
  //Random effects covariates
  RE    ~    normal(   0, 1);
  b_RE  ~ student_t(nu, 0, 1);
  b_HIV ~ student_t(nu, 0, 1);
  
  //1-Specificity of each test
  z_Xpert[1] ~ normal(logit(.005), 1.59);
  z_Mgit[1]  ~ normal(logit(.001),  .7 ); //.82
  z_Smear[1] ~ normal(logit(.001),  .7 ); //.7
  
  //Sensitivity of each test
  z_Xpert[2] ~ normal(0, .6);
  z_Mgit[2]  ~ normal(0, .6);
  z_Smear[2] ~ normal(0, .6);
  
  for (n in 1:N){
    int N_Xd_miss = 3 - sum(obs_Xd[n, 1:3]);

    if (obs_Xd[n,1]){
      Xd[n,1] ~ bernoulli(p_HIV);
      z_mp[n] ~ multi_normal_cholesky(mp_a0 + mp_a[:,1]*Xd[n,1] + mp_a[:,2]*Xc[n,2], L_Omega_mp);
      z_cs[n] ~ multi_normal_cholesky(cs_a0 + cs_a[:,1]*Xd[n,1] + cs_a[:,2]*Xc[n,2], L_Omega_cs);
      Xc_imp[n,1] ~ normal(age_a0 + age_a*Xd[n,1], age_sigma);
      Xc_imp[n,2] ~ normal(id_a0  + id_a*Xd[n,1] , id_sigma );
    } else {
      if (is_nan(
      multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,1] + cs_a[:,2]*Xc[n,2], L_Omega_cs) +
      multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,2]*Xc[n,2], L_Omega_cs) +
      multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,1] + mp_a[:,2]*Xc[n,2], L_Omega_mp) +
      multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + cs_a[:,2]*Xc[n,2], L_Omega_mp)))
      // This is to suppress the program from complaining at the start
      target += not_a_number();
      else {
        target += log_mix(p_HIV,
                          multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,1] + cs_a[:,2]*Xc[n,2], L_Omega_cs),
                          multi_normal_cholesky_lpdf(z_cs[n] | cs_a0             + cs_a[:,2]*Xc[n,2], L_Omega_cs));
        target += log_mix(p_HIV,
                          multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,1] + mp_a[:,2]*Xc[n,2], L_Omega_mp),
                          multi_normal_cholesky_lpdf(z_mp[n] | mp_a0             + mp_a[:,2]*Xc[n,2], L_Omega_mp));
        target += log_mix(p_HIV,
                          normal_lpdf(Xc_imp[n,1]| age_a0 + age_a, age_sigma),
                          normal_lpdf(Xc_imp[n,1]| age_a0        , age_sigma));    
        target += log_mix(p_HIV,
                          normal_lpdf(Xc_imp[n,2]| id_a0  + id_a, id_sigma ),
                          normal_lpdf(Xc_imp[n,2]| id_a0        , id_sigma ));                  
      }
    }
  
    if (N_Xd_miss > 0){ //if there is some discrete variables missing
      int N_pattern = int_power(2, N_Xd_miss);
      vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,:], obs_Xd[n, 1:3], a[1:3]);
      vector[N_pattern] log_liks;
      pat_thetas[2] += a0 + dot_product(a[4:(nX + 1)], X_compl[n]);
    
      //check if HIV is missing
      if (obs_Xd[n,1]){
        for (i in 1:N_pattern){
          real logprob_theta = pat_thetas[1,i];
          real theta = inv_logit(pat_thetas[2,i]);
          
          real bac_load   = b_HIV*Xd_imp[n,1];
          real z_Smear_RE = z_Smear[2] + bac_load + b_RE*RE[n];
          real z_Mgit_RE  = z_Mgit[2]  + bac_load + b_RE*RE[n];
          real z_Xpert_RE = z_Xpert[2] + bac_load + b_RE*RE[n];
          
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
          
          vector[2] z_Smear_RE = z_Smear[2] + bac_load + b_RE*RE[n];
          vector[2] z_Mgit_RE  = z_Mgit[2]  + bac_load + b_RE*RE[n];
          vector[2] z_Xpert_RE = z_Xpert[2] + bac_load + b_RE*RE[n];
          
          log_liks[i] = logprob_theta + log_mix(theta, log_sum_exp(
            logprob_Y[1] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1])),
            logprob_Y[2] + (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2]))
            ),
            bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      } 
      
      target += log_sum_exp(log_liks);
    } else {
      // The normal way
      row_vector[nX + 1] X = append_col(Xd_imp[n,:], X_compl[n,:]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load   = b_HIV*Xd_imp[n, 1];
      real z_Smear_RE = z_Smear[2] + bac_load + b_RE * RE[n];
      real z_Mgit_RE  = z_Mgit [2] + bac_load + b_RE * RE[n];
      real z_Xpert_RE = z_Xpert[2] + bac_load + b_RE * RE[n];
      
      target += log_mix(theta, 
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
    }
  }
}

generated quantities{
  vector<lower=0, upper=1>[N] theta;
  int <lower=0, upper=1> C[N];
  vector[N] log_lik;
  matrix[N, nX + 1] X;
  {
    matrix[N, 3] Xd_imp; //fully imputed discrete X
    vector[nXc + 1] Xc_imp[N]; //fully imputed cont X
    matrix[N, 3] Xd_rng;
    
    Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], rep_array(Phi(HIV_a0), N)); //HIV
    Xd_imp[:,2] = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
    Xd_imp[:,3] = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
    Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
    Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
    Xc_imp[:,3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp)); //CSF Lab tests
    Xc_imp[:,8] = to_array_1d(to_vector(Xc_imp[:,4]) ./ to_vector(Xc_imp[:,3])); //Glucose ratio
    if (nXc > 7) // Other if exists
    for (j in 7:nXc) Xc_imp[:, j + 1] = Xc[:, j];
    
    Xd_rng = binary_2d_rng(Xd_imp, obs_Xd[,1:3]);
    X = append_col(append_col(Xd_rng, to_matrix(Xd[:,4:nXd])), append_all(Xc_imp));
  }
 
  theta = inv_logit(a0 + X*a);
  C = bernoulli_rng(theta);

  for (n in 1:N){
    if (C[n] == 0){
      log_lik[n] = bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]);
    } else {
      real bac_load = b_HIV * X[n,1] + b_RE * RE[n];
      log_lik[n] = bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2] + bac_load) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2] + bac_load) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2] + bac_load);
    }
  }
}
