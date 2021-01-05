functions{
#include /includes/functions.stan
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

transformed data {
  //Train data
  int<lower=1> N = sum(keptin);
  int<lower=1> N_valid = N_all - N;
  int N_miss_gscv_valid = N_valid - sum(obs_Tc_all[which_not(keptin),3]);
  
  int<lower=0, upper=1> Y_Smear[N] = Y_Smear_all[which(keptin)];
  int<lower=0, upper=1> Y_Mgit[N] = Y_Mgit_all[which(keptin)];
  int<lower=0, upper=1> Y_Xpert[N] = Y_Xpert_all[which(keptin)];
  
  real Xc[N, nXc] = Xc_all[which(keptin),:];
  int  Xd[N, nXd] = Xd_all[which(keptin),:];
  real Tc[N, nTc] = Tc_all[which(keptin),:];
  int  Td[N, nTd] = Td_all[which(keptin),:];
  
  int<lower=0, upper=1> obs_Xc[N, nXc] = obs_Xc_all[which(keptin),:];
  int<lower=0, upper=1> obs_Xd[N, nXd] = obs_Xd_all[which(keptin),:];
  int<lower=0, upper=1> obs_Tc[N, nTc] = obs_Tc_all[which(keptin),:];
  int<lower=0, upper=1> obs_Td[N, nTd] = obs_Td_all[which(keptin),:];
  
#include /includes/transform_data.stan
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
  
  // Impute gcs Xc[8]
  // Use auxillary gcs components: Tc[1:3] GCSE, GCSM, GCSV
  // GCSV is missing.
  real<lower=0> gcsv_a0;
  vector[2] gcsv_a;
  real<lower=0> gcsv_sigma;
  real<lower=0, upper=4> gcsv_imp[N - sum(obs_Tc[,3])];
  // --------------------------------------------------------------------------
  // Parameters of the logistics regression -----------------------------------
  // real<lower=3> nu;
  real a0; //intercept
  real<lower=0> a_pos; // assuming HIV must have positive effect. 
  vector[nX + 1 - 1] a_; // Extra 1 is for sqrt(glu_ratio) = sqrt(csf_glu)/sqrt(bld_glu)
  real b_HIV; //adjustment of RE with HIV X[,1]
  real b_glu;
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
  real GCSV_imp[N]; 
  int nu = 4;
 
  Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], rep_array(HIV_a0, N)); //HIV
  Xd_imp[:,2] = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
  Xd_imp[:,3] = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
  Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
  Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
  Xc_imp[:,3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp)); //CSF Lab tests
  Xc_imp[:,8] = to_array_1d(to_vector(Xc_imp[:,4]) ./ to_vector(Xc_imp[:,3])); //Glucose ratio
  // - GCSV
  GCSV_imp = impute_cont_1d(Tc[:,3], obs_Tc[:,3], gcsv_imp); //GCSV
  Xc_imp[:,9] = to_array_1d(log2(to_vector(Tc[:,1]) + to_vector(Tc[:,2]) + to_vector(GCSV_imp) + 1));

 
  if (nXc > 8) // Other if exists
    for (j in 9:nXc) Xc_imp[:, j + 1] = Xc[:, j];
  
  X_compl = append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp)); //complement part of X
  
  // nu ~ gamma(2, .9);
  // Imputation ---------------------------------------------------------------
  // - HIV
  HIV_a0 ~ student_t(nu, 0, 1);
  // - Clinical symptoms
  L_Omega_cs ~ lkj_corr_cholesky(4);
  cs_a0 ~ student_t(nu, 0, 1);
  to_vector(cs_a) ~ student_t(nu, 0, 1);
  // - Motor palsy
  L_Omega_mp ~ lkj_corr_cholesky(4);
  mp_a0 ~ student_t(nu, 0, 1);
  to_vector(mp_a) ~ student_t(nu, 0, 1);
  // - Age
  age_a0 ~ student_t(nu, 0, 1);
  age_a  ~ student_t(nu, 0, 1);
  age_sigma ~ normal(0, 1);
  
  // - Illness day
  id_a0 ~ student_t(nu, 0, 1);
  id_a  ~ student_t(nu, 0, 1);
  id_sigma ~ normal(0, 1);
  
  // - CSF lab tests  // -GCSV
  to_vector(glu_a) ~ student_t(nu, 0, 1);
  L_Omega_csf ~ lkj_corr_cholesky(4);
  L_sigma_csf ~ student_t(nu, 0, 1);
  csf_a0 ~  student_t(nu, 0, 1);
 
  gcsv_a0 ~ student_t(nu, 0, 1);
  gcsv_a  ~ student_t(nu, 0, 1);
  gcsv_sigma ~ normal(0, 1);
  
  {
    vector[5] csf_mu[N];
    matrix[5, 5] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
    real gcsv_x[N];
    for (n in 1:N){
      csf_mu[n, 1:2] = csf_a0[1:2] + glu_a*Td[n, 7];
      csf_mu[n, 3:5] = csf_a0[3:5];
      gcsv_x[n] = gcsv_a0 + dot_product(gcsv_a,(to_vector(Tc[n, 1:2])));
    }
    
    Xc_imp[:,3:7] ~ multi_normal_cholesky(csf_mu, L_Sigma_csf);
    GCSV_imp ~ normal(gcsv_x, gcsv_sigma);
  }
  
  // Main model ---------------------------------------------------------------
  
  //Probs of each test become positive
#include /includes/main_priors.stan
b_glu ~ student_t(nu, 0, 1);
  
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
                          normal_lpdf(Xc_imp[n,2]| id_a0  + id_a , id_sigma ),
                          normal_lpdf(Xc_imp[n,2]| id_a0         , id_sigma ));                  
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
          
          real bac_load   = b_HIV*Xd_imp[n,1] + b_glu*Xc_imp[n,4];
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
          vector[2] bac_load = pat_bac_load[2] + b_RE*RE[n] + b_glu*Xc_imp[n,4];
          
          vector[2] z_Smear_RE = z_Smear[2] + bac_load;
          vector[2] z_Mgit_RE  = z_Mgit[2]  + bac_load;
          vector[2] z_Xpert_RE = z_Xpert[2] + bac_load;
          
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
      real z_Smear_RE = z_Smear[2] + bac_load + b_RE * RE[n] + b_glu*Xc_imp[n,4];
      real z_Mgit_RE  = z_Mgit [2] + bac_load + b_RE * RE[n] + b_glu*Xc_imp[n,4];
      real z_Xpert_RE = z_Xpert[2] + bac_load + b_RE * RE[n] + b_glu*Xc_imp[n,4];
      
      target += log_mix(theta, 
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
      bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
    }
  }
}

generated quantities{
  vector[N_all] log_lik;
  vector[N_all] p_Smear;
  vector[N_all] p_Mgit;
  vector[N_all] p_Xpert;
  vector[N_all] theta;
  
  {
    int  C[N_all];
    matrix[N_all, nX + 1] X;
    
    {
      // Imputation model --------------------------------------------------------
      matrix[N_all, 3] Xd_imp; //fully imputed discrete X
      vector[nXc + 1] Xc_imp[N_all]; //fully imputed cont X
      // matrix[N_all, 3] Xd_rng;
      real age_imp_all[N_all - sum(obs_Xc_all[:,1])];
      real id_imp_all[N_all - sum(obs_Xc_all[:,2])];
      // print(which_not(keptin));
      
      Xd_imp[:,1] = binary_rng(impute_binary(Xd_all[:,1], obs_Xd_all[:,1], rep_array(HIV_a0, N_all)), obs_Xd_all[:,1]); //HIV
      age_imp_all = to_array_1d(normal_rng(age_a0 + age_a*to_vector(Xd_imp[which_not(obs_Xc_all[:,1]),1]), age_sigma));
      id_imp_all  = to_array_1d(normal_rng(id_a0  + id_a*to_vector(Xd_imp[which_not(obs_Xc_all[:,2]),1]) , id_sigma ));
      Xc_imp[:,1] = impute_cont_1d(Xc_all[:,1], obs_Xc_all[:,1], age_imp_all); //Age
      Xc_imp[:,2] = impute_cont_1d(Xc_all[:,2], obs_Xc_all[:,2], id_imp_all); //Illness day
      
      {
        int Td_cs_valid[N_valid,3] = Td_all[which_not(keptin),1:3];
        int obs_cs_valid[N_valid,3] = obs_Td_all[which_not(keptin),1:3];
        vector[3] Mu_cs[N_valid];
        int CS_imp_valid[N_valid, 3];
        int j = 1;
        for (n in which_not(keptin)){
          Mu_cs[j] = cs_a0 + cs_a[:,1]*Xd_all[n,1] + cs_a[:,2]*Xc_all[n,2];
          j += 1;
        }
        CS_imp_valid = multi_probit_partial_rng(Td_cs_valid, obs_cs_valid, Mu_cs, L_Omega_cs);
        Xd_imp[which(keptin),2] = binary_rng(impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_Td[:,1:3]), obs_Xd[:,2]); //Clinical Symptoms
        Xd_imp[which_not(keptin),2] = to_vector(any(CS_imp_valid));
      }
      
      
      {
        int Td_mp_valid[N_valid,3] = Td_all[which_not(keptin),4:6];
        int obs_mp_valid[N_valid,3] = obs_Td_all[which_not(keptin),4:6];
        vector[3] Mu_mp[N_valid];
        int mp_imp_valid[N_valid, 3];
        int j = 1;
        for (n in which_not(keptin)){
          Mu_mp[j] = mp_a0 + mp_a[:,1]*Xd_all[n,1] + mp_a[:,2]*Xc_all[n,2];
          j += 1;
        }
        mp_imp_valid = multi_probit_partial_rng(Td_mp_valid, obs_mp_valid, Mu_mp, L_Omega_mp);
        Xd_imp[which(keptin),3] = binary_rng(impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_mp)), obs_Td[:,1:3]), obs_Xd[:,2]); //Clinical Symptoms
        Xd_imp[which_not(keptin),3] = to_vector(any(mp_imp_valid));
      }
      
      {
        vector[5] csf_mu_valid[N_valid];
        int j = 1;
        for (n in which_not(keptin)){
          csf_mu_valid[j, 1:2] = csf_a0[1:2] + glu_a*Td_all[n, 7];
          csf_mu_valid[j, 3:5] = csf_a0[3:5];
          j += 1;
        }
        Xc_imp[which_not(keptin),3:7] = multi_normal_cholesky_partial_rng(Xc_all[which_not(keptin),3:7], obs_Xc_all[which_not(keptin),3:7], csf_mu_valid, L_Omega_csf);
        Xc_imp[which(keptin),3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp));
      }
      
      Xc_imp[:,8] = to_array_1d(to_vector(Xc_imp[:,4]) ./ to_vector(Xc_imp[:,3])); //Glucose ratio
      
      {
        real GCSV_imp_all[N_all]; // - GCSV
        real gcsv_imp_valid[N_miss_gscv_valid];
        if (N_miss_gscv_valid > 0){
          real gcsv_x_valid[N_miss_gscv_valid];
          int j = 1;
          for (n in which_not(obs_Tc_all[which_not(keptin),3])){
            gcsv_x_valid[j] = gcsv_a0 + dot_product(gcsv_a,(to_vector(Tc_all[n, 1:2])));
            j += 1;
          }
          
          gcsv_imp_valid = to_array_1d(half_normal_rng(gcsv_x_valid, gcsv_sigma));
          // print(gcsv_imp_valid);
        }
        
        GCSV_imp_all[which_not(keptin)] = impute_cont_1d(Tc_all[which_not(keptin),3], obs_Tc_all[which_not(keptin),3], gcsv_imp_valid); 
        GCSV_imp_all[which(keptin)] = impute_cont_1d(Tc[:,3], obs_Tc[:,3], gcsv_imp);
        Xc_imp[:,9] = to_array_1d(log2(to_vector(Tc_all[:,1]) + to_vector(Tc_all[:,2]) + to_vector(GCSV_imp_all) + 1));
      }
      
      if (nXc > 8) // Other if exists
      for (j in 9:nXc) Xc_imp[:, j + 1] = Xc_all[:, j];
      
      // Xd_rng = binary_2d_rng(Xd_imp, obs_Xd_all[,1:3]);
      X = append_col(append_col(Xd_imp, to_matrix(Xd_all[:,4:nXd])), append_all(Xc_imp));
    }
    
    theta = inv_logit(a0 + X*a);
    C = bernoulli_rng(theta);
    
    {
      vector[N_all] RE_all;
      vector[N_all] bac_load;
      RE_all[which(keptin)] = RE;
      for (n in which_not(keptin)) RE_all[n] = normal_rng(0, 1);
      bac_load = b_RE*RE_all + b_HIV*X[:,1] + b_glu*X[:,10];
      p_Smear = (1 - theta) * inv_logit(z_Smear[1]) + theta .* inv_logit(z_Smear[2] + bac_load);
      p_Mgit  = (1 - theta) * inv_logit(z_Mgit[1])  + theta .* inv_logit(z_Mgit[2]  + bac_load);
      p_Xpert = (1 - theta) * inv_logit(z_Xpert[1]) + theta .* inv_logit(z_Xpert[2] + bac_load);
      for (n in 1:N_all){
        if (C[n] == 0){
          log_lik[n] = bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[1]);
        } else {
          log_lik[n] = bernoulli_logit_lpmf(Y_Xpert_all[n] | z_Xpert[2] + bac_load[n]) + bernoulli_logit_lpmf(Y_Mgit_all[n] | z_Mgit[2] + bac_load[n]) + bernoulli_logit_lpmf(Y_Smear_all[n] | z_Smear[2] + bac_load[n]);
        }
      }
    }
  }
}
