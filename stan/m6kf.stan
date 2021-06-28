functions{
#include includes/functions.stan
#include includes/pca_helpers.stan
}

data {
  int<lower=1> N_all;   //Number of patient
  int<lower=1, upper=4> nFA; //Latent factors
  
#include includes/data/X.stan
#include includes/data/Y.stan
#include includes/cross_validation/data.stan
#include includes/data/penalty.stan  
}

transformed data{
  // Penalty term adaptation
  int adapt_penalty[2];
  // CSF PCA
  int P_csf = 5;
  int D_csf = nFA;
  int m_csf = D_csf*(P_csf-D_csf) + (D_csf*(D_csf-1))%/%2;
  // * Global variables -------------------------------------------------------
  int nX = nXc + nXd; // Total number of covariatesi
  // int mFA = 
  real CSF_scaled[N_all, 5];
  // vector[5] mu_csf = rep_vector(0, 5);
#include includes/cross_validation/transform_data_Y.stan
#include includes/cross_validation/transform_data_X.stan
#include includes/impute_model/transform_data.stan

  for (i in 1:2) adapt_penalty[i] = penalty_term[i] == 0 ? 1 : 0;
  for (i in 1:5) {
    vector[sum(obs_Xc_all[:,i+2])] X_unscaled = to_vector(Xc_all[which(obs_Xc_all[:,i+2]),i+2]);
    CSF_scaled[which(obs_Xc_all[:,i+2]), i] = to_array_1d((X_unscaled - mean(X_unscaled)) ./ sd(X_unscaled));
    CSF_scaled[which_not(obs_Xc_all[:,i+2]), i] = rep_array(0.0, N_all - sum(obs_Xc_all[:,i+2]));
  }
}
 
parameters {
    // Parameters of the imputation model ---------------------------------------
// * Impute HIV Xd[1] -------------------------------------------------------
real HIV_a0;
vector[2] HIV_a; //WHITE and LYM

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
// vector[2] glu_a; //For diabetes
real bld_glu_imp[N - obs_bld_glu];
real csf_glu_imp[N - obs_csf_glu]; 
real csf_other_imp[(N*3) - obs_csf_other];

  
  // Parameters of the logistics regression -----------------------------------
  real a0; //intercept
  real<lower=0> a_pos; // assuming HIV must have positive effect. 
  vector[nX - 1 - 5] a_; // Extra 1 is for sqrt(glu_ratio) = sqrt(csf_glu)/sqrt(bld_glu)
  vector[nFA] a_FA; //contraints the order of factor analysis #or ordered
  real b_HIV; //adjustment of RE with HIV Xd[,1]
  vector<lower=0>[3] b_RE;
  vector[N] RE; //base random effect;
  
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
 
  //Penalty terms
  real<lower=0> sp1[adapt_penalty[1]]; //sigma for the penalty
  real<lower=0> sp2[adapt_penalty[2]]; //sigma for the penalty
  
  vector[m_csf] L_t0_csf;   // lower diagonal elements of L
  vector<lower=0>[D_csf] L_d_csf;   // lower diagonal elements of L
  real<lower=0> psi0_csf;         // vector of variances
  real<lower=0>   mu_psi_csf;
  real<lower=0>  sigma_psi_csf;
  real mu_lt_csf;
  real<lower=0>  sigma_lt_csf; 
  vector[D_csf] U0_csf[N];
}


transformed parameters {
  vector[nX - (5-nFA)] a;
#include includes/impute_model/transform_parameters.stan
  // matrix[nFA,nFA] L_Sigma_Ucsf = diag_pre_multiply(L_sigma_Ucsf, L_Omega_Ucsf);
  if (nXc > 7) a = append_row(append_row(append_row(a_pos, a_[1:(nXd - 1 + 2)]), a_FA), a_[(nXd - 1 + 3):]);
  else a = append_row(append_row(a_pos, a_[1:(nXd - 1 + 2)]), a_FA);
  vector[m_csf] L_t_csf = L_t0_csf * sigma_lt_csf + mu_lt_csf;
  real<lower=0> psi_csf = psi0_csf * sigma_psi_csf + mu_psi_csf; 
  cholesky_factor_cov[P_csf,D_csf] L_csf;  //lower triangular factor loadings Matrix 
  cov_matrix[P_csf] Q_csf;   //Covariance mat
  
  {
    int idx2 = 0;
    for(i in 1:(D_csf-1)){
      for(j in (i+1):D_csf){
        L_csf[i,j] = 0; //constrain the upper triangular elements to zero 
      }
    }
    for (j in 1:D_csf) {
        L_csf[j,j] = L_d_csf[j];
      for (i in (j+1):P_csf) {
        idx2 += 1;
        L_csf[i,j] = L_t_csf[idx2];
      } 
    }
  } 
  Q_csf = L_csf*L_csf' + diag_matrix(rep_vector(psi_csf + 1e-14, P_csf)); 
}


model {
  int nu = 4;
  real SP[2] = {adapt_penalty[1] == 1 ? sp1[1] : penalty_term[1], adapt_penalty[2] == 1 ? sp2[1] : penalty_term[2]};
   // Imputation model ---------------------------------------------------------
  vector[N] z_HIV = HIV_a0 + to_matrix(Tc[:,4:5]) * HIV_a; //Populational Probability of having HIVi
  vector[N] p_HIV = Phi(z_HIV);
  matrix[N, 3] Xd_imp; //fully imputed discrete X
  vector[nXc] Xc_imp[N]; //fully imputed cont X
  matrix[N, nX - 3 - (5-nFA)] X_compl;

  Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], to_array_1d(z_HIV)); //HIV
  Xd_imp[:,2] = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
  Xd_imp[:,3] = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
  Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
  Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
  Xc_imp[:,3:7] = impute_cont_2d(CSF_scaled[which(keptin),:], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp)); //CSF Lab tests
  
  if (nXc > 7) // Other if exists
    for (j in 8:nXc) Xc_imp[:, j] = Xc[:, j];
  
 
  // Imputation ---------------------------------------------------------------
  // - HIV
    HIV_a0 ~ normal(0, 2.5);
    HIV_a ~ normal(0, 2.5);
    // - Clinical symptoms
    L_Omega_cs ~ lkj_corr_cholesky(4);
    cs_a0 ~ normal(0, 2.5);
    to_vector(cs_a) ~ normal(0, 2.5);
    // - Motor palsy
    L_Omega_mp ~ lkj_corr_cholesky(4);
    mp_a0 ~ normal(0, 2.5);
    to_vector(mp_a) ~ normal(0, 2.5);
    // - Age
    age_a0 ~ normal(0, 2.5);
    age_a  ~ normal(0, 2.5);
    age_sigma ~ std_normal();
    
    // - Illness day
    id_a0 ~ normal(0, 2.5);
    id_a  ~ normal(0, 2.5);
    id_sigma ~ std_normal();

    // - CSF lab with pca
    mu_psi_csf ~ std_normal();
    sigma_psi_csf ~ std_normal();
    mu_lt_csf ~ std_normal();
    sigma_lt_csf ~ std_normal();
    L_t0_csf ~ std_normal();
    psi0_csf ~ std_normal();
    L_d_csf ~ normal(0,1);
    
    {
      matrix[P_csf, P_csf] L = cholesky_decompose(Q_csf);
      matrix[D_csf, D_csf] G = cholesky_decompose(diag_matrix(rep_vector(1, D_csf)) - L_csf'/Q_csf*L_csf); 
      // vector[D_csf] Mu_Ucsf[N];
      vector[D_csf] U_csf[N];
      for (n in 1:N){
        U0_csf[n] ~ std_normal();
        vector[D_csf] mu_Ucsf = to_vector(L_csf'/Q_csf * to_matrix(Xc_imp[n,3:7]));
        U_csf[n] = G * U0_csf[n] + mu_Ucsf;
      }
      
      Xc_imp[:,3:7] ~ multi_normal_cholesky(rep_vector(0, P_csf), L);
   
      if (nXc > 7)
        X_compl = append_col(append_col(append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp[:,1:2])), append_all(U_csf)), append_all(Xc_imp[:,8:])); //complement part of X
      else
        X_compl = append_col(append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp[:,1:2])), append_all(U_csf)); //complement part of X
    }
  
  // Main model ---------------------------------------------------------------
#include includes/main_prior/m0.stan
#include includes/main_prior/m.stan
#include includes/main_prior/m_RE.stan
#include includes/main_prior/penalty.stan

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
          
          vector[2] z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n]);
          vector[2] z_Mgit_RE  = z_Mgit[2]  + b_RE[2]*(bac_load + RE[n]);
          vector[2] z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n]);
          
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
      row_vector[nX - (5-nFA)] X = append_col(Xd_imp[n,:], X_compl[n,:]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load   = b_HIV*Xd_imp[n, 1];
      real z_Smear_RE = z_Smear[2] + b_RE[1]*(bac_load + RE[n]);
      real z_Mgit_RE  = z_Mgit [2] + b_RE[2]*(bac_load + RE[n]);
      real z_Xpert_RE = z_Xpert[2] + b_RE[3]*(bac_load + RE[n]);
      
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
  matrix[N_all, nX - (5 - nFA)] X;
  vector[nFA] U_csf_all[N_all];

  {
    {
    // Imputation model --------------------------------------------------------
      matrix[N_all, 3] Xd_imp; //fully imputed discrete X
      vector[nXc] Xc_imp[N_all]; //fully imputed cont X
      real age_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),1])];
      real id_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),2])];
      
      {
        vector[N_all] z_HIV = HIV_a0 + to_matrix(Tc_all[:,4:5]) * HIV_a;
        Xd_imp[:,1] = binary_rng(impute_binary(Xd_all[:,1], obs_Xd_all[:,1], to_array_1d(z_HIV)), obs_Xd_all[:,1]);
      }
      
      //Age & illness day
      {
        int j =1; int k= 1;
        for (n in which_not(keptin)){
          if (obs_Xc_all[n,1] == 0) {
            age_imp_valid[j]  = normal_rng(age_a0  + age_a*(Xd_imp[n,1]) , age_sigma );
            j+=1;
          }
          if (obs_Xc_all[n,2] == 0) {
            id_imp_valid[k]  = normal_rng(id_a0  + id_a*(Xd_imp[n,1]) , id_sigma );
            k+=1;
          }
        }
      }
      
      Xc_imp[which_not(keptin),1] = impute_cont_1d(Xc_all[which_not(keptin),1], obs_Xc_all[which_not(keptin),1], age_imp_valid); //Age
      Xc_imp[which_not(keptin),2] = impute_cont_1d(Xc_all[which_not(keptin),2], obs_Xc_all[which_not(keptin),2], id_imp_valid); //illness day
      Xc_imp[which(keptin),1] = impute_cont_1d(Xc[:, 1], obs_Xc[:, 1], age_imp);
      Xc_imp[which(keptin),2] = impute_cont_1d(Xc[:, 2], obs_Xc[:, 2], id_imp);
      
      // Clinical symptoms
      {
        int Td_cs_valid[N_valid,3] = Td_all[which_not(keptin),1:3];
        int obs_cs_valid[N_valid,3] = obs_Td_all[which_not(keptin),1:3];
        vector[3] Mu_cs[N_valid];
        int CS_imp_valid[N_valid, 3];
        int j = 1;
        for (n in which_not(keptin)){
          Mu_cs[j] = cs_a0 + cs_a[:,1]*Xd_imp[n,1] + cs_a[:,2]*Xc_imp[n,2];
          j += 1;
        }
        CS_imp_valid = multi_probit_partial_rng(Td_cs_valid, obs_cs_valid, Mu_cs, L_Omega_cs);
        Xd_imp[which(keptin),2] = binary_rng(impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs), obs_Xd[:,2]); //Clinical Symptoms
        Xd_imp[which_not(keptin),2] = to_vector(any(CS_imp_valid));
      }
      
      // Motor palsy
      {
        int Td_mp_valid[N_valid,3] = Td_all[which_not(keptin),4:6];
        int obs_mp_valid[N_valid,3] = obs_Td_all[which_not(keptin),4:6];
        vector[3] Mu_mp[N_valid];
        int mp_imp_valid[N_valid, 3];
        int j = 1;
        for (n in which_not(keptin)){
          Mu_mp[j] = mp_a0 + mp_a[:,1]*Xd_imp[n,1] + mp_a[:,2]*Xc_imp[n,2];
          j += 1;
        }
        mp_imp_valid = multi_probit_partial_rng(Td_mp_valid, obs_mp_valid, Mu_mp, L_Omega_mp);
        Xd_imp[which(keptin),3] = binary_rng(impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp), obs_Xd[:,3]); //Motor palsy
        Xd_imp[which_not(keptin),3] = to_vector(any(mp_imp_valid));
      }
      
      if (nXc > 7) // Other if exists
      for (j in 8:nXc) Xc_imp[:, j] = Xc_all[:, j];
      
      {
        vector[P_csf] Mu_csf[N] = rep_array(rep_vector(0, P_csf), N);
        matrix[P_csf, P_csf] L = cholesky_decompose(Q_csf);
        matrix[D_csf, D_csf] G = cholesky_decompose(diag_matrix(rep_vector(1, D_csf)) - L_csf'/Q_csf*L_csf); 
        int j = 1;
        vector[D_csf] U_csf[N];
        Xc_imp[which(keptin),3:7] = impute_cont_2d(CSF_scaled[which(keptin),:], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp));
        Xc_imp[which_not(keptin), 3:7] = multi_normal_cholesky_partial_rng(CSF_scaled[which_not(keptin),:], obs_Xc_all[which_not(keptin),3:7], Mu_csf, L);
        // j = 1; //rest j
        for (n in which_not(keptin)) {
          U_csf_all[n, :] = multi_normal_cholesky_rng(to_vector(L_csf'/Q_csf * to_matrix(Xc_imp[n, 3:7])), G);
          j += 1;
        }
        j = 1; //reset j 
        for (n in which(keptin)){
          vector[D_csf] mu_Ucsf = to_vector(L_csf'/Q_csf * to_matrix(Xc_imp[n,3:7]));
          U_csf_all[n] = G * U0_csf[j] + mu_Ucsf;
          j += 1;
        }
        
        
        if (nXc > 7)
           X = append_col(append_col(append_col(append_col(Xd_imp, to_matrix(Xd_all[:,4:nXd])), append_all(Xc_imp[:,1:2])), append_all(U_csf_all)), append_all(Xc_imp[:,8:]));
        else
           X = append_col(append_col(append_col(Xd_imp, to_matrix(Xd_all[:,4:nXd])), append_all(Xc_imp[:,1:2])), append_all(U_csf_all));
        
      }
    }
    
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
