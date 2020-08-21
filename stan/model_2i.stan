functions {
  int sum2d(int[,] a) {
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }
  
  int sum2d_with_missing(int[,] a, int[,] b){
    int s = 0;
    if (size(a) != size(b)) reject("Size a and b must match!");
    for (i in 1:size(a)){
      for (j in 1:size(a[i]))
        s += (b[i,j] == 1) ? a[i,j] : 0;
    }
    return s;
  }
  
  //Append all vector in an array to a matrix
  matrix append_all(vector[] x){
    matrix[size(x), dims(x[1])[1]] M;
    for (i in 1:size(x)) M[i,] = x[i]'; 
    return M;
  }
  
  //Function to sum over 2|3 "or" logic.
  real sum_probs_or(real[] probs){
    int s = size(probs);
    real sum_prob;
    
    if (s>3||s<2) reject("Only size 2 and 3 support");
    sum_prob = (s==2) ? sum(probs) - prod(probs) : sum(probs) - prod(probs[1:2]) - prod(probs[2:3]) - prod(probs[{1,3}]) + prod(probs);
    
    return sum_prob;
  }
  
  // Impute discrete variable. Missing value will be "imputed" by their expected probs
  real[] impute_discrete(int[] raw, int[] obs, real[] z){
    int N = size(raw);
    real imp[N];
    
    if (size(obs)!=N||size(z)!=N) reject("Size mismatched!");
    
    for (n in 1:N)
      imp[n] = obs[n] ? raw[n] : Phi(z[n]);
    
    return imp;
  }
  
  // Impute discrete variable for combined val. Missing value will be "imputed" by their expected probs
  real[] impute_discrete(int[] raw_cmb, int[] obs_cmb, real[,] z_el, int[,] obs_el){
    int N = size(raw_cmb);
    real cmb[N];
    
    if (size(obs_cmb)!=N||size(z_el)!=N||size(obs_el)!=N) reject("Size mismatched!");
    
    for (n in 1:N){
      if (obs_cmb[n]) cmb[n] = raw_cmb[n];
      else {
        int n_obs_el = sum(obs_el[n]);
        if (n_obs_el == 0) cmb[n] = sum_probs_or(Phi(z_el[n]));
        if (n_obs_el == 1) {
          if (obs_el[n,1] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,2:3]));
          if (obs_el[n,2] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,{1,2}]));
          if (obs_el[n,3] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,1:2]));
        }
        if (n_obs_el == 2){
          if (obs_el[n,1] == 0) cmb[n] = Phi(z_el[n,1]);
          if (obs_el[n,2] == 0) cmb[n] = Phi(z_el[n,2]);
          if (obs_el[n,3] == 0) cmb[n] = Phi(z_el[n,3]);
        }
      }
    }
    return cmb;
  }
  
  // Impute continuous variable
  real[] impute_cont(int[] raw, int[] obs, real[] rep){
    int N = size(raw);
    real imp[N];
    int i = 1;
    
    if (size(obs)!=N) reject("Size mismatched!");
    
    for (n in 1:N){
      if (obs[n]) imp[n] = raw[n];
      else {
        imp[n] = rep[i];
        i += 1;
      }
    }
    return imp;
  }
}

data {
  int<lower=1> N; //Number of patient
  int<lower=1> nXc; //Number of continuous X
  int<lower=1> nXd; //Number of discrete X
  int<lower=1> nTd; //Number of disc. aux covariates for imputation model
  int<lower=1> nTc; //Number of cont. aux covariates for imputation model
  
  int<lower=0, upper=1> Y_Smear[N];
  int<lower=0, upper=1> Y_Mgit[N];
  int<lower=0, upper=1> Y_Xpert[N];
  
  real Xc[N, nXc]; //Continuous covariates
  int  Xd[N, nXd]; //Discrete covariates
  int  Td[N, nTd]; //Auxillary covariates - discrete
  int  Td[N, nTc]; //Auxillary covariates - cont
  
  //Matrices of observation
  int<lower=0, upper=1> obs_Xc[N, nXc]; 
  int<lower=0, upper=1> obs_Xd[N, nXd];
  int<lower=0, upper=1> obs_Td[N, nTd]; 
  int<lower=0, upper=1> obs_Tc[N, nTc]; 
}

transformed data {
  // Clinical symptoms Td[1:3] ------------------------------------------------
  int Td_cs[N,3] = Td[,1:3];
  int obs_cs[N,3] = obs_Td[,1:3];
  int<lower=0> N_miss_cs;
  int<lower=1,upper=N> n_miss_cs[(N * 3) - sum2d(obs_cs)];
  int<lower=1,upper=3> d_miss_cs[size(n_miss_cs)];
  int<lower=0> N_pos_cs;
  int<lower=1,upper=N> n_pos_cs[sum2d_with_missing(Td_cs, obs_cs)];
  int<lower=1,upper=3> d_pos_cs[size(n_pos_cs)];
  int<lower=0> N_neg_cs;
  int<lower=1,upper=N> n_neg_cs[sum2d(obs_cs) - size(n_pos_cs)];
  int<lower=1,upper=3> d_neg_cs[size(n_neg_cs)];

  N_pos_cs = size(n_pos_cs);
  N_neg_cs = size(n_neg_cs);
  N_miss_cs = size(n_miss_cs);
  {
    int i = 1;
    int j = 1;
    int k = 1;
    for (n in 1:N) {
      for (d in 1:3) {
        if (obs_cs[n,d] == 0){
          n_miss_cs[k] = n;
          d_miss_cs[k] = d;
          k += 1;
        } else if (Td_cs[n,d] == 1) {
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
  // ---------------------------------------------------------------------------
  
  
  // Motor palsy Td[4:6] ------------------------------------------------
  int Td_mp[N,3] = Td[,4:6];
  int obs_mp[N,3] = obs_Td[,4:6];
  int<lower=0> N_miss_mp;
  int<lower=1,upper=N> n_miss_mp[(N * 3) - sum2d(obs_mp)];
  int<lower=1,upper=3> d_miss_mp[size(n_miss_mp)];
  int<lower=0> N_pos_mp;
  int<lower=1,upper=N> n_pos_mp[sum2d_with_missing(Td_mp, obs_mp)];
  int<lower=1,upper=3> d_pos_mp[size(n_pos_mp)];
  int<lower=0> N_neg_mp;
  int<lower=1,upper=N> n_neg_cs[sum2d(obs_mp) - size(n_pos_mp)];
  int<lower=1,upper=3> d_neg_cs[size(n_neg_mp)];

  N_pos_mp = size(n_pos_mp);
  N_neg_mp = size(n_neg_mp);
  N_miss_mp = size(n_miss_mp);
  {
    int i = 1;
    int j = 1;
    int k = 1;
    for (n in 1:N) {
      for (d in 1:3) {
        if (obs_mp[n,d] == 0){
          n_miss_mp[k] = n;
          d_miss_mp[k] = d;
          k += 1;
        } else if (Td_mp[n,d] == 1) {
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
  // ---------------------------------------------------------------------------
}
 
parameters {
  // For imputation ------------------------------------------------------------
  
  // Impute HIV Xd[1]
  real HIV_a0;
  
  // Impute symptoms Td[1:3] for components variables
  // Xd[2] is the combined one aka clin_symptoms
  vector[3] cs_a;
  cholesky_factor_corr[3] L_Omega_cs;
  vector<lower=0>[N_pos_cs] z_pos_cs;
  vector<upper=0>[N_neg_cs] z_neg_cs;
  vector[N_miss_cs] z_miss_cs;
  
  // Impute motor palsy Td[4:6] for components variables
  // Xd[3] is the combined one aka clin_symptoms
  vector[3] mp_a;
  cholesky_factor_corr[3] L_Omega_mp;
  vector<lower=0>[N_pos_mp] z_pos_mp;
  vector<upper=0>[N_neg_mp] z_neg_mp;
  vector[N_miss_mp] z_miss_mp;
  
  // Impute age Xc[1]
  real mu_age;
  real<lower=0> sigma_age;
  real age_imp[N - sum(obs_Xc[,1])];
  
  // Impute illness day Xc[2]
  real mu_id;
  real<lower=0> sigma_id;
  real id_imp[N - sum(obs_Xc[,2])];
  
  // Impute glucose ratio, blood glucose, csf lymphocytes count, protein, and lactate
  // gluratio: Xc[3], bldglu: Xc[4], loglym: Xc[5], logprotein: Xc[6], loglac: Xc[7]
  real 
  
  // For logit regression ------------------------------------------------------
  real a0; //intercept
  vector[nX] a; //slope
  real b_HIV; //adjustment of RE with HIV X[,1]
  real<lower=0> b[3];
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  
  vector[N] RE; //base random effect;
}

transformed parameters{
  // Create an fully imputed X matrix
  matrix[N, nXd] Xd_imp; //fully imputed discrete X
  matrix[N, nXc] Xc_imp; //fully imputed cont X
  
  int nX = nXc + nXd;
  matrix[N, nX] X;
  matrix[N, nX] obs_X';
  
  // Imputation models ---------------------------------------------------------
  // - Clinical symptoms
  vector[3] z_cs[N];
  
  // - Motor palsy
  vector[3] z_mp[N];
  
  // Main model ----------------------------------------------------------------
  
  //Prevalence of TBM
  real theta[N] = to_array_1d(inv_logit(a0 + X*a));
  
  // z of each test positivitiy with random effect
  matrix[N,2] z_Smear_RE; 
  matrix[N,2] z_Mgit_RE; 
  matrix[N,2] z_Xpert_RE;
  
  //Bacterial load
  vector[N] bac_load;
  bac_load = RE + b_HIV*X[,1];
  
  //Add random effects
  z_Smear_RE = rep_matrix(z_Smear' , N);
  z_Mgit_RE  = rep_matrix(z_Mgit' , N);
  z_Xpert_RE = rep_matrix(z_Xpert', N);
  
  z_Smear_RE[,2] += b[1]*bac_load;
  z_Mgit_RE[,2]  += b[2]*bac_load;
  z_Xpert_RE[,2] += b[3]*bac_load;
  
  // Imputation ---------------------------------------------------------------
  // - HIV
  Xd_imp[,1] = to_vector(impute_discrete(Xd[,1], obs_Xd[,1], rep_vector(Phi(HIV_a0), N)));
  
  // - Clinical symptoms
  for (n in 1:N_miss_cs)
    z_cs[n_miss_cs[n], d_miss_cs[n]] = z_miss_cs[n];
  for (n in 1:N_pos_cs)
    z_cs[n_pos_cs[n], d_pos_cs[n]] = z_pos_cs[n];
  for (n in 1:N_neg_cs)
    z_cs[n_neg_cs[n], d_neg_cs[n]] = z_neg_cs[n];
    
  Xd_imp[,2] = to_vector(impute_discrete(Xd[,2], obs_Xd[,2], to_array_2d(append_all(z_cs)), obs_cs));
  
  // - Motor palsy
  for (n in 1:N_miss_mp)
    z_mp[n_miss_mp[n], d_miss_mp[n]] = z_miss_mp[n];
  for (n in 1:N_pos_mp)
    z_mp[n_pos_mp[n], d_pos_mp[n]] = z_pos_mp[n];
  for (n in 1:N_neg_mp)
    z_mp[n_neg_mp[n], d_neg_mp[n]] = z_neg_mp[n];
    
  Xd_imp[,3] = to_vector(impute_discrete(Xd[,3], obs_Xd[,3], to_array_2d(append_all(z_mp)), obs_mp));
  
  // - Age
  Xc_imp[,1] = impute_cont(Xc[,1], obs_Xc[,1], age_imp);
  
  // - Illness day
  Xc_imp[,2] = impute_cont(Xc[,2], obs_Xc[,2], id_imp);
  
  // Other vars
  Xd_imp[,4] = Xd[,4];
  X = append_col(Xd_imp, Xc_imp);
  obs_X = append_col(to_matrix(obs_Xc), to_matrix(obs_Xd)); 
}

model {
  
  // Imputation ---------------------------------------------------------------
  // - HIV
  real p_HIV = Phi(HIV_a0);
  HIV_a0 ~ student_t(5, 0, 2.5);
  for (n in 1:N) {
    if (obs_Xd[n,1] == 1){
      Xd[n,1] ~ bernoulli(p_HIV);
    }
  }
  
  // - Clinical symptoms
  L_Omega_cs ~ lkj_corr_cholesky(4);
  cs_a ~ student_t(5, 0, 2.5);
  
  {
    vector[3] cs_a_x[N];
    for (n in 1:N)
      if (obs_Xd[n,1] == 1){
        cs_a_x[n] = cs_a * Xd[n,1];
        z_cs[n] ~ multi_normal_cholesky(cs_a_x[n], L_Omega_cs);
      } else {
        target += log_mix(p_HIV,
        multi_normal_cholesky_lpdf(z_cs[n] | cs_a, L_Omega_cs),
        multi_normal_cholesky_lpdf(z_cs[n] | [0, 0, 0], L_Omega_cs));
      }
  }
  
  // - Motor palsy
  L_Omega_mp ~ lkj_corr_cholesky(4);
  mp_a ~ student_t(5, 0, 2.5);
  
  {
    vector[3] mp_a_x[N];
    for (n in 1:N)
      if (obs_Xd[n,1] == 1){
        mp_a_x[n] = mp_a * Xd[n,1];
        z_mp[n] ~ multi_normal_cholesky(mp_a_x[n], L_Omega_mp);
      } else {
        target += log_mix(p_HIV,
        multi_normal_cholesky_lpdf(z_mp[n] | mp_a, L_Omega_mp),
        multi_normal_cholesky_lpdf(z_mp[n] | [0, 0, 0], L_Omega_mp));
      }
  }
  
  // - Age
  mu_age ~ student_t(5, 0, 10);
  sigma_age ~ normal(0, 2.5);
  Xc[,1] ~ normal(mu_age, sigma_age);
  
  // - Illness day
  mu_id ~ student_t(5, 0, 10);
  sigma_id ~ normal(0, 2.5);
  Xc[,2] ~ normal(mu_id, sigma_id);
  
  // --------------------------------------------------------------------------
  
  
  //Probs of each test become positive
  // Priors of covariates
  a0       ~ student_t(5, 0  ,10  );
  a[1]     ~ student_t(5, 0  , 2.5);
  a[2]     ~ student_t(5, 0  , 2.5);
  a[3]     ~ student_t(5, 2  , 1  );
  a[4]     ~ student_t(5, 1.7, 1  );
  a[5]     ~ student_t(5, 1  , 1  );
  a[6]     ~ student_t(5, 1  , 1  );
  a[7]     ~ student_t(5, 1  , 2.5);
  a[8]     ~ student_t(5,-1  , 2.5);
  a[9]     ~ student_t(5,-1  , 2.5);
  a[10:12] ~ student_t(5, 0  , 2.5);
  a[13]    ~ student_t(5, 2  , 1  );
  a[14]    ~ student_t(5, 4  , 1  );
  // a[15]     ~ student_t(5, 0  , 2.5);
  
  
  //Random effects covariates
  RE    ~    normal(   0, 1  );
  b_HIV ~ student_t(5, 0, 2.5);
  b     ~ student_t(5, 0, 1  );
  
  //1-Specificity of each test
  z_Xpert[1] ~ normal(inv_Phi(.005), .7  );
  z_Mgit[1]  ~ normal(-3.023       , .89 );
  z_Smear[1] ~ normal(-3.023       , .89 );
  
  //Sensitivity of each test
  z_Xpert[2] ~ normal(inv_Phi(.593), .117);
  z_Mgit[2]  ~ normal(inv_Phi(.665), .217);
  z_Smear[2] ~ normal(inv_Phi(.786), .405);
  
  for (n in 1:N){
    target += log_mix(theta[n],
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,2]),// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[2]),
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,1]));// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[1]));
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = log_mix(theta[n],
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,2]),// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[2]),
    bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[n,1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[n,1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[n,1]));// + bernoulli_logit_lpmf(Y_Img[n] | z_Img[1]));
  }
}
