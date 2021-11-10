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
    matrix[size(x), num_elements(x[1])] M;
    for (i in 1:size(x)) M[i,] = x[i]'; 
    return M;
  }
  
  //Calculate the sum of size of each vector in an array of vector
  int sum_size(vector[] x){
    int s = 0;
    int N = size(x);
    for (n in 1:N){
      s = s + num_elements(x[n]);
    }
    return s;
  }
  
  //Concat vectors into one
  vector concat(vector[] x){
    int N = size(x);
    int S = sum_size(x);
    vector[S] y;
    int i = 1;
    
    for (n in 1:N){
      y[i:(i+num_elements(x[n])-1)] = x[n];
      i += num_elements(x[n]);
    }
    
    return y;
  }
  
  //Rep vectors into one vector
  vector rep_vectors(vector x, int n){
    vector[num_elements(x) * n] y;
    
    for (i in 1:n){
      y[(num_elements(x)*(i-1) + 1) : num_elements(x)*i] = x;
    }
    
    return y;
  }
  
  //Function to sum over 2|3 "or" logic.
  real sum_probs_or(real[] probs){
    int s = size(probs);
    real sum_prob;
    
    if (s>3||s<2) reject("Only size 2 and 3 support");
    sum_prob = (s==2) ? (sum(probs) - prod(probs)) : (sum(probs) - prod(probs[1:2]) - prod(probs[2:3]) - prod(probs[{1,3}]) + prod(probs));
    
    return sum_prob;
  }
  
  // Impute discrete variable. Missing value will be "imputed" by their expected probs
  real[] impute_discrete_1d(int[] raw, int[] obs, real[] z){
    int N = size(raw);
    real imp[N];
    
    if (size(obs)!=N||size(z)!=N) reject("Size mismatched!");
    
    for (n in 1:N)
    imp[n] = obs[n] ? raw[n] : Phi(z[n]);
    
    return imp;
  }
  
  // Impute discrete variable for combined val. Missing value will be "imputed" by their expected probs
  vector impute_discrete_2d(int[] raw_cmb, int[] obs_cmb, real[,] z_el, int[,] obs_el){
    int N = size(raw_cmb);
    vector[N] cmb;
    
    if (size(obs_cmb)!=N||size(z_el)!=N||size(obs_el)!=N) reject("Size mismatched!");
    
    for (n in 1:N){
      if (obs_cmb[n]) cmb[n] = raw_cmb[n];
      else {
        int n_obs_el = sum(obs_el[n]);
        if (n_obs_el == 0) cmb[n] = sum_probs_or(Phi(z_el[n]));
        else if (n_obs_el == 1) {
          if (obs_el[n,1] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,2:3]));
          if (obs_el[n,2] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,{1,2}]));
          if (obs_el[n,3] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,1:2]));
        }
        else if (n_obs_el == 2){
          if (obs_el[n,1] == 0) cmb[n] = Phi(z_el[n,1]);
          if (obs_el[n,2] == 0) cmb[n] = Phi(z_el[n,2]);
          if (obs_el[n,3] == 0) cmb[n] = Phi(z_el[n,3]);
        }
      }
    }
    return cmb;
  }
  
  // Impute continuous variable
  vector impute_cont_1d(real[] raw, int[] obs, real[] rep){
    int N = size(raw);
    vector[N] imp;
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
  
  matrix impute_cont_2d(real[,] raw, int[,] obs, real[] rep){
    // int N =;
    matrix[dims(raw)[1], dims(raw)[2]] imp;
    int k = 1;
    
    if (size(obs)!=size(raw) || size(obs[1])!=size(raw[1]))
    reject("Size mismatched!");
    
    for (i in 1:dims(raw)[1]){
      for (j in 1:dims(raw)[2]){
        if (obs[i,j]) imp[i,j] = raw[i,j];
        else {
          imp[i,j] = rep[k];
          k += 1;
        }
      }
    }
    return imp;
  }
  
  // This function for rng logit with awarenss of observervation
  vector discrete_rng(vector imputed_1d, int[] obs){
    int N = num_elements(imputed_1d);
    vector[N] val;
    if (size(obs) != N) reject("Size mismatched!");
    
    for (n in 1:N) val[n] = obs[n] ? imputed_1d[n] : bernoulli_rng(imputed_1d[n]);
    
    return val;
  }
  
  matrix discrete_2d_rng(matrix imputed_2d, int[,] obs){
    int M = dims(imputed_2d)[2];
    int N = dims(imputed_2d)[1];
    matrix[N, M] val;
    if (dims(obs)[1] != N || dims(obs)[2] != M) reject("Size mismatched!");
    for (m in 1:M) val[, m] = discrete_rng(imputed_2d[,m], obs[, m]);
    return val;
  }
  
  // This function do the power operator but return int rather than real
  int int_power(int x, int n){
    int X = 1;
    if (n < 0) reject("n must be natural");
    if (n == 0){
      if (X == 0) reject("0^0 is not allowed!");
      return 1;
    } else {
      for (n_ in 1:n){
        X *= x;
      }
      return X;
    }
  }
  
  vector[] get_patterns(row_vector Xd_imp, int[] obs_Xd_imp, vector a_Xd_imp){
    // X_imp are discrete X with missing, X_cmpl are discrete X with no missing and continuous X
    int N_Xd_imp = size(obs_Xd_imp);
    int N_miss = N_Xd_imp - sum(obs_Xd_imp);
    int T = int_power(2, N_miss);
    matrix[T, N_Xd_imp] obs_pattern; //Construction vector
    matrix[T, N_Xd_imp] probs_pattern; //Probability vector
    matrix[T, N_Xd_imp] a_pattern; //Coefficent vector
    
    vector[T] pat = rep_vector(0, T);
    vector[T] log_probs = rep_vector(0, T);
    // vector[T] val[2];
    int j = 0;
    if (size(obs_Xd_imp) != N_Xd_imp) reject("Size mismatched!");
    
    for (i in 1:N_Xd_imp){
      vector[int_power(2, j+1)] V;
      if (obs_Xd_imp[i] == 0){
        V = append_row(rep_vector(0, int_power(2, j)), rep_vector(1, int_power(2, j)));
        obs_pattern[,i] = rep_vectors(V, T / int_power(2, j+1));
        j += 1;
      } else {
        obs_pattern[,i] = rep_vector(Xd_imp[i], T);
      }
    }
    for (i in 1:N_Xd_imp){
      for (t in 1:T){
        if (obs_pattern[t,i] == 1){
          probs_pattern[t,i] = Xd_imp[i];
          a_pattern[t,i] = a_Xd_imp[i];
        } else {
          probs_pattern[t,i] = 1 - Xd_imp[i];
          a_pattern[t,i] = 0;
        }
      }
      pat += a_pattern[,i];
      log_probs += log(probs_pattern[,i]);
    }
    
    return {log_probs, pat};
  }
}

data {
  int<lower=1> N; //Number of patient
  int<lower=1> nXc; //Number of continuous X
  int<lower=1> nXd; //Number of discrete X
  int<lower=1> nTd; //Number of disc. aux covariates for imputation model
  
  int<lower=0, upper=1> Y_Smear0[N];
  int<lower=0, upper=1> Y_Mgit0[N];
  int<lower=0, upper=1> Y_Xpert0[N];
  
  real Xc0[N, nXc]; //Continuous covariates
  int  Xd0[N, nXd]; //Discrete covariates
  int  Td0[N, nTd]; //Auxillary covariates - discrete
  
  //Matrices of observation
  int<lower=0, upper=1> obs_Xc0[N, nXc]; 
  int<lower=0, upper=1> obs_Xd0[N, nXd];
  int<lower=0, upper=1> obs_Td0[N, nTd]; 
  
  //Bootstraps
  int<lower=0, upper=1> resample;
}

transformed data {
  simplex[N] uniform = rep_vector(1.0 / N, N);
  int<lower = 1, upper = N> boot_idxs[N];
 
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
  
  
  // Var declarations ..........................................................
  int nX = nXc + nXd;

  
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
  int<lower=1,upper=N> n_neg_mp[sum2d(obs_mp) - size(n_pos_mp)];
  int<lower=1,upper=3> d_neg_mp[size(n_neg_mp)];

  // For csf mvNormal regression
  // This is a bug in RStan when using custom functions in parameters.
  int obs_csf_glu = sum2d(obs_Xc[,3:4]);
  int obs_csf_other = sum2d(obs_Xc[,5:7]);


  // Var assigments.............................................................
   for (n in 1:N)
    boot_idxs[n] = resample ? categorical_rng(uniform) : n;
    
  Y_Smear = Y_Smear0[boot_idxs];
  Y_Mgit = Y_Mgit0[boot_idxs];
  Y_Xpert = Y_Xpert0[boot_idxs];
  Xc = Xc0[boot_idxs,:];
  Xd = Xd0[boot_idxs,:];
  Td = Td0[boot_idxs,:];
  obs_Xc = obs_Xc0[boot_idxs,:];
  obs_Xd = obs_Xd0[boot_idxs,:];
  obs_Td = obs_Td0[boot_idxs,:];
  
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
  vector[3] cs_a0;
  vector[3] cs_a;
  cholesky_factor_corr[3] L_Omega_cs;
  vector<lower=0>[N_pos_cs] z_pos_cs;
  vector<upper=0>[N_neg_cs] z_neg_cs;
  vector[N_miss_cs] z_miss_cs;
  
  // Impute motor palsy Td[4:6] for components variables
  // Xd[3] is the combined one aka clin_motor_palsy
  vector[3] mp_a0;
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
  
  
  // Impute {log2(bldglu) and [log2(csfglu)} and log_lympho, log_protein, log_lactate]
  // Impute (sqrt(bldglu), sqrt(csfglu), log_lympho, log_protein, log_lactate) as mvNormal
  vector[5] csf_a0;
  vector[2] glu_a_diabete;
  cholesky_factor_corr[5] L_Omega_csf;
  real<lower=0> bld_csf_glu_imp[(N*2) - obs_csf_glu];
  real csf_other_imp[(N*3) - obs_csf_other];
  
  //----------------------------------------------------------------------------
  
  // For logit regression ------------------------------------------------------
  real a0; //intercept
  // vector[nX + 1] a; //slope
  vector<lower=0>[3] a_pos; //1,2, #5, 6, #8
  // vector<upper=0>[2] a_neg;
  vector[nX + 1 - 3] a_;
  real b_HIV; //adjustment of RE with HIV X[,1]
  real<lower=0> b_RE;
  
  //Probability of each vars
  ordered[2] z_Smear; 
  ordered[2] z_Mgit;
  ordered[2] z_Xpert;
  
  vector[N] RE; //base random effect;
}

transformed parameters{
  vector[nX + 1] a = concat({a_pos[1:2], a_[1:3], a_pos[3:3], a_[4:(nX+1-3)]});
  // Create an fully imputed X matrix
  matrix[N, nXd - 3] Xd_compl; 
  matrix[N, 3] Xd_imp; //fully imputed discrete X
  matrix[N, nXc + 1] Xc_imp; //fully imputed cont X
  
  // Imputation models ---------------------------------------------------------
  // - Clinical symptoms
  vector[3] z_cs[N];
  
  // - Motor palsy
  vector[3] z_mp[N];
  
  // - HIV
  Xd_imp[,1] = to_vector(impute_discrete_1d(Xd[,1], obs_Xd[,1], rep_array(Phi(HIV_a0), N)));
  
  // - Clinical symptoms
  for (n in 1:N_miss_cs)
    z_cs[n_miss_cs[n], d_miss_cs[n]] = z_miss_cs[n];
  for (n in 1:N_pos_cs)
    z_cs[n_pos_cs[n], d_pos_cs[n]] = z_pos_cs[n];
  for (n in 1:N_neg_cs)
    z_cs[n_neg_cs[n], d_neg_cs[n]] = z_neg_cs[n];
    
  Xd_imp[,2] = to_vector(impute_discrete_2d(Xd[,2], obs_Xd[,2], to_array_2d(append_all(z_cs)), obs_cs));
  
  // - Motor palsy
  for (n in 1:N_miss_mp)
    z_mp[n_miss_mp[n], d_miss_mp[n]] = z_miss_mp[n];
  for (n in 1:N_pos_mp)
    z_mp[n_pos_mp[n], d_pos_mp[n]] = z_pos_mp[n];
  for (n in 1:N_neg_mp)
    z_mp[n_neg_mp[n], d_neg_mp[n]] = z_neg_mp[n];
    
  Xd_imp[,3] = impute_discrete_2d(Xd[,3], obs_Xd[,3], to_array_2d(append_all(z_mp)), obs_mp);
  
  // - Age
  Xc_imp[,1] = impute_cont_1d(Xc[,1], obs_Xc[,1], age_imp);
  
  // - Illness day
  Xc_imp[,2] = impute_cont_1d(Xc[,2], obs_Xc[,2], id_imp);
  
  // - Blood glucose, CSF glucose, lymphocyte, protein, lactate
  // blood glucose, csf_glucose, csf lymphocytes count, protein, and lactate
  // logbldglu: Xc[4], logcsfglu: Xc[5], loglym: Xc[6], logprotein: Xc[7], loglac: Xc[8]
  // real CSF[N,5] = Xc[,4:8]; //raw BLDGLU, CSFGLU, LYMPH, PROTEIN, LACTATE
  // int obs_csf[N,5] = obs_Xc[,4:8];
  Xc_imp[,3:7] = impute_cont_2d(Xc[,3:7], obs_Xc[,3:7], append_array(bld_csf_glu_imp, csf_other_imp));
  Xc_imp[,8] = Xc_imp[,4] ./ Xc_imp[,3];
  
  
  // Other vars
  Xd_compl = to_matrix(Xd[,4:nXd]); 
}

model {
  matrix[N, nX - 3 + 1] X_compl = append_col(Xd_compl, Xc_imp);
  
  // Imputation ---------------------------------------------------------------
  // - HIV
  real p_HIV = Phi(HIV_a0);
  HIV_a0 ~ student_t(5, 0, 1);
  
  // - Clinical symptoms
  L_Omega_cs ~ lkj_corr_cholesky(4);
  cs_a0 ~ student_t(5, 0, 1);
  cs_a ~ student_t(5, 0, 1);
  
  
  // - Motor palsy
  L_Omega_mp ~ lkj_corr_cholesky(4);
  mp_a0 ~ student_t(5, 0, 1);
  mp_a ~ student_t(5, 0, 1);

  // - Age
  mu_age ~ student_t(5, 0, 1);
  sigma_age ~ normal(0, 1);
  Xc_imp[,1] ~ normal(mu_age, sigma_age);
  
  // - Illness day
  mu_id ~ student_t(5, 0, 1);
  sigma_id ~ normal(0, 1);
  Xc_imp[,2] ~ normal(mu_id, sigma_id);
  
  // - Blood glucose and glucose ratio
  to_vector(glu_a_diabete) ~ student_t(5, 0, 1);
  L_Omega_csf ~ lkj_corr_cholesky(4);
  csf_a0 ~  student_t(5, 0, 1);
  
  
  {
    vector[5] csf_a_x[N];
    vector[5] csf_imp[N];
    for (n in 1:N){
      csf_imp[n] = Xc_imp[n,3:7]';
      csf_a_x[n, 1:2] = csf_a0[1:2] + glu_a_diabete * Td[n, 7];
      csf_a_x[n, 3:5] = csf_a0[3:5];
    }
    csf_imp ~ multi_normal_cholesky(csf_a_x, L_Omega_csf);
  }

// --------------------------------------------------------------------------
  
  //Probs of each test become positive
// Priors of covariates
a0       ~ student_t(5, 0  , 2.5);
a        ~ student_t(5, 0  , 2.5);

//Random effects covariates
RE    ~    normal(   0, 1);
b_RE  ~ student_t(5, 0, 1);
b_HIV ~ student_t(5, 0, 1);
// b     ~ student_t(5, 0, 1);

//1-Specificity of each test
z_Xpert[1] ~ normal(logit(.005), 1.59);
z_Mgit[1]  ~ normal(logit(.001),  .7 ); //.82
z_Smear[1] ~ normal(logit(.001),  .7 );

//Sensitivity of each test
z_Xpert[2] ~ normal(0, .6);
z_Mgit[2]  ~ normal(0, .6);
z_Smear[2] ~ normal(0, .6);

for (n in 1:N){
  int N_Xd_miss = 3 - sum(obs_Xd[n, 1:3]);

  if (obs_Xd[n,1] == 1){
    Xd[n,1] ~ bernoulli(p_HIV);
    z_mp[n] ~ multi_normal_cholesky(mp_a0 + mp_a * Xd[n,1], L_Omega_mp);
    z_cs[n] ~ multi_normal_cholesky(cs_a0 + cs_a * Xd[n,1], L_Omega_cs);
  } else {
  if (is_nan(
    multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a, L_Omega_cs) +
    multi_normal_cholesky_lpdf(z_cs[n] | cs_a0, L_Omega_cs) +
    multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a, L_Omega_mp) +
    multi_normal_cholesky_lpdf(z_mp[n] | mp_a0, L_Omega_mp)))
    target += not_a_number();
    else {
  target += log_mix(p_HIV,
                    multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a, L_Omega_cs),
                    multi_normal_cholesky_lpdf(z_cs[n] | cs_a0, L_Omega_cs));
  target += log_mix(p_HIV,
                    multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a, L_Omega_mp),
                    multi_normal_cholesky_lpdf(z_mp[n] | mp_a0, L_Omega_mp));
  }
  }
  
  if (N_Xd_miss){
    int N_pattern = int_power(2, N_Xd_miss);
    vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,], obs_Xd[n, 1:3], a[1:3]);
    vector[N_pattern] log_liks;
    pat_thetas[2] += a0 + dot_product(a[4:15], X_compl[n]);
    
    //check if HIV is missing
    if (obs_Xd[n,1]){
      for (i in 1:N_pattern){
        real logprob_theta = pat_thetas[1][i];
        real theta = inv_logit(pat_thetas[2][i]);
        
        real bac_load = b_HIV * Xd_imp[n,1];
        real z_Smear_RE = z_Smear[2] + bac_load + b_RE * RE[n];
        real z_Mgit_RE  = z_Mgit [2] + bac_load + b_RE * RE[n];
        real z_Xpert_RE = z_Xpert[2] + bac_load + b_RE * RE[n];
        
        log_liks[i] = logprob_theta + log_mix(theta, 
                                              bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
                                              bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
      }
    } else {
      for (i in 1:N_pattern){
        real logprob_theta = pat_thetas[1][i];
        real theta = inv_logit(pat_thetas[2][i]);
        // real theta = (i % 2 == 1) ? inv_logit(pat_thetas[2][i]) : inv_logit(pat_thetas[2][i]);
        
        vector[2] pat_bac_load[2] = get_patterns([Xd_imp[n,1]], {0}, [b_HIV]');
        vector[2] logprob_Y = pat_bac_load[1];
        vector[2] bac_load = pat_bac_load[2];
        
        vector[2] z_Smear_RE = z_Smear[2] + bac_load + b_RE * RE[n];
        vector[2] z_Mgit_RE  = z_Mgit [2] + bac_load + b_RE * RE[n];
        vector[2] z_Xpert_RE = z_Xpert[2] + bac_load + b_RE * RE[n];
        
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
      row_vector[nX + 1] X = append_col(Xd_imp[n,], X_compl[n,]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load = b_HIV*Xd_imp[n, 1];
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
  
  matrix<lower=0, upper=1>[N, 3] Xd_rng = discrete_2d_rng(Xd_imp, obs_Xd[,1:3]);
  matrix[N, nX + 1] X = append_col(append_col(Xd_rng, Xd_compl), Xc_imp);
  
  theta = inv_logit(a0 + X*a[1:15]); 
  C = bernoulli_rng(theta);
  
  for (n in 1:N){
    if (C[n]){
      log_lik[n] = bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]);
    } else {
      real bac_load = b_HIV * X[n,1] + b_RE * RE[n];
      log_lik[n] = bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[2] + bac_load) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[2] + bac_load) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[2] + bac_load);
    }
  }
}

