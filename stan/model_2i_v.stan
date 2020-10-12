functions {
  // int[] equal_vec(vector x, real y){
  //   int eq[num_elements(x)];
  //   for (i in 1:num_elements(x)){
  //     eq[i] = x[i] == y ? 1 : 0;
  //   }
  //   
  //   return eq;
  // }
  
  int[] greater_vect(vector x, real y){
    int gr[num_elements(x)];
    for (i in 1:num_elements(x)){
      gr[i] = x[i] >= y ? 1 : 0;
    }
    
    return gr;
  }
  
  // int[] lesser_vect(vector x, real y){
  //   int ls[num_elements(x)];
  //   for (i in 1:num_elements(x)){
  //     ls[i] = x[i] <= y ? 1 : 0;
  //   }
  //   
  //   return ls;
  // }
  
  int[] which(vector x, real threshold){
    int x_int[num_elements(x)] = greater_vect(x, threshold);
    int w[sum(x_int)];
    int k=1;
    for (i in 1:num_elements(x)){
      if (x_int[i] == 1){
        w[k] = i;
        k += 1;
      } 
    }
    return w;
  }
  
  int[] anti_which(vector x, real threshold){
    int x_int[num_elements(x)] = greater_vect(x, threshold);
    int w[sum(x_int)];
    int k=1;
    for (i in 1:num_elements(x)){
      if (x_int[i] == 0){
        w[k] = i;
        k += 1;
      } 
    }
    return w;
  }
  
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
    vector[T] prod_probs = rep_vector(1, T);
    // vector[T] val[2];
    int j = 0;
    if (size(obs_Xd_imp) != N_Xd_imp) reject("Size mismatched!");
    
    for (i in 1:N_Xd_imp){
      vector[int_power(2, j+1)] V;
      if (obs_Xd_imp[i] == 0){
        V = append_row(rep_vector(0, int_power(2, j)), rep_vector(1, int_power(2, j)));
        obs_pattern[,i] = rep_vectors(V, T/int_power(2, j+1));
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
      prod_probs .*= probs_pattern[,i];
    }
   
    return {prod_probs, pat};
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
  real Tc[N, nTc]; //Auxillary covariates - cont
  
  //Matrices of observation
  int<lower=0, upper=1> obs_Xc[N, nXc]; 
  int<lower=0, upper=1> obs_Xd[N, nXd];
  int<lower=0, upper=1> obs_Td[N, nTd]; 
  int<lower=0, upper=1> obs_Tc[N, nTc]; 
}

transformed data {
  
  // Var declarations ..........................................................
  int nX = nXc + nXd;
  int i_obs_HIV[sum(obs_Xd[,1])] = which(to_vector(obs_Xd[,1]), 1);
  int i_miss_HIV[N-sum(obs_Xd[,1])] = anti_which(to_vector(obs_Xd[,1]), 1);
  int i_obs_all[N*3-sum(to_array_1d(obs_Xd[,1:3]))] = which([sum(obs_Xd[, 1]), sum(obs_Xd[,2]), sum(obs_Xd[,3])]', 3);
  
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

  // Var assigments.............................................................
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
  
  // Impute gcs Xc[3]
  // Use auxillary gcs components: Tc[1:3] GCSE, GCSM, GCSV
  // GCSV is missing.
  real<lower=0> gcsv_a0;
  vector[2] gcsv_a;
  real<lower=0> gcsv_sigma;
  real<lower=0, upper=4> gcsv_imp[N - sum(obs_Tc[,3])];
  
  // Impute {log2(bldglu) and [log2(csfglu)} and log_lympho, log_protein, log_lactate]
  // csfglu stays in both groups as to form a link
  vector[2] glu_a0;
  matrix[2,2] glu_a; //two covs are "diabete" and "age"
  vector[4] csf_a0; // no covs for the other 4
  cholesky_factor_corr[2] L_Omega_glu;
  cholesky_factor_corr[4] L_Omega_csf;
  real csf_glu_imp[(N*5) - sum2d(obs_Xc[,4:8])];
  
  
  //----------------------------------------------------------------------------
  
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
  matrix[N, nXd - 3] Xd_compl; 
  matrix[N, 3] Xd_imp; //fully imputed discrete X
  matrix[N, nXc] Xc_imp; //fully imputed cont X
  
  // Imputation models ---------------------------------------------------------
  // - Clinical symptoms
  vector[3] z_cs[N];
  
  // - Motor palsy
  vector[3] z_mp[N];
  
  // - GCSV
  vector<lower=0, upper=5>[N] GCSV_imp;
  
  // - HIV
  Xd_imp[,1] = to_vector(impute_discrete_1d(Xd[,1], obs_Xd[,1], rep_array(Phi(HIV_a0), N)));
  
  // - Clinical symptoms
  for (n in 1:N_miss_cs)
    z_cs[n_miss_cs[n], d_miss_cs[n]] = z_miss_cs[n];
  for (n in 1:N_pos_cs)
    z_cs[n_pos_cs[n], d_pos_cs[n]] = z_pos_cs[n];
  for (n in 1:N_neg_cs)
    z_cs[n_neg_cs[n], d_neg_cs[n]] = z_neg_cs[n];
    
  Xd_imp[,2] = to_vector(impute_discrete_2d(Xd[,2], obs_Xd[,2], to_array_2d(Phi(append_all(z_cs))), obs_cs));
  
  // - Motor palsy
  for (n in 1:N_miss_mp)
    z_mp[n_miss_mp[n], d_miss_mp[n]] = z_miss_mp[n];
  for (n in 1:N_pos_mp)
    z_mp[n_pos_mp[n], d_pos_mp[n]] = z_pos_mp[n];
  for (n in 1:N_neg_mp)
    z_mp[n_neg_mp[n], d_neg_mp[n]] = z_neg_mp[n];
    
  Xd_imp[,3] = impute_discrete_2d(Xd[,3], obs_Xd[,3], to_array_2d(Phi(append_all(z_mp))), obs_mp);
  
  // - Age
  Xc_imp[,1] = impute_cont_1d(Xc[,1], obs_Xc[,1], age_imp);
  
  // - Illness day
  Xc_imp[,2] = impute_cont_1d(Xc[,2], obs_Xc[,2], id_imp);
  
  // - GCSV
  GCSV_imp = impute_cont_1d(Tc[,3], obs_Tc[,3], gcsv_imp);
  Xc_imp[,3] = impute_cont_1d(Xc[,3], obs_Xc[,3], 
                         to_array_1d(to_vector(Tc[,1]) + to_vector(Tc[,2]) + to_vector(GCSV_imp)));
  
  // - Blood glucose, CSF glucose, lymphocyte, protein, lactate
  // Glucose ratio, blood glucose, csf lymphocytes count, protein, and lactate
  // gluratio: Xc[4], bldglu: Xc[5], loglym: Xc[6], logprotein: Xc[7], loglac: Xc[8]
  // real CSF[N,5] = Xc[,4:8]; //raw BLDGLU, GLURATIO, LYMPH, PROTEIN, LACTATE
  // int obs_csf[N,5] = obs_Xc[,4:8];
  Xc_imp[,4:8] = impute_cont_2d(Xc[,4:8], obs_Xc[,4:8], csf_glu_imp);
  
  // Other vars
  Xd_compl = to_matrix(Xd[,4:6]); 
}

model {
  matrix[N, nX - 3] X_compl = append_col(Xd_compl, Xc_imp);
  
  // Imputation ---------------------------------------------------------------
  // - HIV
  real p_HIV = Phi(HIV_a0);
  HIV_a0 ~ student_t(5, 0, 1);
  Xd[i_obs_HIV, 1] ~ bernoulli(p_HIV);
  
  // - Clinical symptoms
  L_Omega_cs ~ lkj_corr_cholesky(4);
  cs_a0 ~ student_t(5, 0, 2.5);
  cs_a ~ student_t(5, 0, 2.5);
  
  // - Motor palsy
  L_Omega_mp ~ lkj_corr_cholesky(4);
  mp_a0 ~ student_t(5, 0, 1);
  mp_a ~ student_t(5, 0, 2.5);
  
  {
    vector[3] mp_a_x[size(i_obs_HIV)];
    vector[3] cs_a_x[size(i_obs_HIV)];
    int k = 1;
    for (n in i_obs_HIV){
      mp_a_x[k] = mp_a0 + mp_a * Xd[n,1];
      cs_a_x[k] = cs_a0 + cs_a * Xd[n,1];
      k += 1;
    }
    z_mp[i_obs_HIV] ~ multi_normal_cholesky(mp_a_x, L_Omega_mp);
    z_cs[i_obs_HIV] ~ multi_normal_cholesky(cs_a_x, L_Omega_cs);
  }
  
  
  for (n in i_miss_HIV){
    target += log_mix(p_HIV,
                      multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a, L_Omega_cs),
                      multi_normal_cholesky_lpdf(z_cs[n] | cs_a0, L_Omega_cs));
    target += log_mix(p_HIV,
                      multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a, L_Omega_mp),
                      multi_normal_cholesky_lpdf(z_mp[n] | mp_a0, L_Omega_mp));
  }
   
  // - Age
  mu_age ~ student_t(5, 0, 1);
  sigma_age ~ normal(0, 2.5);
  Xc_imp[,1] ~ normal(mu_age, sigma_age);
  
  // - Illness day
  mu_id ~ student_t(5, 0, 1);
  sigma_id ~ normal(0, 1);
  Xc_imp[,2] ~ normal(mu_id, sigma_id);
  
    // GCSV
  gcsv_a0 ~ student_t(5, 0, 1);
  gcsv_a  ~ student_t(5, 0, 2.5);
  gcsv_sigma ~ normal(0, 1);
  
  {
    vector[N] gcsv_x;
    gcsv_x = gcsv_a0 + to_matrix(Tc[,1:2])*gcsv_a;
    GCSV_imp ~ normal(gcsv_x, gcsv_sigma);
  }
  
  // - Blood glucose and glucose ratio
  // Td[7] is diabetes status
  // Xc_imp[1] is imputed age
  L_Omega_glu ~ lkj_corr_cholesky(4);
  glu_a0 ~ student_t(5, 0, 1);
  to_vector(glu_a) ~ student_t(5, 0, 2.5);
  // Other csf lab values
  L_Omega_csf ~ lkj_corr_cholesky(4);
  csf_a0 ~  student_t(5, 0, 1);
  
  {
    vector[2] glu_a_x[N];
    vector[2] glu_imp[N];
    vector[4] csf_imp[N];
    for (n in 1:N){
      glu_imp[n] = [Xc_imp[n,3], Xc_imp[n,4]]';
      glu_a_x[n] = glu_a0 + glu_a[,1]*Xc_imp[n,1] + glu_a[,2]*Td[n,7];
      csf_imp[n] = [Xc_imp[n,4], Xc_imp[n,5], Xc_imp[n,6], Xc_imp[n,7]]';
    }
      
      glu_imp ~ multi_normal_cholesky(glu_a_x, L_Omega_glu);
      csf_imp ~ multi_normal_cholesky(csf_a0, L_Omega_csf);
  }
  
  // --------------------------------------------------------------------------
  
  //Probs of each test become positive
  // Priors of covariates
  a0       ~ student_t(5, 0  ,10  );
  
  // Discrete var
  a[1]     ~ student_t(5, 0  , 2.5);
  a[2]     ~ student_t(5, 2  , 1  );
  a[3]     ~ student_t(5, 1  , 1  );
  a[4]     ~ student_t(5, 1  , 1  );
  a[5]     ~ student_t(5, 2  , 1  );
  a[6]     ~ student_t(5, 4  , 1  );
  
  a[7]     ~ student_t(5, 0  , 2.5);
  a[8]     ~ student_t(5, 1.7, 1  );
  a[9]     ~ student_t(5, 1  , 2.5);
  a[10]    ~ student_t(5,-1  , 2.5);
  a[11]    ~ student_t(5,-1  , 2.5);
  a[12:14] ~ student_t(5, 0  , 2.5);
  
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
    int N_Xd_miss = 3 - sum(obs_Xd[n, 1:3]);
    
    if (N_Xd_miss){
      int N_pattern = int_power(2, N_Xd_miss);
      vector[N_pattern] pat_thetas[2] = get_patterns(Xd_imp[n,], obs_Xd[n, 1:3], a[1:3]);
      vector[N_pattern] log_liks;
      pat_thetas[2] += a0 + dot_product(a[4:nX], X_compl[n]);
      
      //check if HIV is missing
      if (obs_Xd[n,1]){
        for (i in 1:N_pattern){
          real prob_theta = pat_thetas[1][i];
          real theta = inv_logit(pat_thetas[2][i]);
          
          real bac_load = RE[n] + b_HIV * Xd_imp[n,1];
          real z_Smear_RE = z_Smear[2] + b[1]*bac_load;
          real z_Mgit_RE  = z_Mgit [2] + b[2]*bac_load;
          real z_Xpert_RE = z_Xpert[2] + b[3]*bac_load;
          
          log_liks[i] = prob_theta*log_mix(theta, 
                                           bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
                                           bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      } else {
        for (i in 1:N_pattern){
          real prob_theta = pat_thetas[1][i];
          real theta = inv_logit(pat_thetas[2][i]);
          
          vector[2] pat_bac_load[2] = get_patterns([Xd_imp[n,1]], {0}, [b_HIV]');
          vector[2] prob_Y = pat_bac_load[1];
          vector[2] bac_load = pat_bac_load[2] + RE[n];
          
          vector[2] z_Smear_RE = z_Smear[2] + b[1]*bac_load;
          vector[2] z_Mgit_RE  = z_Mgit [2] + b[2]*bac_load;
          vector[2] z_Xpert_RE = z_Xpert[2] + b[3]*bac_load;
          
          log_liks[i] = prob_theta*log_mix(theta, 
                                           log_sum_exp(
                                             prob_Y[1] * (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[1])),
                                             prob_Y[2] * (bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE[2]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE[2]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE[2]))
                                           ),
                                           bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
        }
      } 
      
      target += log_sum_exp(log_liks);
    } else {
      // The normal way
      row_vector[nX] X = append_col(Xd_imp[n,], X_compl[n,]);
      real theta = inv_logit(a0 + dot_product(a, X));
      
      real bac_load = RE[n] + b_HIV*Xd_imp[n, 1];
      real z_Smear_RE = z_Smear[2] + b[1]*bac_load;
      real z_Mgit_RE  = z_Mgit [2] + b[2]*bac_load;
      real z_Xpert_RE = z_Xpert[2] + b[3]*bac_load;
      
      target += log_mix(theta, 
                        bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert_RE) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit_RE) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear_RE),
                        bernoulli_logit_lpmf(Y_Xpert[n] | z_Xpert[1]) + bernoulli_logit_lpmf(Y_Mgit[n] | z_Mgit[1]) + bernoulli_logit_lpmf(Y_Smear[n] | z_Smear[1]));
    }
  }
}

