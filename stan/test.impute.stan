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
  
  matrix append_all(vector[] x){
    matrix[size(x), dims(x[1])[1]] M;
    for (i in 1:size(x)) M[i,] = x[i]'; 
    return M;
  }
  
  real sum_probs_or(real[] probs){
    int s = size(probs);
    real sum_prob;
    
    if (s>3||s<2) reject("Only size 2 and 3 support");
    sum_prob = (s==2) ? sum(probs) - prod(probs) : sum(probs) - prod(probs[1:2]) - prod(probs[2:3]) - prod(probs[{1,3}]) + prod(probs);
    
    return sum_prob;
  }
  
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
}

data {
  int<lower = 1> N; //Number of patient
  
  int<lower = 0, upper = 1> HIV[N];
  int<lower = 0, upper = 1> cs[N,3];
  int<lower = 0, upper = 1> cbm_cs[N];
  
  int<lower=0, upper=1> obs_cs[N, 3]; //Matrix of observation
  int<lower=0, upper=1> obs_HIV[N]; //Matrix of observation
  int<lower=0, upper=1> obs_cbm_cs[N];
}

transformed data {
  int<lower=0> N_miss_cs;
  int<lower=1,upper=N> n_miss_cs[(N * 3) - sum2d(obs_cs)];
  int<lower=1,upper=3> d_miss_cs[size(n_miss_cs)];
  int<lower=0> N_pos_cs;
  int<lower=1,upper=N> n_pos_cs[sum2d_with_missing(cs, obs_cs)];
  int<lower=1,upper=3> d_pos_cs[size(n_pos_cs)];
  int<lower=0> N_neg_cs;
  int<lower=1,upper=N> n_neg_cs[sum2d(obs_cs) - size(n_pos_cs)];
  int<lower=1,upper=3> d_neg_cs[size(n_neg_cs)];

  N_pos_cs = size(n_pos_cs);
  N_neg_cs = size(n_neg_cs);
  N_miss_cs = size(n_miss_cs);
  
  {
    int i;
    int j;
    int k;
    i=1;j=1;k=1;
    for (n in 1:N) {
      for (d in 1:3) {
        if (obs_cs[n,d] == 0){
          n_miss_cs[k] = n;
          d_miss_cs[k] = d;
          k += 1;
        } else if (cs[n,d] == 1) {
          n_pos_cs[i] = n;
          d_pos_cs[i] = d;
          i += 1;
        } else {
          n_neg_cs[j] = n;
          d_neg_cs[j] = d;
          print(j);
          j += 1;
        }
      }
    }
  }
}

parameters{
  real HIV_a0;
  
  vector[3] cs_a;
  cholesky_factor_corr[3] L_Omega_cs;
  vector<lower=0>[N_pos_cs] z_pos_cs;
  vector<upper=0>[N_neg_cs] z_neg_cs;
  vector[N_miss_cs] z_miss_cs;
}

transformed parameters{
  vector[3] z_cs[N];
  real clin_sym[N];
   for (n in 1:N_miss_cs)
    z_cs[n_miss_cs[n], d_miss_cs[n]] = z_miss_cs[n];
  for (n in 1:N_pos_cs)
    z_cs[n_pos_cs[n], d_pos_cs[n]] = z_pos_cs[n];
  for (n in 1:N_neg_cs)
    z_cs[n_neg_cs[n], d_neg_cs[n]] = z_neg_cs[n];
  clin_sym = impute_discrete(cbm_cs, obs_cbm_cs, to_array_2d(append_all(z_cs)), obs_cs);
}

model {
  
  // -- HIV
  real p_HIV = Phi(HIV_a0);
  HIV_a0 ~ student_t(5, 0, 2.5);
  for (n in 1:N) {
    if (obs_HIV[n] == 1){
      HIV[n] ~ bernoulli(p_HIV);
    }
  }
  
    // -- Clinical symptoms
  L_Omega_cs ~ lkj_corr_cholesky(4);
  cs_a ~ student_t(5, 0, 2.5);
  
  {
    vector[3] cs_a_x[N];
    for (n in 1:N)
      if (obs_HIV[n] == 1){
        cs_a_x[n] = cs_a * HIV[n];
        z_cs[n] ~ multi_normal_cholesky(cs_a_x[n], L_Omega_cs);
      } else {
        target += log_mix(p_HIV,
        multi_normal_cholesky_lpdf(z_cs[n] | cs_a, L_Omega_cs),
        multi_normal_cholesky_lpdf(z_cs[n] | [0, 0, 0], L_Omega_cs));
      }
  }
}

generated quantities {
  int imp_clin_symptoms[N];
  for (n in 1:N){
    if (obs_cbm_cs[n]) imp_clin_symptoms[n] = cbm_cs[n];
    else imp_clin_symptoms[n] =  bernoulli_rng(clin_sym[n]);
  }
}
