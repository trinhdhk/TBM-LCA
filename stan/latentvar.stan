functions{
  matrix append_all(vector[] x){
    int M = size(x);
    int N = num_elements(x[1]);
    matrix[M, N] mat;
    for (m in 1:M) mat[m,:] = x[m]';
    return mat;
  }
}
data{
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  matrix[N,M] X;
}

/*
transformed data{
  vector[M] Xcenter[N];
  for (m in 1:M) Xcenter[:,m] = to_array_1d(X[:,m] - mean(X[:,m]));
}*/

parameters{
  ordered[K] Y[N];
  matrix<lower=0, upper=1>[M, K] W;
  real<lower=0> s;
  // vector[K] Mu;
  // cholesky_factor_corr[M] L_Omega;
  // vector<lower=0>[M] L_sigma;
}


model{
  // matrix[M, M] L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  s ~ normal(0, 1);
  to_vector(W) ~ normal(0,1);
  // Mu ~ normal(0, 1);
  // L_Omega ~ lkj_corr_cholesky(4);
  // L_sigma ~ normal(0,1);
  
  for (n in 1:N){
    Y[n] ~ student_t(4, 0, 1);
    X[n] ~ normal(W * (Y[n]), s);
    // Xcenter[n] ~ normal(W * Y[n], s);
    // Xcenter[n] ~ multi_normal_cholesky(W * Y[n], L_Sigma);
  }
}

generated quantities{
  vector[M] X_mu[N];
  vector[K] Y2[N]; 
  vector[N] log_lik = rep_vector(0, N);
  {
    matrix[K, M] invW;
    // matrix[M, M] L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (n in 1:N){
      X_mu[n] = W * Y[n];
      {
        matrix[M, K] Q = qr_thin_Q(W);
        matrix[K, K] R = qr_thin_R(W);
        invW = R \ Q';
      }
      Y2[n] = invW * X_mu[n];
      log_lik[n] += normal_lpdf(X[n] | W * (Y[n]), s);
      // log_lik[n] += multi_normal_cholesky_lpdf(Xcenter[n] | W * Y[n], L_Sigma);
    } 
  }
}
