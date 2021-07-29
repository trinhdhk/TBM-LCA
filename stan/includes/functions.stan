/* --------------
Helper functions for model
Trinh Dong, 2021
---------------*/
 
  // Summing over an int array
  int sum2d(int[,] a) {
    int s = 0;
    for (i in 1:size(a))
    s += sum(a[i]);
    return s;
  }
  
  // Summing over an int array, with awareness of missing values
  int sum2d_with_missing(int[,] a, int[,] b){
    int s = 0;
    if (size(a) != size(b)) reject("Size a and b must match!");
    for (i in 1:size(a)){
      for (j in 1:size(a[i]))
      s += (b[i,j] == 1) ? a[i,j] : 0;
    }
    return s;
  }
  
  //Mimick the any function in R
  int[] any(int[,] X){
    int N = dims(X)[1];
    int M = dims(X)[2];
    int S[N];
    
    for (n in 1:N){
      if (min(X[n]) < 0 || max(X[n]) > 1) {
        reject("X must be of type binary!");
        print("Min X: ", min(X[n]));
        print("Max X: ", max(X[n]));
      }
      S[n] = sum(X[n]) > 0 ? 1 : 0;      
    }
    return S;
  }
  
  // Mimick the which function in R
  // might not work for !v, use which_not instead;
  int[] which(int[] v){
    int W[sum(v)];
    int w = 1;
    if ((min(v) < 0) || (max(v) > 1)) reject("v must be of type binary!");
    for (i in 1:size(v)){
      if (v[i] == 1) {
        W[w] = i;
        w += 1;
      }
    }
    return W;
  }
  
  // Equivalent to which(!v) in R
  int[] which_not(int[] v){
    int W[size(v) - sum(v)];
    int w = 1;
    if ((min(v) < 0) || (max(v) > 1)) reject("v must be of type binary!");
    for (i in 1:size(v)){
      if (v[i] == 0) {
        W[w] = i;
        w += 1;
      }
    }
    return W;
  }
  
   //Append all vector in an array to a matrix
  matrix append_all(vector[] x){
    matrix[size(x), num_elements(x[1])] M;
    if (num_elements(x) % size(x)) 
      reject("Problem with vector sizes. They must be of the same length.");
    for (i in 1:size(x)) M[i,:] = (x[i])'; 
    return M;
  }
  
  //Calculate the sum of size of each vector in an array of vector
  int sum_size(vector[] x){
    int s = 0;
    int N = size(x);
    for (n in 1:N) s += num_elements(x[n]);
    return s;
  }
  
  //Concat vectors into one
  vector concat(vector[] x){
    int N = size(x);
    int S = sum_size(x);
    vector[S] y;
    int i = 1;
    
    for (n in 1:N){
      y[i:(i + num_elements(x[n])-1)] = x[n];
      i += num_elements(x[n]);
    }
    
    return y;
  }
  
  //Rep vectors into one vector
  //n in the number of replica
  vector rep_vectors(vector x, int n){
    vector[num_elements(x) * n] y;
    for (i in 1:n) y[(num_elements(x)*(i-1) + 1) : num_elements(x)*i] = x;
    return y;
  }
  
  // Function to sum over 2||3 "or" logical variables.
  // this can be generalise but I am lazy so meh
  real sum_probs_or(real[] probs){
    int s = size(probs);
    real sum_prob;
    
    if (s>3||s<2) reject("Only size 2 and 3 are supported");
    sum_prob = (s==2) ? (sum(probs) - prod(probs)) : (sum(probs) - prod(probs[1:2]) - prod(probs[2:3]) - prod(probs[{1,3}]) + prod(probs));
    
    return sum_prob;
  }
  
  // This function do the power operator but return int rather than real
  int int_power(int x, int n){
    int X = 1;
    
    if (n < 0) reject("n must be natural");
    if (n == 0){
      if (X == 0) reject("0^0 is not allowed!");
      return 1;
    }
    
    for (n_ in 1:n) X *= x;
    return X;
  }
  
  // For logistic regressions with missing binary covariates, we have to marginalise/intergrate over the different mixtures. 
  // Patterns returned are different combination of those covariates.
  // For example of 3 missing predictors: 000, 001, 010, 011, 100, 101, 110, 111
  // Currently only binary is supported.
  vector[] get_patterns(row_vector Xd_imp, int[] obs_Xd_imp, vector a_Xd_imp){
    // X_imp are discrete X with missing,
    int N_Xd_imp = size(obs_Xd_imp);
    int N_miss = N_Xd_imp - sum(obs_Xd_imp);
    int T = int_power(2, N_miss);
    matrix[T, N_Xd_imp] obs_pattern; //Construction vector
    matrix[T, N_Xd_imp] probs_pattern; //Probability vector
    matrix[T, N_Xd_imp] a_pattern; //Coefficent vector
    
    vector[T] pat = rep_vector(0, T);
    vector[T] log_probs = rep_vector(0, T);
    vector[T] out[2];
    // vector[T] val[2];
    int j = 0;
    if (size(obs_Xd_imp) != N_Xd_imp) reject("Size mismatched!");
    
    for (i in 1:N_Xd_imp){
      vector[int_power(2, j+1)] V;
      
      if (obs_Xd_imp[i] == 0){
        V = append_row(rep_vector(0, int_power(2, j)), rep_vector(1, int_power(2, j)));
        obs_pattern[:,i] = rep_vectors(V, T %/% int_power(2, j+1));
        j += 1;
      } else {
        obs_pattern[:,i] = rep_vector(Xd_imp[i], T);
        // This is a trick that only works for binary variables.
        // Basically if X[i] is observed, then that X[i] is either 1 or 0
        // if X[i] == 1, setting obs_patterns all 1 will get probs_pattern[,i] == 1 and a_pattern[,i] = a_X[i], which is what we want
        // if X[i] == 0, settings obs_patterns all with 0 wiil get probs_pattern[,i] == 1 - 0 = 1 and a_pattern[,i] = 0, which is also what we want
      }
      
      for (t in 1:T){
        if (obs_pattern[t,i] == 1){
          probs_pattern[t,i] = Xd_imp[i];
          a_pattern[t,i] = a_Xd_imp[i];
        } else {
          probs_pattern[t,i] = 1 - Xd_imp[i];
          a_pattern[t,i] = 0;
        }
      }
      
      pat += a_pattern[:,i];
      log_probs += log(probs_pattern[:,i]);
    }
    
    out[1] = log_probs;
    out[2] = pat;
    return(out);
  }
  
  // Imputation functions -----------------------------------------------------
  // Impute binary discrete variable. Missing value will be "imputed" by their expected probability
  // z is the logistic value of P
  vector impute_binary(int[] raw, int[] obs, real[] z){
    int N = size(raw);
    vector[N] imp;
    
    if (size(obs)!=N||size(z)!=N) reject("Size mismatched!");
    
    for (n in 1:N) imp[n] = (obs[n] == 1) ? raw[n] : Phi(z[n]);
    
    return imp;
  }
  
  // Impute discrete binary combined variable, i.e. cmb = any(el).
  // Missing value will be "imputed" by their expected probs
  vector impute_binary_cmb(int[] raw_cmb, int[] obs_cmb, real[,] z_el, int[,] obs_el){
    // raw_cmb is the raw binary data
    // obs_cmb is the vector of observation for raw_cmb
    // z_el is the matrix of Logistics(prob) for each elements
    // obs_el is the matrix of observations for each elements
    int N = size(raw_cmb);
    vector[N] cmb;
    
    if (size(obs_cmb)!=N||size(z_el)!=N||size(obs_el)!=N) reject("Size mismatched!");
    
    for (n in 1:N){
      if (obs_cmb[n]) cmb[n] = raw_cmb[n];
      else {
        int n_obs_el = sum(obs_el[n,:]);
        if (n_obs_el == 0) cmb[n] = sum_probs_or(Phi(z_el[n,:]));
        if (n_obs_el == 1) {
          if (obs_el[n,1] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,2:3]));
          if (obs_el[n,2] == 1) cmb[n] = sum_probs_or(Phi(z_el[n,{1,3}]));
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
  // raw is the raw vector with missing
  // obs is the vector of observation for raw
  // rep is the impute value for missing values
  real[] impute_cont_1d(real[] raw, int[] obs, real[] rep){
    int N = size(raw);
    real imp[N];
    int i = 1;
    
    if (size(obs)!=N) reject("Size mismatched!");
    if (size(rep)!=N-sum(obs)) reject("Replacement size mismatched!");
    
    for (n in 1:N){
      if (obs[n]) imp[n] = raw[n];
      else {
        imp[n] = rep[i];
        i += 1;
      }
    }
    return imp;
  }
  
  
  // Impute continuous variables 2d
  // raw is the raw matrix with missing
  // obs is the vector of observation for raw
  // rep is the impute value for missing values
  vector[] impute_cont_2d(real[,] raw, int[,] obs, real[] rep){
    int N = dims(raw)[1];
    int M = dims(raw)[2];
    vector[M] imp[N];
    int k = 1;
    
    if (M <= 0) reject("M must be > 0");
    if (dims(obs)[1]!=N || dims(obs)[2]!=M) reject("Size mismatched!");
    
    for (i in 1:N){
      for (j in 1:M){
        if (obs[i,j] == 1) imp[i,j] = raw[i,j];
        else {
          imp[i,j] = rep[k];
          k += 1;
        }
      }
    }
    return imp;
  }
  
   // This function for rng binary variables with awareness of observervation
  vector binary_rng(vector imputed_1d, int[] obs){
    int N = num_elements(imputed_1d);
    vector[N] val;
    if (size(obs) != N) reject("Size mismatched!");
    
    for (n in 1:N) {
      while(is_nan(val[n])){
        val[n] = (obs[n] == 1) ? imputed_1d[n] : bernoulli_rng(imputed_1d[n]);
        if (is_nan(val[n])) print(imputed_1d[n]);
      }
    }
    
    return val;
  }
  
  // This function for rng binary variables with awarenss of observervation - 2d
  matrix binary_2d_rng(matrix imputed_2d, int[,] obs){
    int M = dims(imputed_2d)[2];
    int N = dims(imputed_2d)[1];
    matrix[N, M] val;
    if (dims(obs)[1] != N || dims(obs)[2] != M) reject("Size mismatched!");
    for (m in 1:M) val[:, m] = binary_rng(imputed_2d[,m], obs[, m]);
    return val;
  }
  
  // Rng a partially missing value from multinormal distribution
  /**
  X: an array of vector to impute
  obs_X: matrix of the same size with X, 1 for observed, 0 for missing
  Mu: vector of mean of X
  L: cholesky decomposition of X
  **/
  
  vector[] multi_normal_partial_rng(real[,] X, int[,] obs_X, vector[] Mu, matrix S){
    int N = dims(X)[1]; //all observation
    int M = dims(X)[2]; //all value
    vector[M] new_X[N]; //fully imputed X 
    
    if (dims(X)[1] == 0) return {[]'};
    if (min(to_matrix(obs_X)) < 0 || max(to_matrix(obs_X)) > 1) reject("obs_X must be positive");
    if (dims(obs_X)[1] != N || dims(obs_X)[2] != M) reject("Size mismatched!");
    if (dims(S)[1] != M || dims(S)[2] != M) reject("Number of elements in L mismatched!");
    if (num_elements(Mu[1]) != M) reject("Number of elements in Mu mismatched!");
    
    for (n in 1:N){
      // vector[M] X_n = X[n];
      int obs[M] = obs_X[n];
      new_X[n] = to_vector(X[n]);
      
      if (sum(obs)<M && sum(obs)>0) { //any missing
        /**
        We assume the missing value in the vector is x,
        and observed values in the vector is y.
        x will follow a conditional (multivariate) normal distribution 
        Normal(E_x, S_x), in which:
        E_x = Mu_x + S * inverse(Syy) * (y - Mu_y);
        S_x = Sxx - S * inverse(Syy) * S';
        Since normal_lpdf receives sigma aka sqrt(S_x) in which S_x is a scalar,
        we will have to divide into 2 conditions: 
        x is a scalar and x is a vector.
        The same applies with y.
        **/
        int N_obs  = sum(obs);
        int N_miss = M - N_obs;
        
        int id_obs[N_obs]          = which(obs);
        matrix[N_obs, N_obs] Syy   = S[id_obs, id_obs]; 
        vector[N_obs] Mu_y         = Mu[n,id_obs];
        vector[N_obs] y            = to_vector(X[n,id_obs]);
        
        int id_miss[N_miss]        = which_not(obs);
        matrix[N_miss, N_miss] Sxx = S[id_miss, id_miss];
        vector[N_miss] Mu_x        = Mu[n,id_miss];
        
        matrix[N_miss, N_obs] Sxy  = S[id_miss, id_obs];
        
        if (N_miss == 1){
          real x[1]         = normal_rng(Mu_x + Sxy*(inverse(Syy)*(y - Mu_y)),
                                         sqrt(Sxx - Sxy*(inverse(Syy)*(Sxy')))[1,1]); //[1,1] otherwise will be recognised as matrices
          new_X[n, id_miss] = to_vector(x);
        } else {
          vector[N_miss] x  = multi_normal_rng(Mu_x + Sxy*(inverse(Syy)*(y - Mu_y)),
                                              Sxx - Sxy*(inverse(Syy)*(Sxy')));
          new_X[n, id_miss] = x;                                              
        }
      } else if (sum(obs) == 0){ //Nothing is observed then normal rng is performed
        vector[M] Mu_x = Mu[n];
        new_X[n] = multi_normal_rng(Mu_x, S);
      }
    }
    
    return new_X;
  } 
  
  vector[] multi_normal_cholesky_partial_rng(real[,] X, int[,] obs_X, vector[] Mu, matrix L){
    int M = dims(X)[2]; //all value
    matrix[M, M] S = tcrossprod(L); // reconstruct the cov matrix from cholesky
    if (dims(X)[1] == 0) return {[]'};
    return(multi_normal_partial_rng(X, obs_X, Mu, S));
  }
  
  /**
  Rng the multi probit outcome.
  Because the multi probit 1993 technique only put the constraints on the SUR
  The idea is to rng the observed basing on the constraints
  And use that rng to sampling the missing using the multi_normal_cholesky_partial_rng
  ** Due to the limitation in Stan, X must be constrained before fed into the function.
  **/
  int[,] multi_probit_partial_rng(int[,] X, int[,] obs_X, vector[] Mu, matrix L){
    int empty[0,0];
    int N = dims(X)[1];
    int M = dims(X)[size(dims(X))];
    
    int sign_X[N,M];
    vector[M] rng_X[N];
    int result[N,M];
    if (N == 0) {
      return empty;
    };
    
    if (min(to_matrix(obs_X)) < 0 || max(to_matrix(obs_X)) > 1) reject("obs_X must be positive");
    if (dims(obs_X)[1] != N || dims(obs_X)[2] != M) reject("Size mismatched!");
    if (dims(L)[1] != M || dims(L)[2] != M) reject("Number of elements in L mismatched!");
    if (num_elements(Mu[1]) != M) reject("Number of elements in Mu mismatched!");
    
    for (n in 1:N)
      for (m in 1:M){
        sign_X[n,m] = (obs_X[n,m]==0) ? 0 : ((X[n,m]==1) ? 1 : -1);
      }
    
    for (n in 1:N){
      result[n, which(obs_X[n,:])] = X[n,which(obs_X[n,:])];
      if (sum(obs_X[n]) < M){
        int x[M] = sign_X[n,:];
        int sign_agree = 0;
        int iter = 1;
        while (sign_agree == 0 && iter < 100000){
          int m      = 1;
          rng_X[n]   = multi_normal_cholesky_rng(Mu[n], L);
          sign_agree = 1;
          while (sign_agree == 1 && m <= M){
            if (rng_X[n, m]*x[m] < 0) sign_agree = 0;
            m += 1;
          }
          iter += 1;
        }
        
        for (k in which_not(obs_X[n,:])){
          result[n, k] = bernoulli_rng(Phi(rng_X[n,k]));
        }
      }
    }
    
    return result;
  }
  
  // Function for truncated normal rng only accept pos value
  real[] pos_normal_rng(real[] mu, real sigma) {
    int  N = size(mu);
    real p_0[N];
    real u[N];
    real y[N];
    for (n in 1:N){
      p_0[n] = normal_cdf(0, mu[n], sigma);
      u[n]   = uniform_rng(p_0[n], 1);
      y[n]   = mu[n] + sigma * Phi(u[n]);
    } 
    return y;
  }

