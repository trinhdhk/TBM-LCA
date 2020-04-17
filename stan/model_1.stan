data {
    int<lower=1> N_pt; // number of observations
    // For logistic regression
    // logit(C) = a0 + X*a
    int<lower=0> N_cov;
    matrix[N_pt, N_cov] X; //covariate design matrix
    
    // For latent class model
    // Currently have no covariate for manifest variables
    // For binary manifest var (mv): P(U[i] == 1|c) = b[c,i] + Z[i]e[i] 
    // where Z[i] is cov design matrix for indicator i and e[i] is param vector for indicator i.
    // For continuous mv: mean[U[i]] = b[c,i] + Z[i]e[i] and U[i] ~ normal(mean[U[i]], s) 
    int<lower=1> N_ind;
    int<lower=0, upper=1> U[N_pt, N_ind]; //indicator design matrix
}

parameters {
    vector[N_pt] logitPc; //logit-probabiliy of latent class
    // real<lower=0> sigma_logitPc;
    // For logistic regression
    real a0; //intercept
    vector[N_cov] a; //slope
    
    // For latent model
    matrix[2, N_ind] b;
    real b_re; //random effect param
    vector[N_pt] re; //random effect value
}

transformed parameters{
    // Latent class probability
    vector<lower=0, upper=1>[N_pt] Pc = exp(logitPc)./(1+exp(logitPc));
    matrix<lower=0, upper=1>[N_pt,2] C = [(1-Pc)', Pc']';
    // matrix[N_pt,2] RE = [rep_vector(0, N_pt)', re']'; 
}

model {
    // Latent class model sampling
    b[1,] ~ normal(0,.5);
    b[2,1] ~ normal(2.886, 7.923);//to be change. for xpert.
    b[2,2] ~ normal(3.023, 7.923);//to be change. for culture.
    b[2,3] ~ normal(0,1);
    b_re ~ uniform(0,5);
    
    re ~ normal(0,.5);
    for (i in 1:N_ind){
        U[:,i] ~ bernoulli_logit(C*(b[:,i]) + b_re * C[:,2] .* re);
    }
    
    // Logistic regression sampling
    a0 ~ normal(0,1);
    a[1:N_cov-1] ~ normal(0,.1);
    a[N_cov] ~ normal(2.5, 1); //for xray
    // sigma_logitPc ~ normal(0, 5);
    // for (n in 1:N_pt){
    //     logitPc[n] ~ normal(a0 + X[n]*a, 10);
    // }
    
    logitPc ~ normal(a0 + X*a, 0.5);
}
