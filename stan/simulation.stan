data {
  int<lower=1> N;
  int<lower=0> nX;
  matrix[N,nX] X;
  int<lower=0, upper=1> Y1[N];
  int<lower=0, upper=1> Y2[N];
  int<lower=0, upper=1> Y3[N];
  int<lower=0, upper=1> obs[N];
  int<lower=0, upper=2> good_prior[2]; 
  //0: weak, 1:moderate, 2: strong. First position is for obs_rate, second is for tests
}

parameters {
  ordered[2] z1;
  ordered[2] z2;
  ordered[2] z3;
  ordered[2] zo;
  vector[nX] a;
  real a0;
  real<lower=0> s;
}

model{
  a  ~ normal(0, s);
  a0 ~ normal(0, s);
  s  ~ normal(0, 2.5);
  
  if (good_prior[1]==0){
    zo[1] ~ logistic(0, .7);
    zo[2] ~ logistic(0, .7);
  } else if (good_prior[1]==1) {
    zo[1] ~ logistic(logit(0.01), .7);
    zo[2] ~ logistic(logit(0.95), .7);
  } else {
    zo[1] ~ logistic(logit(0.01), .25);
    zo[2] ~ logistic(logit(0.95), .25);
  }
  
  if (good_prior[2]==0){
    z1    ~ logistic(0,.7);
    z2    ~ logistic(0,.7);
    z3    ~ logistic(0,.7);
  } else if (good_prior[2]==1){
    z1[1] ~ logistic(-2.19,.7);
    z2[1] ~ logistic(-4.59,.7);
    z3[1] ~ logistic(-2.94,.7);
    z1[2] ~ logistic(0,.7);
    z2[2] ~ logistic(0,.7);
    z3[2] ~ logistic(0,.7);
  } else {
    z1[1] ~ logistic(-2.19,.2);
    z2[1] ~ logistic(-4.59,.2);
    z3[1] ~ logistic(-2.94,.2);
    z1[2] ~ logistic(0,.7);
    z2[2] ~ logistic(0,.7);
    z3[2] ~ logistic(0,.7);
  }
  for (n in 1:N){
      real z = a0 + X[n,:]*a;
      real p = inv_logit(z);
      real ll1[2] = 
        obs[n] == 1 ? 
        {bernoulli_logit_lpmf(Y1[n]|z1[1]), bernoulli_logit_lpmf(Y1[n]|z1[2])} :
        {0, 0};
      real ll2[2] = 
        obs[n] == 1 ? 
        {bernoulli_logit_lpmf(Y2[n]|z2[1]), bernoulli_logit_lpmf(Y2[n]|z2[2])} :
        {0, 0};
      real ll3[2] = 
        obs[n] == 1 ? 
        {bernoulli_logit_lpmf(Y3[n]|z3[1]), bernoulli_logit_lpmf(Y3[n]|z3[2])} :
        {0, 0};
      target += log_mix(p,
        bernoulli_logit_lpmf(obs[n]|zo[2]) + 
        ll1[2] + ll2[2] + ll3[2],
        bernoulli_logit_lpmf(obs[n]|zo[1]) + 
        ll1[1] + ll2[1] + ll3[1]);
  }
}
