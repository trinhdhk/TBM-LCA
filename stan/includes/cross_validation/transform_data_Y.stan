int<lower=1> N = sum(keptin);
int<lower=0> N_valid = N_all - N;

int<lower=0, upper=1> Y_Smear[N] = Y_Smear_all[which(keptin)];
int<lower=0, upper=1> Y_Mgit[N] = Y_Mgit_all[which(keptin)];
int<lower=0, upper=1> Y_Xpert[N] = Y_Xpert_all[which(keptin)];
