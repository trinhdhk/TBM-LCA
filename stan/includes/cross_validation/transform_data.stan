//Train data
int<lower=1> N = sum(keptin);
int<lower=0> N_valid = N_all - N;
int N_miss_gscv_valid = N_valid - sum(obs_Tc_all[which_not(keptin),3]);

int<lower=0, upper=1> Y_Smear[N] = Y_Smear_all[which(keptin)];
int<lower=0, upper=1> Y_Mgit[N] = Y_Mgit_all[which(keptin)];
int<lower=0, upper=1> Y_Xpert[N] = Y_Xpert_all[which(keptin)];

real Xc[N, nXc] = Xc_all[which(keptin),:];
int  Xd[N, nXd] = Xd_all[which(keptin),:];
real Tc[N, nTc] = Tc_all[which(keptin),:];
int  Td[N, nTd] = Td_all[which(keptin),:];

int<lower=0, upper=1> obs_Xc[N, nXc] = obs_Xc_all[which(keptin),:];
int<lower=0, upper=1> obs_Xd[N, nXd] = obs_Xd_all[which(keptin),:];
int<lower=0, upper=1> obs_Tc[N, nTc] = obs_Tc_all[which(keptin),:];
int<lower=0, upper=1> obs_Td[N, nTd] = obs_Td_all[which(keptin),:];