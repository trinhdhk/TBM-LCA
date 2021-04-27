int N_miss_gscv_valid = N_valid - sum(obs_Tc_all[which_not(keptin),3]);

real Xc[N, nXc] = Xc_all[which(keptin),:];
int  Xd[N, nXd] = Xd_all[which(keptin),:];
real Tc[N, nTc] = Tc_all[which(keptin),:];
int  Td[N, nTd] = Td_all[which(keptin),:];

int<lower=0, upper=1> obs_Xc[N, nXc] = obs_Xc_all[which(keptin),:];
int<lower=0, upper=1> obs_Xd[N, nXd] = obs_Xd_all[which(keptin),:];
int<lower=0, upper=1> obs_Tc[N, nTc] = obs_Tc_all[which(keptin),:];
int<lower=0, upper=1> obs_Td[N, nTd] = obs_Td_all[which(keptin),:];