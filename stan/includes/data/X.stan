  int<lower=1> nXc;     //Number of continuous X
  int<lower=1> nXd;     //Number of discrete X
  int<lower=1> nTd;     //Number of disc. aux covariates for imputation model
  int<lower=1> nTc;     //Number of cont. aux covariates for imputation model
  
  real Xc_all[N_all, nXc]; //Continuous covariates
  int  Xd_all[N_all, nXd]; //Discrete covariates
  int  Td_all[N_all, nTd]; //Auxillary covariates - discrete
  real Tc_all[N_all, nTc]; //Auxillary covariates - continuous
  
  //Matrices of observation
  int<lower=0, upper=1> obs_Xc_all[N_all, nXc]; 
  int<lower=0, upper=1> obs_Xd_all[N_all, nXd];
  int<lower=0, upper=1> obs_Td_all[N_all, nTd]; 
  int<lower=0, upper=1> obs_Tc_all[N_all, nTc]; 
