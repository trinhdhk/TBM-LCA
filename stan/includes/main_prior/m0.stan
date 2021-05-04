// Priors of covariates
a0       ~ student_t(nu, 0, 2.5);

//1-Specificity of each test
z_Xpert[1] ~ normal(logit(.005), 1.59);
z_Mgit[1]  ~ normal(logit(.001),  .5 ); //.82
z_Smear[1] ~ normal(logit(.001),  .5 ); //.7

//Sensitivity of each test
z_Xpert[2] ~ normal(0, .7);
z_Mgit[2]  ~ normal(0, .7);
z_Smear[2] ~ normal(0, .7);
