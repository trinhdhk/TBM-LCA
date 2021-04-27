// Priors of covariates
a0       ~ student_t(nu, 0, 5);
a        ~ student_t(nu, 0, 2.5);

//Random effects covariates
RE    ~ normal(0, 1);
b_RE  ~ student_t(nu, 0, 1);
b_HIV ~ student_t(nu, 0, 1);

//1-Specificity of each test
z_Xpert[1] ~ normal(logit(.005), 1.59);
z_Mgit[1]  ~ normal(logit(.001), 1.  ); //.82
z_Smear[1] ~ normal(logit(.001), 1.  ); //.7

//Sensitivity of each test
z_Xpert[2] ~ normal(0, 1);
z_Mgit[2]  ~ normal(0, 1);
z_Smear[2] ~ normal(0, 1);