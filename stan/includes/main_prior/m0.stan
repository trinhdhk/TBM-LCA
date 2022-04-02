// Priors of covariates
a0       ~ student_t(nu, 0, 3);

//1-Specificity of each test
if (unsure_spc == 1){
  z_Xpert[1] ~ logistic(logit(.005),1.1);
  z_Mgit[1]  ~ logistic(logit(.001),1.3); //.82
  z_Smear[1] ~ logistic(logit(.001),1.3); //.7
} else {
  z_Xpert[1] ~ logistic(logit(.005), .7);
  z_Mgit[1]  ~ logistic(logit(.001),1  ); //.82
  z_Smear[1] ~ logistic(logit(.001),1  ); //.7
}

//Sensitivity of each test
z_Xpert[2] ~ logistic(0, .38);
z_Mgit[2]  ~ logistic(0, .38);
z_Smear[2] ~ logistic(0, .38);
