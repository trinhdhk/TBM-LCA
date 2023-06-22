{
  vector[N_all] pY000_theta = inv_logit(-z_Smear_RE) .* inv_logit(-z_Mgit_RE) .* inv_logit(-z_Xpert_RE);
  real pY000_1mtheta = inv_logit(-z_Smear[1]) * inv_logit(-z_Mgit[1]) * inv_logit(-z_Xpert[1]);
  vector[N_all] pY00_theta[3];
  real pY00_1mtheta[3];
  
  pY00_theta[1] = inv_logit(-z_Smear_RE) .* inv_logit(-z_Mgit_RE);
  pY00_theta[2] = inv_logit(-z_Smear_RE) .* inv_logit(-z_Xpert_RE);
  pY00_theta[3] = inv_logit(-z_Mgit_RE) .* inv_logit(-z_Xpert_RE);
  
  pY00_1mtheta[1] = inv_logit(-z_Smear[1]) * inv_logit(-z_Mgit[1]) ;
  pY00_1mtheta[2] = inv_logit(-z_Smear[1]) * inv_logit(-z_Xpert[1]) ;
  pY00_1mtheta[3] = inv_logit(-z_Mgit[1]) * inv_logit(-z_Mgit[1]) ;
  
  Y0_theta[1] = (theta .* inv_logit(-z_Smear_RE)) ./ (1 - p_Smear);
  Y0_theta[2] = (theta .* inv_logit(-z_Mgit_RE)) ./ (1 - p_Mgit);
  Y0_theta[3] = (theta .* inv_logit(-z_Xpert_RE)) ./ (1 - p_Xpert);
  
  for (i in 1:3){
    Y00_theta[i] = (theta .* pY00_theta[i]) ./ (theta .* pY00_theta[i] + (1 - theta) * pY00_1mtheta[i]);
  }
  Y000_theta = (theta .* pY000_theta) ./ (theta .* pY000_theta + (1 - theta) * pY000_1mtheta);
}
