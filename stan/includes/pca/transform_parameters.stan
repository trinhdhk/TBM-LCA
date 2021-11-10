// PCA ------------------------------------------------------
  vector[m_csf] L_t_csf = L_t0_csf * sigma_lt_csf + mu_lt_csf;
  real<lower=0> psi_csf = psi0_csf * sigma_psi_csf + mu_psi_csf; 
  cholesky_factor_cov[P_csf,D_csf] L_csf;  //lower triangular factor loadings Matrix 
  matrix[P_csf, P_csf] Q_csf;   //Covariance mat
  
  {
    int idx2 = 0;
    for(i in 1:(D_csf-1)){
      for(j in (i+1):D_csf){
        L_csf[i,j] = 0; //constrain the upper triangular elements to zero 
      }
    }
    for (j in 1:D_csf) {
        L_csf[j,j] = L_d_csf[j];
      for (i in (j+1):P_csf) {
        idx2 += 1;
        L_csf[i,j] = L_t_csf[idx2];
      } 
    }
  } 
  
  Q_csf = L_csf*L_csf' + diag_matrix(rep_vector(psi_csf + 1e-14, P_csf)); 
 //----------------------------------------------------------- 
