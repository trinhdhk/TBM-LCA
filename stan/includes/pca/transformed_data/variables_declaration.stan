  // CSF PCA
  int P_csf = 5;
  int D_csf = nFA;
  int m_csf = D_csf*(P_csf-D_csf) + (D_csf*(D_csf-1))%/%2; 
  real CSF_scaled[N_all, 5];
