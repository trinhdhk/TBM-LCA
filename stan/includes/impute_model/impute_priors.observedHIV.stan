  if (Td[n,7] == 1) Xd[n,1] ~ bernoulli(p_HIV[n]);
  z_mp[n] ~ multi_normal_cholesky(mp_a0 + mp_a[:,1]*Xd[n,1] + mp_a[:,2]*Xc_imp[n,1] + mp_a[:,3]*obs_test[n] + mp_a[:,4]*test[n], L_Omega_mp);
  z_cs[n] ~ multi_normal_cholesky(cs_a0 + cs_a[:,1]*Xd[n,1] + cs_a[:,2]*Xc_imp[n,1] + cs_a[:,3]*obs_test[n] + cs_a[:,4]*test[n], L_Omega_cs);
  // Xc_imp[n,1] ~ normal(age_a0 + age_a*Xd[n,1], age_sigma);
  Xc_imp[n,1] ~ normal(id_a0  + id_a[1]*Xd[n,1] + id_a[2]*obs_test[n] + id_a[3]*test[n], id_sigma );
  {
    matrix[6,6] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
    Xc_imp[n,2:7] ~ multi_normal_cholesky(csf_a0 + csf_a[:,1]*Xd[n,1] + csf_a[:,2]*obs_test[n] + csf_a[:,3]*test[n], L_Sigma_csf);
  }

  {
    matrix[3,3] L_Sigma_gcs = diag_pre_multiply(L_sigma_gcs, L_Omega_gcs);
    GCS_imp[n,:] ~ multi_normal_cholesky(gcs_a0 + gcs_a[:,1]*Xd[n,1] + gcs_a[:,2]*obs_test[n] + gcs_a[:,3]*test[n], L_Sigma_gcs);
  }
