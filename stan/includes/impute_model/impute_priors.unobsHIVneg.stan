
  ll_z_mp[2] = multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,2]*Xc[n,2], L_Omega_mp);
  ll_z_cs[2] = multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,2]*Xc[n,2], L_Omega_cs);
  ll_Xc_imp_1[2] = normal_lpdf(Xc_imp[n,1] | age_a0, age_sigma);
  ll_Xc_imp_2[2] = normal_lpdf(Xc_imp[n,2] | id_a0 , id_sigma );
