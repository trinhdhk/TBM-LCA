  
  ll_z_mp[1] = multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,1] + mp_a[:,2]*Xc[n,2], L_Omega_mp);
  ll_z_cs[1] = multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,1] + cs_a[:,2]*Xc[n,2], L_Omega_cs);
  ll_Xc_imp_1[1] = normal_lpdf(Xc_imp[n,1] | age_a0 + age_a, age_sigma);
  ll_Xc_imp_2[1] = normal_lpdf(Xc_imp[n,2] | id_a0  + id_a , id_sigma );
