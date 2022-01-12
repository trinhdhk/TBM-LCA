  if (Td[n,7]) Xd[n,1] ~ bernoulli(p_HIV[n]);
  z_mp[n] ~ multi_normal_cholesky(mp_a0 + mp_a[:,1]*Xd[n,1] + mp_a[:,2]*Xc_imp[n,1], L_Omega_mp);
  z_cs[n] ~ multi_normal_cholesky(cs_a0 + cs_a[:,1]*Xd[n,1] + cs_a[:,2]*Xc_imp[n,1], L_Omega_cs);
  // Xc_imp[n,1] ~ normal(age_a0 + age_a*Xd[n,1], age_sigma);
  Xc_imp[n,1] ~ normal(id_a0  + id_a*Xd[n,1] , id_sigma );
