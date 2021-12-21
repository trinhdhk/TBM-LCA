if (obs_Xd[n,1]){
  if (Td[n,7]) Xd[n,1] ~ bernoulli(p_HIV[n]);
  z_mp[n] ~ multi_normal_cholesky(mp_a0 + mp_a[:,1]*Xd[n,1] + mp_a[:,2]*Xc[n,2], L_Omega_mp);
  z_cs[n] ~ multi_normal_cholesky(cs_a0 + cs_a[:,1]*Xd[n,1] + cs_a[:,2]*Xc[n,2], L_Omega_cs);
  Xc_imp[n,1] ~ normal(age_a0 + age_a*Xd[n,1], age_sigma);
  Xc_imp[n,2] ~ normal(id_a0  + id_a*Xd[n,1] , id_sigma );
} else {
  if (is_nan(
    multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,1] + cs_a[:,2]*Xc[n,2], L_Omega_cs) +
    multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,2]*Xc[n,2], L_Omega_cs) +
    multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,1] + mp_a[:,2]*Xc[n,2], L_Omega_mp) +
    multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + cs_a[:,2]*Xc[n,2], L_Omega_mp)))
    // This is to suppress the program from complaining at the start
    target += not_a_number();
    else {
      target += log_mix(p_HIV[n],
        multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,1] + cs_a[:,2]*Xc[n,2], L_Omega_cs),
        multi_normal_cholesky_lpdf(z_cs[n] | cs_a0             + cs_a[:,2]*Xc[n,2], L_Omega_cs));
      target += log_mix(p_HIV[n],
        multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,1] + mp_a[:,2]*Xc[n,2], L_Omega_mp),
        multi_normal_cholesky_lpdf(z_mp[n] | mp_a0             + mp_a[:,2]*Xc[n,2], L_Omega_mp));
      target += log_mix(p_HIV[n],
        normal_lpdf(Xc_imp[n,1]| age_a0 + age_a, age_sigma),
        normal_lpdf(Xc_imp[n,1]| age_a0        , age_sigma));
      target += log_mix(p_HIV[n],
        normal_lpdf(Xc_imp[n,2]| id_a0  + id_a , id_sigma ),
        normal_lpdf(Xc_imp[n,2]| id_a0         , id_sigma ));     
    }
}
