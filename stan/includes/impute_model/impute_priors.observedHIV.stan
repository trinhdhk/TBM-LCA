  if (Td[n,7] == 1) Xd[n,1] ~ bernoulli(p_HIV[n]);
  // z_mp[n] ~ multi_normal_cholesky(mp_a0 + mp_a[:,1]*Xd[n,1] + mp_a[:,2]*Xc_imp[n,1] + mp_a[:,3]*obs_test[n] + mp_a[:,4]*test[n], L_Omega_mp);
  // z_cs[n] ~ multi_normal_cholesky(cs_a0 + cs_a[:,1]*Xd[n,1] + cs_a[:,2]*Xc_imp[n,1] + cs_a[:,3]*obs_test[n] + cs_a[:,4]*test[n], L_Omega_cs);
  
  {
    theta_cs[n] = cs_a0 + cs_a[1]*Xd[n,1] + cs_a[2]*Xc_imp[n,1] + cs_a[3]*obs_test[n] + cs_a[4]*test[n];
    real l1 = 1;
    real l2 = 1;
    for (i in 1:3) {
      if (obs_Td[n,i]==1) {
        l1 *=  exp(bernoulli_lpmf(Td[n,i] | cs_p[i]));
        if (Td[n,i] == 1) l2 = 0;
      }
    }
    if (sum(Td[n,1:3]) == 0) {
      l1 = (sum(obs_Td[n,1:3]) == 3) ? 0 : (l1 - prod(1-cs_p));
    }
    l1 /= (1-prod(1 - cs_p));
    if (l2 == 1) {
       if (l1 == 0) {
        target += log(1 - inv_logit(theta_cs[n]));
      } else target += log_mix(inv_logit(theta_cs[n]), log(l1), 0);
    } else {
      target += log(inv_logit(theta_cs[n])) + log(l1);
    }
   if (obs_Xd[n,2] == 0) Xd_imp[n,2] = inv_logit(theta_cs[n]) * l1 / (inv_logit(theta_cs[n])*l1 + (1-inv_logit(theta_cs[n]))*l2);
  }
  
  {
    theta_mp[n] = mp_a0 + mp_a[1]*Xd[n,1] + mp_a[2]*Xc_imp[n,1] + mp_a[3]*obs_test[n] + mp_a[4]*test[n];
    real l1 = 1;
    real l2 = 1;
    for (i in 4:6) {
      if (obs_Td[n,i]==1) {
        l1 *=  exp(bernoulli_lpmf(Td[n,i] | mp_p[i-3]));
        if (Td[n,i] == 1) l2 = 0;
      }
    }
    if (sum(Td[n,4:6]) == 0) {
      l1 = (sum(obs_Td[n,4:6]) == 3) ? 0 : (l1 - prod(1-mp_p));
    }
    l1 /= 1 -(prod(1 - mp_p));
    if (l2 == 1) {
      if (l1 == 0) {
        target += log(1 - inv_logit(theta_mp[n]));
      } else target += log_mix(inv_logit(theta_mp[n]), log(l1), 0);
    } else {
      target += log(inv_logit(theta_mp[n])) + log(l1);
    }
    
    if (obs_Xd[n,3] == 0) Xd_imp[n,3] = inv_logit(theta_mp[n]) * l1 / (inv_logit(theta_mp[n])*l1 + (1-inv_logit(theta_mp[n]))*l2);
  }
  
  
  Xc_imp[n,1] ~ normal(id_a0  + id_a[1]*Xd[n,1] + id_a[2]*obs_test[n] + id_a[3]*test[n], id_sigma );
  {
    matrix[6,6] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
    Xc_imp[n,2:7] ~ multi_normal_cholesky(csf_a0 + csf_a[:,1]*Xd[n,1] + csf_a[:,2]*obs_test[n] + csf_a[:,3]*test[n], L_Sigma_csf);
  }

  {
    if ((Xc[n,8] > -1) || (obs_Xc[n,8] == 0)){
      matrix[3,3] L_Sigma_gcs = diag_pre_multiply(L_sigma_gcs, L_Omega_gcs);
      GCS_imp[n,:] ~ multi_normal_cholesky(gcs_a0 + gcs_a[:,1]*Xd[n,1] + gcs_a[:,2]*obs_test[n] + gcs_a[:,3]*test[n], L_Sigma_gcs);
    }
  }
