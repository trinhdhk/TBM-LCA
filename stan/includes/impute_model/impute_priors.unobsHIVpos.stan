  
  {
    theta_cs_hiv[1] = cs_a0 + cs_a[1] + cs_a[2]*Xc_imp[n,1] + cs_a[3]*obs_test[n] + cs_a[4]*test[n];
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
    l1 /= 1 - prod(1 - cs_p);
    if (l2 == 1) {
       if (l1==0) {
        ll_z_cs[1] =  log(1 - inv_logit(theta_cs_hiv[1]));
      } else ll_z_cs[1] = log_mix(inv_logit(theta_cs_hiv[1]), log(l1), 0);
    } else {
      ll_z_cs[1] = log(inv_logit(theta_cs_hiv[1])) + log(l1);
    }
    theta_cs_hiv2[1] = inv_logit(theta_cs_hiv[1]) * l1 / (inv_logit(theta_cs_hiv[1])*l1 + (1-inv_logit(theta_cs_hiv[1]))*l2);
  }
  
  {
    theta_mp_hiv[1] = mp_a0 + mp_a[1] + mp_a[2]*Xc_imp[n,1] + mp_a[3]*obs_test[n] + mp_a[4]*test[n];
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
    l1 /= 1 - prod(1 - mp_p);
    if (l2 == 1) {
       if (l1==0) {
        ll_z_mp[1] = log(1 - inv_logit(theta_mp_hiv[1]));
      } else ll_z_mp[1] = log_mix(inv_logit(theta_mp_hiv[1]), log(l1), 0);
    } else {
      ll_z_mp[1] = log(inv_logit(theta_mp_hiv[1])) + log(l1);
    }
    theta_mp_hiv2[1] = inv_logit(theta_mp_hiv[1]) * l1 / (inv_logit(theta_mp_hiv[1])*l1 + (1-inv_logit(theta_mp_hiv[1]))*l2);
  }
  // ll_z_mp[1] = multi_normal_cholesky_lpdf(z_mp[n] | mp_a0 + mp_a[:,1] + mp_a[:,2]*Xc_imp[n,1] + mp_a[:,3]*obs_test[n] + mp_a[:,4]*test[n], L_Omega_mp);
  // ll_z_cs[1] = multi_normal_cholesky_lpdf(z_cs[n] | cs_a0 + cs_a[:,1] + cs_a[:,2]*Xc_imp[n,1] + cs_a[:,3]*obs_test[n] + cs_a[:,4]*test[n], L_Omega_cs);
  // ll_Xc_imp_1[1] = normal_lpdf(Xc_imp[n,1] | age_a0 + age_a, age_sigma);
  ll_Xc_imp_2[1] = normal_lpdf(Xc_imp[n,1] | id_a0  + id_a[1]  + id_a[2]*obs_test[n] + id_a[3]*test[n], id_sigma );
  {
    matrix[6,6] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
    ll_csf_imp[1] = multi_normal_cholesky_lpdf(Xc_imp[n,2:7] | csf_a0 + csf_a[:,1] + csf_a[:,2]*obs_test[n] + csf_a[:,3]*test[n], L_Sigma_csf);
  }
  {
    if ((Xc[n,8] > -1) || (obs_Xc[n,8] == 0)){
      matrix[3,3] L_Sigma_gcs = diag_pre_multiply(L_sigma_gcs, L_Omega_gcs);
      ll_gcs_imp[1] = multi_normal_cholesky_lpdf(GCS_imp[n,:] | gcs_a0 + gcs_a[:,1] + gcs_a[:,2]*obs_test[n] + gcs_a[:,3]*test[n], L_Sigma_gcs);
    }
  }
