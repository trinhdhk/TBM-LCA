{
  // Imputation model --------------------------------------------------------
  matrix[N_all, 3] Xd_imp; //fully imputed discrete X
  vector[nXc+nQ] Xc_imp[N_all]; //fully imputed cont X
  // real age_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),1])];
  real id_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),1])];
  
  //HIV and Blood
  {
    real z_HIV[N_all];
    for (n in 1:N_all){
      if (Td_all[n,7]==0) {
        if (obs_Td_all[n,8]==1){
          z_HIV[n] = (Td_all[n,8]==1) ? inv_Phi(0.0043) : inv_Phi(0.0021);
        } else {
          z_HIV[n] = inv_Phi(0.0032);
        }
      } else{
        z_HIV[n] = HIV_a0;
      }
    }
    Xd_imp[:,1] = binary_rng(impute_binary(Xd_all[:,1], obs_Xd_all[:,1], z_HIV), obs_Xd_all[:,1]);
  }
  
  //Age & illness day
  {
    int j =1; int k= 1;
    for (n in which_not(keptin)){
      // if (obs_Xc_all[n,1] == 0) {
      //   age_imp_valid[j]  = normal_rng(age_a0  + age_a*(Xd_imp[n,1]) , age_sigma );
      //   j+=1;
      // }
      if (obs_Xc_all[n,1] == 0) {
        id_imp_valid[k]  = normal_rng(id_a0  + id_a*(Xd_imp[n,1]) , id_sigma );
        k+=1;
      }
    }
  }
  
  if (sum(keptin) < N_all){
    // Xc_imp[which_not(keptin),1] = impute_cont_1d(Xc_all[which_not(keptin),1], obs_Xc_all[which_not(keptin),1], age_imp_valid); //Age
    Xc_imp[which_not(keptin),1] = impute_cont_1d(Xc_all[which_not(keptin),1], obs_Xc_all[which_not(keptin),1], id_imp_valid); //illness day
  }
  // Xc_imp[which(keptin),1] = impute_cont_1d(Xc[:, 1], obs_Xc[:, 1], age_imp);
  Xc_imp[which(keptin),1] = impute_cont_1d(Xc[:, 1], obs_Xc[:, 1], id_imp);
  
  // Clinical symptoms
  {
    int Td_cs_valid[N_valid,3] = Td_all[which_not(keptin),1:3];
    int obs_cs_valid[N_valid,3] = obs_Td_all[which_not(keptin),1:3];
    vector[3] Mu_cs[N_valid];
    int CS_imp_valid[N_valid, 3];
    int j = 1;
    if (sum(keptin) < N_all){
      for (n in which_not(keptin)){
        Mu_cs[j] = cs_a0 + cs_a[:,1]*Xd_imp[n,1] + cs_a[:,2]*Xc_imp[n,1];
        j += 1;
      }
      CS_imp_valid = multi_probit_partial_rng(Td_cs_valid, obs_cs_valid, Mu_cs, L_Omega_cs);
      Xd_imp[which_not(keptin),2] = to_vector(any(CS_imp_valid));
    }
    Xd_imp[which(keptin),2] = binary_rng(impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs), obs_Xd[:,2]); //Clinical Symptoms
  }
  
  // Motor palsy
  {
    int Td_mp_valid[N_valid,3] = Td_all[which_not(keptin),4:6];
    int obs_mp_valid[N_valid,3] = obs_Td_all[which_not(keptin),4:6];
    vector[3] Mu_mp[N_valid];
    int mp_imp_valid[N_valid, 3];
    int j = 1;
    if (sum(keptin) < N_all){
      for (n in which_not(keptin)){
        Mu_mp[j] = mp_a0 + mp_a[:,1]*Xd_imp[n,1] + mp_a[:,2]*Xc_imp[n,1];
        j += 1;
      }
      mp_imp_valid = multi_probit_partial_rng(Td_mp_valid, obs_mp_valid, Mu_mp, L_Omega_mp);
      Xd_imp[which_not(keptin),3] = to_vector(any(mp_imp_valid));
    }
    Xd_imp[which(keptin),3] = binary_rng(impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp), obs_Xd[:,3]); 
  }
  
  // CSF lab
  {
    matrix[6,6] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
    if (sum(keptin) < N_all){
      vector[6] csf_mu_valid[N_valid] = rep_array(rep_vector(0, 6), N_valid) ;
      Xc_imp[which_not(keptin),2:7] = multi_normal_cholesky_partial_rng(Xc_all[which_not(keptin),2:7], obs_Xc_all[which_not(keptin),2:7], csf_mu_valid, L_Sigma_csf);
    }
    
    Xc_imp[which(keptin),2:7] = impute_cont_2d(Xc[:,2:7], obs_Xc[:,2:7], csf_imp);
  }
  
  // GCS
  {
    vector[3] GCS_imp[N_all];
    matrix[3,3] L_Sigma_gcs = diag_pre_multiply(L_sigma_gcs, L_Omega_gcs);
    if (sum(keptin) < N_all){
      vector[3] gcs_mu_valid[N_valid] = rep_array(gcs_a0, N_valid) ;
      GCS_imp[which_not(keptin),:] = multi_normal_cholesky_partial_unit_rng(Tc_all[which_not(keptin),1:3], obs_Tc_all[which_not(keptin),1:3], gcs_mu_valid, L_Sigma_gcs);
    }
    GCS_imp[which(keptin),:] = impute_cont_2d(Tc[:,1:3], obs_Tc[:,1:3], gcs_imp);
    for (n in 1:N_all){
      real tmp;
      tmp = obs_Xc_all[n,8] ? (Xc_all[n,8]*3) : (GCS_imp[n,1]*3 + GCS_imp[n,2]*5 + GCS_imp[n,3]*4 - 3);
      // Xc_imp[n,9] = obs_Xc_all[n,9] ? Xc_all[n,9] : to_array_1d(to_vector(GCS_imp[:,1])*3 + to_vector(GCS_imp[:,2])*5 + to_vector(GCS_imp[:,3])*4)[n];
      Xc_imp[n,8] = round(tmp)/3;
    }
  }
  
  if (nXc > 8) // Other if exists
    for (j in 9:nXc) Xc_imp[:, j] = Xc_all[:, j];
    
  {
    int j = nXc + 1;
    for (q in quad_Xc_idx){
      Xc_imp[:,j] = square(Xc_imp[:,q]);
      j += 1;
    }
  }
  
  X = append_col(append_col(Xd_imp, to_matrix(Xd_all[:,4:nXd])), append_all(Xc_imp));
}
