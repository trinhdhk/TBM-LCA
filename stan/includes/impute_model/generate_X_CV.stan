{
  // Imputation model --------------------------------------------------------
  matrix[N_all, 3] Xd_imp; //fully imputed discrete X
  vector[nXc] Xc_imp[N_all]; //fully imputed cont X
  real age_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),1])];
  real id_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),2])];
  
  {
    vector[N_all] z_HIV = HIV_a0 + to_matrix(Tc_all[:,4:5]) * HIV_a;
    Xd_imp[:,1] = binary_rng(impute_binary(Xd_all[:,1], obs_Xd_all[:,1], to_array_1d(z_HIV)), obs_Xd_all[:,1]);
  }
  
  //Age & illness day
  {
    int j =1; int k= 1;
    for (n in which_not(keptin)){
      if (obs_Xc_all[n,1] == 0) {
        age_imp_valid[j]  = normal_rng(age_a0  + age_a*(Xd_imp[n,1]) , age_sigma );
        j+=1;
      }
      if (obs_Xc_all[n,2] == 0) {
        id_imp_valid[k]  = normal_rng(id_a0  + id_a*(Xd_imp[n,1]) , id_sigma );
        k+=1;
      }
    }
  }
  
  if (sum(keptin) < N_all){
    Xc_imp[which_not(keptin),1] = impute_cont_1d(Xc_all[which_not(keptin),1], obs_Xc_all[which_not(keptin),1], age_imp_valid); //Age
    Xc_imp[which_not(keptin),2] = impute_cont_1d(Xc_all[which_not(keptin),2], obs_Xc_all[which_not(keptin),2], id_imp_valid); //illness day
  }
  Xc_imp[which(keptin),1] = impute_cont_1d(Xc[:, 1], obs_Xc[:, 1], age_imp);
  Xc_imp[which(keptin),2] = impute_cont_1d(Xc[:, 2], obs_Xc[:, 2], id_imp);
  
  // Clinical symptoms
  {
    int Td_cs_valid[N_valid,3] = Td_all[which_not(keptin),1:3];
    int obs_cs_valid[N_valid,3] = obs_Td_all[which_not(keptin),1:3];
    vector[3] Mu_cs[N_valid];
    int CS_imp_valid[N_valid, 3];
    int j = 1;
    if (sum(keptin) < N_all){
      for (n in which_not(keptin)){
        Mu_cs[j] = cs_a0 + cs_a[:,1]*Xd_imp[n,1] + cs_a[:,2]*Xc_imp[n,2];
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
        Mu_mp[j] = mp_a0 + mp_a[:,1]*Xd_imp[n,1] + mp_a[:,2]*Xc_imp[n,2];
        j += 1;
      }
      mp_imp_valid = multi_probit_partial_rng(Td_mp_valid, obs_mp_valid, Mu_mp, L_Omega_mp);
      Xd_imp[which_not(keptin),3] = to_vector(any(mp_imp_valid));
    }
    Xd_imp[which(keptin),3] = binary_rng(impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp), obs_Xd[:,3]); 
  }
  
  // CSF lab
  {
    if (sum(keptin) < N_all){
      vector[6] csf_mu_valid[N_valid] = rep_array(rep_vector(0, 6), N_valid) ;
      Xc_imp[which_not(keptin),3:8] = multi_normal_cholesky_partial_rng(Xc_all[which_not(keptin),3:8], obs_Xc_all[which_not(keptin),3:8], csf_mu_valid, L_Omega_csf);
    }
    
    Xc_imp[which(keptin),3:8] = impute_cont_2d(Xc[:,3:8], obs_Xc[:,3:8], csf_imp);
  }
  
  if (nXc > 8) // Other if exists
  for (j in 9:nXc) Xc_imp[:, j] = Xc_all[:, j];
  
  X = append_col(append_col(Xd_imp, to_matrix(Xd_all[:,4:nXd])), append_all(Xc_imp));
}
