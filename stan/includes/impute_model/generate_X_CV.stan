{
  // Imputation model --------------------------------------------------------
  matrix[N_all, 3] Xd_imp; //fully imputed discrete X
  vector[nXc + 1] Xc_imp[N_all]; //fully imputed cont X
  real age_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),1])];
  real id_imp_valid[N_valid - sum(obs_Xc_all[which_not(keptin),2])];
  
  
  Xd_imp[:,1] = binary_rng(impute_binary(Xd_all[:,1], obs_Xd_all[:,1], rep_array(HIV_a0, N_all)), obs_Xd_all[:,1]); //HIV
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
  
  Xc_imp[which_not(keptin),1] = impute_cont_1d(Xc_all[which_not(keptin),1], obs_Xc_all[which_not(keptin),1], age_imp_valid); //Age
  Xc_imp[which_not(keptin),2] = impute_cont_1d(Xc_all[which_not(keptin),2], obs_Xc_all[which_not(keptin),2], id_imp_valid); //illness day
  Xc_imp[which(keptin),1] = impute_cont_1d(Xc[:, 1], obs_Xc[:, 1], age_imp);
  Xc_imp[which(keptin),2] = impute_cont_1d(Xc[:, 2], obs_Xc[:, 2], id_imp);
  
  {
    int Td_cs_valid[N_valid,3] = Td_all[which_not(keptin),1:3];
    int obs_cs_valid[N_valid,3] = obs_Td_all[which_not(keptin),1:3];
    vector[3] Mu_cs[N_valid];
    int CS_imp_valid[N_valid, 3];
    int j = 1;
    for (n in which_not(keptin)){
      Mu_cs[j] = cs_a0 + cs_a[:,1]*Xd_imp[n,1] + cs_a[:,2]*Xc_imp[n,2];
      j += 1;
    }
    CS_imp_valid = multi_probit_partial_rng(Td_cs_valid, obs_cs_valid, Mu_cs, L_Omega_cs);
    Xd_imp[which(keptin),2] = binary_rng(impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs), obs_Xd[:,2]); //Clinical Symptoms
    Xd_imp[which_not(keptin),2] = to_vector(any(CS_imp_valid));
  }
  
  
  {
    int Td_mp_valid[N_valid,3] = Td_all[which_not(keptin),4:6];
    int obs_mp_valid[N_valid,3] = obs_Td_all[which_not(keptin),4:6];
    vector[3] Mu_mp[N_valid];
    int mp_imp_valid[N_valid, 3];
    int j = 1;
    for (n in which_not(keptin)){
      Mu_mp[j] = mp_a0 + mp_a[:,1]*Xd_imp[n,1] + mp_a[:,2]*Xc_imp[n,2];
      j += 1;
    }
    mp_imp_valid = multi_probit_partial_rng(Td_mp_valid, obs_mp_valid, Mu_mp, L_Omega_mp);
    Xd_imp[which(keptin),3] = binary_rng(impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp), obs_Xd[:,3]); //Motor palsy
    Xd_imp[which_not(keptin),3] = to_vector(any(mp_imp_valid));
  }
  
  {
    vector[5] csf_mu_valid[N_valid];
    int j = 1;
    for (n in which_not(keptin)){
      csf_mu_valid[j, 1:2] = csf_a0[1:2] + glu_a*Td_all[n, 7];
      csf_mu_valid[j, 3:5] = csf_a0[3:5];
      j += 1;
    }
    Xc_imp[which_not(keptin),3:7] = multi_normal_cholesky_partial_rng(Xc_all[which_not(keptin),3:7], obs_Xc_all[which_not(keptin),3:7], csf_mu_valid, L_Omega_csf);
    Xc_imp[which(keptin),3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp));
  }
  
  Xc_imp[:,8] = to_array_1d(to_vector(Xc_imp[:,4]) ./ to_vector(Xc_imp[:,3])); //Glucose ratio
  
  {
    real GCSV_imp_all[N_all]; // - GCSV
    real gcsv_imp_valid[N_miss_gscv_valid];
    if (N_miss_gscv_valid > 0){
      real gcsv_x_valid[N_miss_gscv_valid];
      int j = 1;
      for (n in which_not(obs_Tc_all[which_not(keptin),3])){
        gcsv_x_valid[j] = gcsv_a0 + dot_product(gcsv_a,(to_vector(Tc_all[n, 1:2])));
        j += 1;
      }
      
      gcsv_imp_valid = to_array_1d(half_normal_rng(gcsv_x_valid, gcsv_sigma));
    }
    
    GCSV_imp_all[which_not(keptin)] = impute_cont_1d(Tc_all[which_not(keptin),3], obs_Tc_all[which_not(keptin),3], gcsv_imp_valid); 
    GCSV_imp_all[which(keptin)] = impute_cont_1d(Tc[:,3], obs_Tc[:,3], gcsv_imp);
    Xc_imp[:,9] = to_array_1d(log2(to_vector(Tc_all[:,1]) + to_vector(Tc_all[:,2]) + to_vector(GCSV_imp_all) + 1));
  }
  
  if (nXc > 8) // Other if exists
  for (j in 9:nXc) Xc_imp[:, j + 1] = Xc_all[:, j];
  
  X = append_col(append_col(Xd_imp, to_matrix(Xd_all[:,4:nXd])), append_all(Xc_imp));
}