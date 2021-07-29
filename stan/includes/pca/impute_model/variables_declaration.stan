  vector[N] z_HIV = HIV_a0 + to_matrix(Tc[:,4:5]) * HIV_a; //Populational Probability of having HIVi
  vector[N] p_HIV = Phi(z_HIV);
  matrix[N, 3] Xd_imp; //fully imputed discrete X
  vector[nXc] Xc_imp[N]; //fully imputed cont X
  matrix[N, nX - 3 - (5-nFA)] X_compl;

  Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], to_array_1d(z_HIV)); //HIV
  Xd_imp[:,2] = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
  Xd_imp[:,3] = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
  Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
  Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
  Xc_imp[:,3:7] = impute_cont_2d(CSF_scaled[which(keptin),:], obs_Xc[:,3:7], csf_imp); //CSF Lab tests
  
  if (nXc > 7) // Other if exists
    for (j in 8:nXc) Xc_imp[:, j] = Xc[:, j];
  
 
  // Imputation ---------------------------------------------------------------
  // - HIV
    HIV_a0 ~ normal(0, 2.5);
    HIV_a ~ normal(0, 1);
    // - Clinical symptoms
    L_Omega_cs ~ lkj_corr_cholesky(4);
    cs_a0 ~ normal(0, 2.5);
    to_vector(cs_a) ~ normal(0, 1);
    // - Motor palsy
    L_Omega_mp ~ lkj_corr_cholesky(4);
    mp_a0 ~ normal(0, 2.5);
    to_vector(mp_a) ~ normal(0, 1);
    // - Age
    age_a0 ~ normal(0, 2.5);
    age_a  ~ normal(0, 1);
    age_sigma ~ std_normal();
    
    // - Illness day
    id_a0 ~ normal(0, 2.5);
    id_a  ~ normal(0, 1);
    id_sigma ~ std_normal();

    // - CSF lab with pca
    mu_psi_csf ~ std_normal();
    sigma_psi_csf ~ std_normal();
    mu_lt_csf ~ std_normal();
    sigma_lt_csf ~ std_normal();
    L_t0_csf ~ std_normal();
    psi0_csf ~ std_normal();
    L_d_csf ~ std_normal();
    
    {
      matrix[P_csf, P_csf] L = cholesky_decompose(Q_csf);
      matrix[D_csf, D_csf] G = cholesky_decompose(diag_matrix(rep_vector(1, D_csf)) - L_csf'/Q_csf*L_csf); 
      // vector[D_csf] Mu_Ucsf[N];
      vector[D_csf] U_csf[N];
      for (n in 1:N){
        U0_csf[n] ~ std_normal();
        vector[D_csf] mu_Ucsf = to_vector(L_csf'/Q_csf * to_matrix(Xc_imp[n,3:7]));
        U_csf[n] = G * U0_csf[n] + mu_Ucsf;
      }
      
      Xc_imp[:,3:7] ~ multi_normal_cholesky(rep_vector(0, P_csf), L);
   
      if (nXc > 7)
        X_compl = append_col(append_col(append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp[:,1:2])), append_all(U_csf)), append_all(Xc_imp[:,8:])); //complement part of X
      else
        X_compl = append_col(append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp[:,1:2])), append_all(U_csf)); //complement part of X
    }
