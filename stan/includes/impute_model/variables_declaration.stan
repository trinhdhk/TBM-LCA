vector[N] z_HIV = HIV_a0 + to_matrix(Tc[:,4:5]) * HIV_a; //Populational Probability of having HIVi
vector[N] p_HIV = Phi(z_HIV);
matrix[N, 3] Xd_imp; //fully imputed discrete X
vector[nXc] Xc_imp[N]; //fully imputed cont X
matrix[N, nX] X_compl;

Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], to_array_1d(z_HIV)); //HIV
Xd_imp[:,2] = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
Xd_imp[:,3] = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
Xc_imp[:,3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp)); //CSF Lab tests

if (nXc > 7) // Other if exists
for (j in 8:nXc) Xc_imp[:, j] = Xc[:, j];

X_compl = append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp)); //complement part of X
