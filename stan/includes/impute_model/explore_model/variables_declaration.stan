real p_HIV = Phi(HIV_a0); //Populational Probability of having HIV
matrix[N, 7] Xd_imp; //fully imputed discrete X
vector[nXc2 + 1 /* GluRatio */] Xc_imp[N]; //fully imputed cont X
matrix[N, nX + 1 - dims(Xd_imp)[2]] X_compl;
real GCSV_imp[N]; 

Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], rep_array(HIV_a0, N)); //HIV
Xd_imp[:,2] = impute_binary(Td[:,1], obs_cs[:,1], z_cs[:,1]); //Weight
Xd_imp[:,3] = impute_binary(Td[:,2], obs_cs[:,2], z_cs[:,2]); //Sweats
Xd_imp[:,4] = impute_binary(Td[:,3], obs_cs[:,3], z_cs[:,3]); //Cough
Xd_imp[:,5] = impute_binary(Td[:,4], obs_mp[:,1], z_mp[:,1]); //Hemi
Xd_imp[:,6] = impute_binary(Td[:,5], obs_mp[:,2], z_mp[:,2]); //Para
Xd_imp[:,7] = impute_binary(Td[:,6], obs_mp[:,3], z_mp[:,3]); //Tetra
Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
Xc_imp[:,3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp)); //CSF Lab tests
Xc_imp[:,8] = to_array_1d(to_vector(Xc_imp[:,4]) ./ to_vector(Xc_imp[:,3])); //Glucose ratio
// - GCSV
GCSV_imp = impute_cont_1d(Tc[:,3], obs_Tc[:,3], gcsv_imp); //GCSV
Xc_imp[:,9]  = to_array_1d(log2(to_vector(Tc[:,1]) + 1));
Xc_imp[:,10] = to_array_1d(log2(to_vector(Tc[:,2]) + 1));
Xc_imp[:,11] = to_array_1d(log2(to_vector(GCSV_imp) + 1));

if (nXc > 8) // Other if exists
  for (j in 9:nXc) Xc_imp[:, j + 1 - 1 + nTc /*GCS*/] = Xc[:, j];

X_compl = append_col(append_col(to_matrix(Xd[:,4:nXd]), to_vector(Td[:, nTd])), append_all(Xc_imp)); //complement part of X
