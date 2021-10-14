// vector[N] z_HIV = HIV_a0 + to_matrix(Tc[:,4:5]) * HIV_a; //Populational Probability of having HIVi
vector[N] z_HIV;
vector[N] p_HIV;
matrix[N, 3] Xd_imp; //fully imputed discrete X
vector[nXc+nQ] Xc_imp[N]; //fully imputed cont X
matrix[N, nA - 3] X_compl;
vector[3] GCS_imp[N];
vector[2] Bld_imp[N];
vector[nXc+nQ] sd_X;

Bld_imp       = impute_cont_2d(Tc[:,4:5], obs_Tc[:,4:5], bld_imp);
z_HIV         = HIV_a0 + append_all(Bld_imp)*HIV_a;
p_HIV         = Phi(z_HIV);
GCS_imp       = impute_cont_2d(Tc[:,1:3], obs_Tc[:,1:3], gcs_imp);

Xd_imp[:,1]   = impute_binary(Xd[:,1], obs_Xd[:,1], to_array_1d(z_HIV)); //HIV
Xd_imp[:,2]   = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
Xd_imp[:,3]   = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
Xc_imp[:,1]   = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
Xc_imp[:,2]   = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
Xc_imp[:,3:8] = impute_cont_2d(Xc[:,3:8], obs_Xc[:,3:8], csf_imp); //CSF Lab tests
for (n in 1:N){
  Xc_imp[n,9] = obs_Xc[n,9] ? Xc[n,9] : to_array_1d(to_vector(GCS_imp[:,1])*3 + to_vector(GCS_imp[:,2])*5 + to_vector(GCS_imp[:,3])*4)[n];
}

if (nXc > 9) // Other if exists
for (j in 10:nXc) Xc_imp[:,j] = Xc[:,j];

{
  int j = nXc + 1;
  for (q in quad_Xc_idx){
    Xc_imp[:,j] = square(Xc_imp[:,q]);
    j += 1;
  }
}

for (ii in 1:(nXc+nQ)){
  sd_X[ii] = sd(Xc_imp[:,ii]);
}

X_compl = append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp)); //complement part of X
