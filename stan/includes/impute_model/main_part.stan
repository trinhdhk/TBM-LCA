real p_HIV = Phi(HIV_a0); //Populational Probability of having HIV
matrix[N, 3] Xd_imp; //fully imputed discrete X
vector[nXc + 1] Xc_imp[N]; //fully imputed cont X
matrix[N, nX + 1 - 3] X_compl;
real GCSV_imp[N]; 

Xd_imp[:,1] = impute_binary(Xd[:,1], obs_Xd[:,1], rep_array(HIV_a0, N)); //HIV
Xd_imp[:,2] = impute_binary_cmb(Xd[:,2], obs_Xd[:,2], to_array_2d(append_all(z_cs)), obs_cs); //Clinical Symptoms
Xd_imp[:,3] = impute_binary_cmb(Xd[:,3], obs_Xd[:,3], to_array_2d(append_all(z_mp)), obs_mp); //Motor palsy
Xc_imp[:,1] = impute_cont_1d(Xc[:,1], obs_Xc[:,1], age_imp); //Age
Xc_imp[:,2] = impute_cont_1d(Xc[:,2], obs_Xc[:,2], id_imp); //Illness day
Xc_imp[:,3:7] = impute_cont_2d(Xc[:,3:7], obs_Xc[:,3:7], append_array(append_array(bld_glu_imp, csf_glu_imp), csf_other_imp)); //CSF Lab tests
Xc_imp[:,8] = to_array_1d(to_vector(Xc_imp[:,4]) ./ to_vector(Xc_imp[:,3])); //Glucose ratio
// - GCSV
GCSV_imp = impute_cont_1d(Tc[:,3], obs_Tc[:,3], gcsv_imp); //GCSV
Xc_imp[:,9] = to_array_1d(log2(to_vector(Tc[:,1]) + to_vector(Tc[:,2]) + to_vector(GCSV_imp) + 1));

if (nXc > 8) // Other if exists
for (j in 9:nXc) Xc_imp[:, j + 1] = Xc[:, j];

X_compl = append_col(to_matrix(Xd[:,4:nXd]), append_all(Xc_imp)); //complement part of X

// Imputation ---------------------------------------------------------------
// - HIV
HIV_a0 ~ normal(0, 2.5);
// - Clinical symptoms
L_Omega_cs ~ lkj_corr_cholesky(4);
cs_a0 ~ normal(0, 2.5);
to_vector(cs_a) ~ normal(0, 2.5);
// - Motor palsy
L_Omega_mp ~ lkj_corr_cholesky(4);
mp_a0 ~ normal(0, 2.5);
to_vector(mp_a) ~ normal(0, 2.5);
// - Age
age_a0 ~ normal(0, 2.5);
age_a  ~ normal(0, 2.5);
age_sigma ~ normal(0, 1);

// - Illness day
id_a0 ~ normal(0, 2.5);
id_a  ~ normal(0, 2.5);
id_sigma ~ normal(0, 1);

// - CSF lab tests  // -GCSV
to_vector(glu_a) ~ normal(0, 2.5);
L_Omega_csf ~ lkj_corr_cholesky(4);
L_sigma_csf ~ normal(0, 2.5);
csf_a0 ~  normal(0, 2.5);

gcsv_a0 ~ normal(0, 2.5);
gcsv_a  ~ normal(0, 2.5);
gcsv_sigma ~ normal(0, 1);

{
  vector[5] csf_mu[N];
  matrix[5, 5] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
  real gcsv_x[N];
  for (n in 1:N){
    csf_mu[n, 1:2] = csf_a0[1:2] + glu_a*Td[n, 7];
    csf_mu[n, 3:5] = csf_a0[3:5];
    gcsv_x[n] = gcsv_a0 + dot_product(gcsv_a,(to_vector(Tc[n, 1:2])));
  }
  
  Xc_imp[:,3:7] ~ multi_normal_cholesky(csf_mu, L_Sigma_csf);
  GCSV_imp ~ normal(gcsv_x, gcsv_sigma);
}