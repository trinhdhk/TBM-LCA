// Imputation ---------------------------------------------------------------
// - HIV
HIV_a0 ~ normal(0, 2.5);
HIV_a ~ normal(0, 2.5);
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

//GCS
to_vector(gcs_a) ~ normal(0, 2.5);
L_Omega_gcs ~ lkj_corr_cholesky(4);
L_sigma_csf ~ normal(0, 2.5);
gcs_a0 ~ normal(0, 2.5);

/* old code
gcsv_a0 ~ normal(0, 2.5);
gcsv_a  ~ normal(0, 2.5);
gcsv_sigma ~ normal(0, 1);
*/


{
  vector[5] csf_mu[N];
  matrix[5, 5] L_Sigma_csf = diag_pre_multiply(L_sigma_csf, L_Omega_csf);
  
  vector[3] gcs_mu[N];
  matrix[3, 3] L_Sigma_gcs = diag_pre_multiply(L_sigma_gcs, L_Omega_gcs);
  //old code: real gcsv_x[N];
  
  for (n in 1:N){
    csf_mu[n, 1:2] = csf_a0[1:2] + glu_a*Td[n, 7];
    csf_mu[n, 3:5] = csf_a0[3:5];
    
    gcs_mu[n, :] = gcs_a0 + gcs_a*Xc_imp[n,1];
    //gcsv_x[n] = gcsv_a0 + dot_product(gcsv_a,(to_vector(Tc[n, 1:2])));
  }
  
  Xc_imp[:,3:7] ~ multi_normal_cholesky(csf_mu, L_Sigma_csf);
  GCS_imp ~ multi_normal_cholesky(gcs_mu, L_Sigma_gcs);
  //old code GCSV_imp ~ normal(gcsv_x, gcsv_sigma);
}