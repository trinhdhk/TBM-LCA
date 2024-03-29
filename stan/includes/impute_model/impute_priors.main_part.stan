// Imputation ---------------------------------------------------------------
// - HIV
HIV_a0 ~ normal(0, 2.5);
HIV_a ~ normal(0, 2.5);

// - Clinical symptoms
{
  cs_p ~ beta(10, 10);
  cs_a0 ~ normal(0, 2.5);
  cs_a ~ normal(0, 2.5);
}


// - Motor palsy
{
  mp_p ~ beta(10, 10);
  mp_a0 ~ normal(0, 2.5);
  mp_a ~ normal(0, 2.5);
}

// - Illness day
id_a0 ~ normal(0, 2.5);
id_a  ~ normal(0, 2.5);
id_sigma ~ normal(0, 1);

// - CSF lab tests  
L_Omega_csf ~ lkj_corr_cholesky(4);
L_sigma_csf ~ normal(0, 1);
csf_a0 ~ normal(0, 2.5);
to_vector(csf_a) ~ normal(0,2.5);

// - GCS
gcs_a0      ~ double_exponential(0, 0.05);
to_vector(gcs_a) ~ normal(0, .25);
L_Omega_gcs ~ lkj_corr_cholesky(4);
L_sigma_gcs ~ normal(0, 0.1);



// - Blood
// L_Omega_bld ~ lkj_corr_cholesky(4);
// L_sigma_bld ~ normal(0, 1);

// {
  // matrix[2,2] L_Sigma_bld = diag_pre_multiply(L_sigma_bld, L_Omega_bld);
  // Bld_imp ~ multi_normal_cholesky(bld_a0, L_Sigma_bld);
// }
