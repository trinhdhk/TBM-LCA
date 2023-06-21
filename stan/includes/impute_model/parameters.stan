  // * Impute HIV Xd[1] -------------------------------------------------------
  real HIV_a0;
  vector[2] HIV_a; // obs_test and test_pos;
  
  // Impute symptoms Td[1:3] for components variables -------------------------
    // Xd[2] is the combined one aka clin_symptoms
  real cs_a0;
  vector[4] cs_a; //For HIV, illness days, obs_test, test_pos
  //vector<lower=0, upper=1>[2] cs_p12; //True positive rate for first two manifest variables
  vector<lower=0, upper=1>[3] cs_p; //True positive rate for manifest variables
  
  // Impute motor palsy Td[4:6] for components variables ----------------------
    // Xd[3] is the combined one aka clin_motor_palsy
  real mp_a0;
  vector[4] mp_a; //For HIV, illness days, obs_test, test_pos
  // vector<lower=0, upper=1>[2] mp_p12; //True positive rate for first two manifest variables
  vector<lower=0, upper=1>[3] mp_p; //True positive rate for first two manifest variables
  
  // Impute illness day Xc[1]
  real id_a0;
  vector[3] id_a; //for HIV, obs_test, test_pos
  real<lower=0> id_sigma;
  real id_imp[N - sum(obs_Xc[:,1])];
  
  // Impute csf lab tests as mvNormal -----------------------------------------
  cholesky_factor_corr[6] L_Omega_csf;
  vector<lower=0>[6] L_sigma_csf;
  real csf_imp[N*6 - obs_csf];
  vector[6] csf_a0; 
  matrix[6,3] csf_a; // obs_test, test_pos
  
  //Impute GCS
  // vector<lower=0, upper=1>[3] gcs_a0;
  cholesky_factor_corr[3] L_Omega_gcs;
  vector<lower=0>[3] L_sigma_gcs;
  real<lower=0, upper=1> gcs_imp[N*3 - obs_gcs_compartments];
  vector<lower=-1, upper=1>[3] gcs_a0;
  matrix<lower=-1, upper=1>[3, 3] gcs_a; //GCS
  
  //Lymph and WHITE in blood
  // vector<lower=0>[2] bld_a0;
  // cholesky_factor_corr[2] L_Omega_bld;
  // vector<lower=0>[2] L_sigma_bld;
  // real<lower=0> bld_imp[N*2 - obs_bld];
