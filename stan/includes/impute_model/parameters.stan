// * Impute HIV Xd[1] -------------------------------------------------------
real HIV_a0;
vector[2] HIV_a; //WHITE and LYM

// Impute symptoms Td[1:3] for components variables -------------------------
  // Xd[2] is the combined one aka clin_symptoms
vector[3] cs_a0;
matrix[3,2] cs_a; //For HIV & illness days
cholesky_factor_corr[3] L_Omega_cs;
vector<lower=0>[N_pos_cs] z_pos_cs;
vector<upper=0>[N_neg_cs] z_neg_cs;
vector[N_miss_cs] z_miss_cs;

// Impute motor palsy Td[4:6] for components variables ----------------------
  // Xd[3] is the combined one aka clin_motor_palsy
vector[3] mp_a0;
matrix[3, 2] mp_a; //For HIV & illness days
cholesky_factor_corr[3] L_Omega_mp;
vector<lower=0>[N_pos_mp] z_pos_mp;
vector<upper=0>[N_neg_mp] z_neg_mp;
vector[N_miss_mp] z_miss_mp;

// Impute age Xc[1] ---------------------------------------------------------
real age_a0;
real age_a; //for HIV
real<lower=0> age_sigma;
real age_imp[N - sum(obs_Xc[:,1])];

// Impute illness day Xc[2]
real id_a0;
real id_a; //for HIV
real<lower=0> id_sigma;
real id_imp[N - sum(obs_Xc[:,2])];

// Impute csf lab tests as mvNormal -----------------------------------------
vector[5] csf_a0;
vector[2] glu_a; //For diabetes
cholesky_factor_corr[5] L_Omega_csf;
vector<lower=0>[5] L_sigma_csf;
real bld_glu_imp[N - obs_bld_glu]; //lower = 1 for numerical stability
real<lower=0> csf_glu_imp[N - obs_csf_glu]; 
real csf_other_imp[(N*3) - obs_csf_other];

// Impute gcs Xc[8]
// Use auxillary gcs components: Tc[,1:3] GCSE, GCSM, GCSV
// GCSV is missing.
vector<lower=0>[3] gcs_a0;
vector[3] gcs_a; //For age
cholesky_factor_corr[3] L_Omega_gcs;
vector<lower=0>[3] L_sigma_gcs;
//Different ranges for different GCS compartment
real<lower=0, upper=4> gcsv_imp[N - sum(obs_Tc[,3])]; 
real<lower=0, upper=5> gcsm_imp[N - sum(obs_Tc[,2])];
real<lower=0, upper=3> gcse_imp[N - sum(obs_Tc[,1])];
