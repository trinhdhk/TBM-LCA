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
// vector[6] csf_a0;
// vector[2] glu_a; //For diabetes //Second thought: this does not work.
cholesky_factor_corr[6] L_Omega_csf;
vector<lower=0>[6] L_sigma_csf;
real csf_imp[N*6 - obs_csf];

