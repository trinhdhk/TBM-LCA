// Priors for model != m0
a        ~ student_t(nu, 0, 2.5);

// Random effects covariates
b_HIV    ~ normal(0, 1);