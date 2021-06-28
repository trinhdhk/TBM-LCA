// Priors for model != m0
if (penalty_family == 0) a ~ student_t(nu, 0, SP[1]);
if (penalty_family == 1) a ~ double_exponential(0, SP[1]);
if (penalty_family == 2) a ~ normal(0, SP[1]);
