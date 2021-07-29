  // Priors for model != m0
  if (penalty_family == 0) a_raw ~ student_t(nu, 0, 1);
  if (penalty_family == 1) a_raw ~ double_exponential(0, 1);
  if (penalty_family == 2) a_raw ~ normal(0, 1);
