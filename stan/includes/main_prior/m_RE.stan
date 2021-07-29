//Random effects covariates
RE    ~ normal(0, 1);
if (penalty_family == 0){
  b_RE_raw  ~ student_t(nu, 0, 1);
  b_HIV_raw ~ student_t(nu, 0, 1);
  b ~ student_t(nu, 0, 1);
}
if (penalty_family == 1){
  b_RE_raw  ~ double_exponential(0, 1);
  b_HIV ~ double_exponential(0, 1);
  b ~ double_exponential(0, 1);
}
if (penalty_family == 2){
  b_RE_raw  ~ normal(0, 1);
  b_HIV_raw ~ normal(0, 1);
  b ~ normal(0, 1);
}

