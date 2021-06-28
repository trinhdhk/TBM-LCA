//Random effects covariates
RE    ~ normal(0, 1);
if (penalty_family == 0){
  b_RE  ~ student_t(nu, 0, SP[2]);
  b_HIV ~ student_t(nu, 0, SP[2]);
}
if (penalty_family == 1){
  b_RE  ~ double_exponential(0, SP[2]);
  b_HIV ~ double_exponential(0, SP[2]);
}
if (penalty_family == 2){
  b_RE  ~ normal(0, SP[2]);
  b_HIV ~ normal(0, SP[2]);
}

