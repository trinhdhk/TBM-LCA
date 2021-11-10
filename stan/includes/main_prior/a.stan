  // Priors for model != m0
  // if (penalty_family == 0) a_raw ~ student_t(nu, 0, 1);
  // if (penalty_family == 1) a_raw ~ double_exponential(0, 1);
  // if (penalty_family == 2) a_raw ~ normal(0, 1);
  
  if (penalty_family == 0) {
    // for (i in 1:nA) a_raw[i] ~ student_t(nu, 0, inv(sd_X[i]));
    a_raw[1:nXd] ~ student_t(nu, 0, 2);
    for (i in 1:(nXc+nQ)){
      a_raw[nXd+i] ~ student_t(nu, 0, inv(sd_X[i]));
    }
  }
  
  if (penalty_family == 1) {
    // for (i in 1:nA) a_raw[i] ~ double_exponential(0, inv(sd_X[i]));
    a_raw[1:nXd] ~ double_exponential(0, 2);
    for (i in 1:(nXc+nQ)){
      a_raw[nXd+i] ~ double_exponential(0, inv(sd_X[i]));
    }
  }
  
   if (penalty_family == 2) {
     // for (i in 1:nA) a_raw[i] ~ normal(0, inv(sd_X[i]));
     a_raw[1:nXd] ~ normal(0, 2);
    for (i in 1:(nXc+nQ)){
      a_raw[nXd+i] ~ normal(0, inv(sd_X[i]));
    }
  }
