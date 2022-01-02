  if (nB > 0){
    if (penalty_family == 1){
      int j = 1;
      for (i in B){
        // b_raw[j]  ~ student_t(nu, 1, inv(sd_X[i+nXc]));
        b_raw[j]  ~ student_t(nu, 0, inv(sd_X[i]));
        j += 1;
      }
    }
    
    if (penalty_family == 1){
      int j = 1;
      for (i in B){
        // b_raw[j] ~ double_exponential(0, inv(sd_X[i+nXc]));
        b_raw[j] ~ double_exponential(0, inv(sd_X[i]));
        j += 1;
      }
    }
    
    if (penalty_family == 2){
      int j = 1;
      for (i in B){
        // b_raw[j]  ~ normal(0, inv(sd_X[i+nXc]));
        b_raw[j]  ~ normal(0, inv(sd_X[i]));
        j += 1;
      }
    }
  }
