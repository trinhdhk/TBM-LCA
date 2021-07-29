  real a0; //intercept
  vector<lower=0>[nA_pos] a_pos; 
  vector<upper=0>[nA_neg] a_neg;
  vector[nA - nA_pos - nA_neg] a_;
