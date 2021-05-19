//Customisable penalty
  /*
  Different penalty values will have different prior families:
    0 = student_t(4);
    1 = double_exponential();
    2 = normal();
  
  penalty: if == 0 will perform auto adaptation for penalty terms,
  currently, RE model and prevalence will have different penalty if adpatation
  is turned on.
  if != 0 will create manual penalty.
  */
  int<lower=0, upper=2> penalty_family;
  real<lower=0> penalty_term;
