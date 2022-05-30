{
  vector[3] mu_i;
  vector[3] mu_ij;
  vector[N_all] p_Smear_Mgit;
  vector[N_all] p_Smear_Xpert;
  vector[N_all] p_Mgit_Xpert;
  {
    matrix[N_all,3] y;
    for (n in 1:N_all){
      y[n,1] = bernoulli_rng(p_Smear[n]);
      y[n,2] = bernoulli_rng(p_Mgit[n]);
      y[n,3] = bernoulli_rng(p_Xpert[n]);
    }
    for (i in 1:3) mu_i[i] = mean(y[:,i]);
  }
  
  p_Smear_Mgit = (1 - theta)*inv_logit(z_Smear[1])*inv_logit(z_Mgit[1]) +
    theta .* inv_logit(z_Smear_RE) .* inv_logit(z_Mgit_RE); 
  p_Smear_Xpert= (1 - theta)*inv_logit(z_Smear[1])*inv_logit(z_Xpert[1]) +
    theta .* inv_logit(z_Smear_RE) .* inv_logit(z_Xpert_RE); 
  p_Mgit_Xpert = (1 - theta)*inv_logit(z_Mgit[1])*inv_logit(z_Xpert[1]) +
    theta .* inv_logit(z_Mgit_RE) .* inv_logit(z_Xpert_RE); 
  {
    matrix[N_all,3] y;
    for (n in 1:N_all){
      y[n,1] = bernoulli_rng(p_Smear_Mgit[n]);
      y[n,2] = bernoulli_rng(p_Smear_Xpert[n]);
      y[n,3] = bernoulli_rng(p_Mgit_Xpert[n]);
    }
    for (i in 1:3) mu_ij[i] = mean(y[:,i]);
  }
  
  {
    int k = 1;
    for (i in 1:2){
      for (j in (i+1):3){
        pairwise_corr[k] = (mu_ij[k] - mu_i[i]*mu_i[j])/sqrt(mu_i[i]*(1-mu_i[i])*mu_i[j]*(1-mu_i[j]));
        k += 1;
      }
    }
  }
}
