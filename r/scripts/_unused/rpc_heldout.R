`%*.%` = function(l, r) t(apply(r, 1, \(.) l * .))
`%.*.%` = function(l, r) t(sapply(seq_len(dim(l)[[1]]), 
                                  function(j) l[j,] * r[j,]))
load('data/cleaned/data_input.Rdata')

# pairwise_corr = matrix(NA, ncol=660, nrow = 2667)
mm = m3
p_Smear = mm$extract_heldout('p_Smear')$p_Smear 
p_Smear_Mgit = matrix(NA, ncol=660, nrow=nrow(p_Smear))
p_Smear_Xpert = matrix(NA, ncol=660, nrow=nrow(p_Smear))
p_Xpert_Mgit = matrix(NA, ncol=660, nrow=nrow(p_Smear))

for (i in seq_along(mm$model)){
  model = mm$model[[i]]
  holdout = mm$folds$holdout[[i]] |> as.logical()
  theta = rstan::extract(model, 'theta')$theta
  
  z_Smear = rstan::extract(model, 'z_Smear')$z_Smear
  z_Mgit = rstan::extract(model, 'z_Mgit')$z_Mgit
  z_Xpert = rstan::extract(model, 'z_Xpert')$z_Xpert
  
  p_Smear = rstan::extract(model, 'p_Smear')$p_Smear
  p_Mgit = rstan::extract(model, 'p_Mgit')$p_Mgit
  p_Xpert = rstan::extract(model, 'p_Xpert')$p_Xpert
  
  b_HIV = rstan::extract(model, 'b_HIV')$b_HIV
  b = rstan::extract(model, 'b')$b
  b_RE = rstan::extract(model, 'b_RE')$b_RE
  if (ncol(b_RE)==1) b_RE = cbind(b_RE, b_RE, b_RE)
  
  n_iters = dim(theta)[[1]]
  
  theta = theta[, holdout]
  p_Smear = p_Smear[, holdout]
  p_Mgit = p_Mgit[, holdout]
  p_Xpert = p_Xpert[, holdout]
  # browser()
  
  hiv = Xd[holdout, 1]
  V = Xc[holdout, 3:8]
  RE = matrix(rnorm(n_iters * sum(holdout)), nrow = n_iters)
  B = hiv %*.% b_HIV + b %*% t(V) + RE
  z_Smear_RE = z_Smear[,2] + B %.*.% b_RE[,1, drop=F]
  z_Mgit_RE = z_Mgit[,2] + B %.*.% b_RE[,2,drop=F]
  z_Xpert_RE = z_Mgit[,2] + B %.*.% b_RE[,3, drop=F]
  
  p_Smear_Mgit[, holdout] = theta * plogis(z_Smear_RE) * plogis(z_Mgit_RE) + (1 - theta) * plogis(z_Smear[,1]) * plogis(z_Mgit[,1])
  p_Smear_Xpert[, holdout] = theta * plogis(z_Smear_RE) * plogis(z_Xpert_RE) + (1 - theta) * plogis(z_Smear[,1]) * plogis(z_Xpert[,1])
  p_Xpert_Mgit[, holdout] = theta * plogis(z_Xpert_RE) * plogis(z_Mgit_RE) + (1 - theta) * plogis(z_Xpert[,1]) * plogis(z_Mgit[,1])
}

p_Smear = mm$extract_heldout('p_Smear')$p_Smear 
p_Mgit = mm$extract_heldout('p_Mgit')$p_Mgit
p_Xpert = mm$extract_heldout('p_Xpert')$p_Xpert

mu_i = cbind(
  rowMeans(p_Smear), 
  rowMeans(p_Mgit), 
  rowMeans(p_Xpert)
)

mu_ij = cbind(
  rowMeans(p_Smear_Mgit), 
  rowMeans(p_Smear_Xpert), 
  rowMeans(p_Xpert_Mgit)
)

pairwise_corr = matrix(NA, ncol=3, nrow=nrow(p_Smear))
{
  k = 1;
  for (i in 1:2){
    for (j in (i+1):3){
      pairwise_corr[,k] = (mu_ij[,k] - mu_i[,i]*mu_i[,j])/sqrt(mu_i[,i]*(1-mu_i[,i])*mu_i[,j]*(1-mu_i[,j]));
      k = k+1;
    }
  }
}

empirical_pairwise_corr = data_19EI[, c(pair_corr(csf_smear, csf_mgit), pair_corr(csf_smear, csf_xpert), pair_corr(csf_mgit, csf_xpert))] |> unlist()
residual_pairwise_corr = {
  a = matrix(rep(empirical_pairwise_corr, nrow(pairwise_corr)), byrow=TRUE, ncol = 3)
  pairwise_corr - a
}
