`%*.%` = function(l, r) t(apply(r, 1, \(.) l * .))
`%.*.%` = function(l, r) t(sapply(seq_len(dim(l)[[1]]), 
                                  function(j) l[j,] * r[j,]))
load('data/cleaned/data_input.Rdata')

# pairwise_corr = matrix(NA, ncol=660, nrow = 2667)
mm = m5
model_id = 5
p_Smear = mm$extract_heldout('p_Smear')$p_Smear 
patterns = c('000', '001', '010', '011', '100', '101', '110', '111')
pat_prob = sapply(patterns, \(i) matrix(NA, nrow=nrow(p_Smear), ncol=660), simplify=FALSE)

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
  
  if (model_id > 1){
    b_HIV = rstan::extract(model, 'b_HIV')$b_HIV
    b = rstan::extract(model, 'b')$b
    b_RE = rstan::extract(model, 'b_RE')$b_RE
    
    if (ncol(b_RE)==1) b_RE = cbind(b_RE, b_RE, b_RE)
  }
 
  
  n_iters = dim(theta)[[1]]
  
  theta = theta[, holdout]
  p_Smear = p_Smear[, holdout]
  p_Mgit = p_Mgit[, holdout]
  p_Xpert = p_Xpert[, holdout]
  # browser()
  
  hiv = Xd[holdout, 1]
  V = Xc[holdout, 3:8]
  RE = matrix(rnorm(n_iters * sum(holdout)), nrow = n_iters)
  # browser()
  if (model_id %in% 2:4){
    B = hiv %*.% b_HIV + b %*% t(V) + RE
    z_Smear_RE = z_Smear[,2] + B %.*.% b_RE[,1, drop=F]
    z_Mgit_RE = z_Mgit[,2] + B %.*.% b_RE[,2,drop=F]
    z_Xpert_RE = z_Mgit[,2] + B %.*.% b_RE[,3, drop=F]
  } else if (model_id == 5){
    z_Smear_RE = z_Smear[,2] + hiv %*.% b_HIV[,1, drop=F] + b[,,1] %*% t(V) + RE %.*.% b_RE[,1, drop=F]
    z_Mgit_RE = z_Mgit[,2] + hiv %*.% b_HIV[,2, drop=F] + b[,,2] %*% t(V) + RE %.*.% b_RE[,2,drop=F]
    z_Xpert_RE = z_Mgit[,2] + hiv %*.% b_HIV[,3, drop=F] + b[,,3] %*% t(V) + RE %.*.% b_RE[,3, drop=F]
  } else {
    z_Smear_RE = z_Smear[,2]
    z_Mgit_RE = z_Mgit[,2]
    z_Xpert_RE = z_Xpert[,2]
  }
  
  
  for (pattern in patterns){
    pat = as.numeric(strsplit(pattern, '')[[1]])
    p11 = if (pat[[1]]==1) plogis(z_Smear_RE) else (1 - plogis(z_Smear_RE)) 
    p12 = if (pat[[2]]==1) plogis(z_Mgit_RE) else (1 - plogis(z_Mgit_RE) )
    p13 = if (pat[[3]]==1) plogis(z_Xpert_RE) else (1 - plogis(z_Xpert_RE)) 
    p01 = if (pat[[1]]==1) plogis(z_Smear[,1]) else (1 - plogis(z_Smear[,1]) )
    p02 = if (pat[[2]]==1) plogis(z_Mgit[,1]) else (1 - plogis(z_Mgit[,1])) 
    p03 = if (pat[[3]]==1) plogis(z_Xpert[,1]) else (1 - plogis(z_Xpert[,1]) )
    
    pat_prob[[pattern]][, holdout] = theta * p11 * p12 * p13 + (1-theta) * p01 * p02 * p03
  }
  
}
empirical_pat_prob = data_19EI[, paste0(as.numeric(csf_smear), as.numeric(csf_mgit), as.numeric(csf_xpert))] |> table() |> as.vector()
sapply(pat_prob, rowMeans) -> mean_pat_prob
expected_pat_prob = apply(mean_pat_prob,2,mean)*660
# chisq.test(empirical_pat_prob, expected_pat_prob, simulate.p.value = T)
