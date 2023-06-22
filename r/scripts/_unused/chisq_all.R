`%*.%` = function(l, r) t(apply(r, 1, \(.) l * .))
`%.*.%` = function(l, r) t(sapply(seq_len(dim(l)[[1]]), 
                                  function(j) l[j,] * r[j,]))
load('data/cleaned/data_input.Rdata')
patterns = c('000', '001', '010', '011', '100', '101', '110', '111')
model = m3$outputs
theta = rstan::extract(model, 'theta')$theta

z_Smear = rstan::extract(model, 'z_Smear')$z_Smear
z_Mgit = rstan::extract(model, 'z_Mgit')$z_Mgit
z_Xpert = rstan::extract(model, 'z_Xpert')$z_Xpert

p_Smear = rstan::extract(model, 'p_Smear')$p_Smear
p_Mgit = rstan::extract(model, 'p_Mgit')$p_Mgit
p_Xpert = rstan::extract(model, 'p_Xpert')$p_Xpert

pat_prob = sapply(patterns, \(i) matrix(NA, nrow=nrow(p_Smear), ncol=660), simplify=FALSE)
b_HIV = rstan::extract(model, 'b_HIV')$b_HIV
b = rstan::extract(model, 'b')$b
b_RE = rstan::extract(model, 'b_RE')$b_RE
# b_RE2 = rstan::extract(model, 'b_RE2')$b_RE2

if (ncol(b_RE)==1) b_RE = cbind(b_RE, b_RE, b_RE)

n_iters = dim(theta)[[1]]

hiv = Xd[, 1]
# V = Xc[, 3:8]
V = rstan::extract(model, 'X')$X[,,13:18]
# RE = matrix(rnorm(n_iters * 660), nrow = n_iters)
RE = rstan::extract(model, 'RE')$RE
# B = hiv %*.% b_HIV + b %*% t(V) + RE
B = hiv %*.% b_HIV + t(sapply(1:10668, function(i) b[i, ,drop=F] %*% t(V[i,,]))) + RE
z_Smear_RE = z_Smear[,2] + B %.*.% b_RE[,1, drop=F] #+ B^2 %.*.% b_RE[,1, drop=F]
z_Mgit_RE = z_Mgit[,2] + B %.*.% b_RE[,2,drop=F] #+ B^2 %.*.% b_RE[,1, drop=F]
z_Xpert_RE = z_Mgit[,2] + B %.*.% b_RE[,3, drop=F] #+ B^2 %.*.% b_RE[,1, drop=F]
# z_Xpert_RE = rep(-1, length(z_Xpert_RE))

for (pattern in patterns){
  pat = as.numeric(strsplit(pattern, '')[[1]])
  p11 = if (pat[[1]]) plogis(z_Smear_RE) else (1 - plogis(z_Smear_RE) )
  p12 = if (pat[[2]]) plogis(z_Mgit_RE) else (1 - plogis(z_Mgit_RE) )
  p13 = if (pat[[3]]) plogis(z_Xpert_RE) else (1 - plogis(z_Xpert_RE) )
  p01 = if (pat[[1]]) plogis(z_Smear[,1]) else (1 - plogis(z_Smear[,1]) )
  p02 = if (pat[[2]]) plogis(z_Mgit[,1]) else (1 - plogis(z_Mgit[,1]) )
  p03 = if (pat[[3]]) plogis(z_Xpert[,1]) else (1 - plogis(z_Xpert[,1]) )
  
  pat_prob[[pattern]] = theta * p11 * p12 * p13 + (1-theta) * p01 * p02 * p03
}
empirical_pat_prob = data_19EI[, paste0(csf_smear*1, csf_mgit*1, csf_xpert*1)] |> table()
empirical_pat_prob2 = c(empirical_pat_prob[1], 0, empirical_pat_prob[2:7])
sapply(pat_prob, function(x) rowMeans(x)) -> mean_pat_prob
exp_p = colMeans(mean_pat_prob)
exp_p = c(exp_p[1]+exp_p[2], exp_p[3:8])
Mean_pat_prob2 = sapply(pat_prob, rowMeans)
Chisq=apply(Mean_pat_prob2[1:1000,],1, function(m) chisq.test(empirical_pat_prob2, p=m, simulate.p.value = T)$statistic)
