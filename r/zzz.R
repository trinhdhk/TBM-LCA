p = rstan::extract(m3s2_n00_b45678_q8_r1_k1$outputs, pars=c('p_Smear', 'p_Mgit', 'p_Xpert'))
data_19EI[,hiv_stat2:=ifelse(is.na(hiv_stat), 0, hiv_stat)]
# hiv = data_19EI[,!hiv_stat2]
hiv=TRUE

Y = sapply(p, function(pp) apply(pp, 2, function(r) sapply(r, rbinom, n=500, size=1)), simplify = F, USE.NAMES = T)
# Y=cbind(Y$p_Smear,Y$p_Mgit,Y$p_Xpert)
pair_corr <- function(p, corr_list){
  pair <- pair.matrix[p,]
  i <- pair[[1]]; j <- pair[[2]]
  mu_i <- if (is.matrix(corr_list$mu_i)) corr_list$mu_i[,i] else corr_list$mu_i[i]
  mu_j <- if (is.matrix(corr_list$mu_i)) corr_list$mu_i[,j] else corr_list$mu_i[j]
  
  mu_ij <- if (is.matrix(corr_list$mu_ij)) corr_list$mu_ij[,p] else corr_list$mu_ij[p]
  (mu_ij  - mu_i*mu_j)/sqrt(mu_i*(1-mu_i)*mu_j*(1-mu_j))
}
com.matrix <- function(M)
  do.call(expand.grid,split(M,rep(1:nrow(M),ncol(M))))
pair.matrix <- com.matrix(rbind(seq_len(3), seq_len(3)))
pair.matrix <- subset(pair.matrix, (pair.matrix[,1]<pair.matrix[,2]))

Y2 = sapply(Y, \(y) y[,hiv], simplify=F)
corr_pred = list()
corr_pred$mu_i = sapply(Y2, apply, 1, mean)
corr_pred$mu_ij = apply(pair.matrix, 1,
                        function(pair) {
                          M = t(sapply(1:6400, \(j) Y2[[pair[[1]]]][j,] * Y2[[pair[[2]]]][j,]))
                          # print(dim(M))
                          # print(apply(M, 1, mean) |> length())
                          apply(M, 1, mean) 
                        })
corr_pred$corr_ij <- sapply(seq_len(nrow(pair.matrix)), pair_corr, corr_list = corr_pred)



obs =data_19EI[hiv==T,.(csf_smear, csf_mgit, csf_xpert)]|>as.matrix()
corr_obs = list()
corr_obs$mu_i <- apply(obs, 2, mean, na.rm=T)
corr_obs$mu_ij <- apply(pair.matrix, 1,
                        function(pair)
                          mean(obs[,pair[1]]*obs[,pair[2]], na.rm=T))
corr_obs$corr_ij <-
  sapply(seq_len(nrow(pair.matrix)), pair_corr, corr_list = corr_obs)

library(patchwork)

plt  = vector('list',3)

fn = function(){
  i = i
  ggplot() + geom_density(aes(x=corr_pred$corr_ij)) + geom_vline(aes(xintercept=corr_obs$corr_ij[i]))
} 

for (i in 1:3) plt[[i]] = fn()
patchwork::wrap_plots(plt, guides = 'collect')
