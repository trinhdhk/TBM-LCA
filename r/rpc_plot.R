models = list(m1=m1, m2=m2, m3=m3, m4=m4, m5=m5)
pairwise_corr = sapply(models, \(x) rstan::extract(x$outputs, pars='pairwise_corr')$pairwise_corr,
                       simplify = FALSE)
load('data/cleaned/data_input.Rdata')
pair_corr = function(x, y){
  mu_i = mean(x)
  mu_j = mean(y)
  mu_ij = mean(x*y)
  (mu_ij - mu_i * mu_j) / sqrt(mu_i * (1-mu_i) * mu_j * (1-mu_j))
}
empirical_pairwise_corr = data_19EI[, c(pair_corr(csf_smear, csf_mgit), pair_corr(csf_smear, csf_xpert), pair_corr(csf_mgit, csf_xpert))] |> unlist()
residual_pairwise_corr = sapply(pairwise_corr,
                                \(pc) {
                                  a = matrix(rep(empirical_pairwise_corr, nrow(pc)), byrow=TRUE, ncol = 3)
                                  pc - a
                                  }, simplify=FALSE)
rpc_summary = sapply(residual_pairwise_corr,
                    \(x) {
                      z = do.call(rbind,
                              apply(x, 2, 
                                    \(y) data.frame(mean = mean(y), ll = quantile(y, .025),  l = quantile(y, .25), h = quantile(y, .75), hh = quantile(y, .975)))
                      )
                      row.names(z) = NULL
                      z$Pair = c('Smear-MGIT', 'Smear-Xpert', 'MGIT-Xpert')
                      # z$Pair = factor(levels=c('Smear-MGIT', 'Smear-Xpert', 'MGIT-Xpert'))
                      z
                    }, simplify = FALSE)


library(ggplot2)
for (i in seq_len(5)){
  rpc_summary[[i]]$Model = i
}

rpc_summary = do.call(rbind, rpc_summary)
rpc_summary$Model = as.factor(rpc_summary$Model)
rpc_plot = ggplot(data=rpc_summary, aes(x=Pair, color=Model, group=Model)) + 
  geom_line(aes(y=mean), stat='summary', size = .7)+
  # stat_summary(aes(y=mean), geom='line', fun=sum) +
  scale_y_continuous(limits=c(-.15, .15)) +
  scale_x_discrete(expand=c(.07,.07)) +
  scale_color_brewer(palette='Set3') +
  theme_bw() + 
  theme(axis.title = element_blank())

saveRDS(rpc_plot, 'export/rpc_plot.RDS')
