# Script to create necessity for posterior report for model 3
# Author: Trinh Dong
##################################################################
library(bayesplot)
library(ggplot2)
library(tidyverse)
load('data/cleaned/data_input.Rdata')
scale = purrr::transpose(scale_Xc)
model = m3

#### p_test for C=1
p_test_RE =rstan::extract(model$outputs, c('z_Xpert_RE', 'z_Smear_RE', 'z_Mgit_RE'))
hiv = rstan::extract(model$outputs, 'X')$X[,,1] |> apply(2, as.logical)
p_test_re_all =  lapply(p_test_RE, \(x) apply(x, 1, \(.x) mean(plogis(.x))))
p_test_re_hiv = lapply(p_test_RE, \(x) {
  x[!hiv] = NA
  apply(x, 1, \(.x) mean(plogis(.x), na.rm=T))
  })
p_test_re_nonhiv = lapply(p_test_RE, \(x) {
  x[hiv] = NA
  apply(x, 1, \(.x) mean(plogis(.x), na.rm=T))
})

summp = function(p) {
  summ = sapply(p, function(pp) c(mean = mean(pp),
                           quantile(pp, c(.5, 0.025, .975)))) |>
    t()
  rownames(summ) <- gsub('_RE', '[2]', rownames(summ))
  cbind(as.data.frame(summ), parameter=rownames(summ))
}
p_re_summary <- summp(p_test_re_all)
p_re_summary_hiv <- summp(p_test_re_hiv)
p_re_summary_nohiv <- summp(p_test_re_nonhiv)

summp = function(p) {
  summ = sapply(p, function(pp) c(mean = mean(pp),
                                  quantile(pp, c(.5, 0.025, .975)))) |>
    t()
  rownames(summ) <- gsub('_RE', '[2]', rownames(summ))
  cbind(as.data.frame(summ), parameter=rownames(summ))
}
p_re_summary <- summp(p_test_re_all)
p_hiv_pos_summary <- summp(p_test_re_hiv)
p_hiv_neg_summary <- summp(p_test_re_nonhiv)


## p_test for spc

#### p_test for C=1
p_test0 =rstan::extract(model$outputs, c('z_Xpert', 'z_Smear', 'z_Mgit')) |> lapply("[",,1) |> lapply(plogis)
p_test0_summary = sapply(p_test0, function(pp) c(mean = mean(pp),
                                                 quantile(pp, c(.5, 0.025, .975)))) |> t() |> as.data.frame()
rownames(p_test0_summary) = paste0(rownames(p_test0_summary), '[1]')
p_test0_summary$parameter = rownames(p_test0_summary)
p_summary = rbind(p_test0_summary, p_re_summary)

a = rstan::extract(model$outputs, pars=c('a0', 'a'))
# a$a[,1:ncol(Xd)] = a$a[,1:ncol(Xd)] * 2
a$a[,18] = -a$a[,18]
a = cbind(a$a0, a$a) 
colnames(a) = c('a0', paste0('a[',1:(ncol(a)-1),']'))
a_orig = a
for (i in 12:(ncol(a_orig)-1)) a_orig[,i] = a_orig[,i]/scale$`scaled:scale`[[i-11]]
# intercept_translation = unlist(scale$`scaled:center`) * a_orig[, 11:(ncol(a_orig)-2)]
# intercept_translation[, 'a[18]'] = -intercept_translation[, 'a[18]'] + 15 * a_orig[, 'a[18]']
# intercept_translation[, 'a[18]'] = intercept_translation[, 'a[18]']  (15 * (-a_orig[,'a[18]']))
# a_orig[,1] = a_orig[,1] - apply(intercept_translation,1,sum)
scale_wbc2 = sd(Xc[,7]^2)
a_orig[,ncol(a_orig)]=a_orig[,ncol(a_orig)]/scale_wbc2
labs = c(
  'Intercept', 
  'HIV +',
  'TB-suggestive symptoms',
  'Focal neurological deficit',
  'Cranial nerve palsy',
  'Past noticed TB contact',
  'Glasgow coma score',
  'Pulmonary TB/X-ray',
  'Miliary TB/X-ray',
  '*log<sub>2</sub>* (Symptom duration, days)',
  '*log<sub>2</sub>* (Paired blood glucose)',
  '*log<sub>2</sub>* (CSF glucose)',
  '*log<sub>2</sub>* (CSF protein)',
  '*log<sub>2</sub>* (CSF lactate)',
  '*log<sub>10</sub>* (CSF lymphocyte)',
  # '*log<sub>10</sub>* (CSF WBC count)',
  # '*log<sub>10</sub>* (CSF WBC count)<sup>2</sup>',
  'CSF eosinophil > 0',
  '*log<sub>10</sub>* (CSF eosinophil)',
  '*log<sub>10</sub>* (CSF RBC)',
  'Evidence of Cryptococcus',
  'Positive CSF Gram stain'
) 

f2 = function(x) formatC(x, digits=2, format='f')

add_annotate = function(tab, text = '{f2(mean)}; {f2(`50%`)} [{f2(`2.5%`)}; {f2(`97.5%`)}]'){
  tab |>
    as.data.frame() %>%
    dplyr::mutate(
      parameter = rownames(.),
      annotate = glue::glue(text)
    )
}


z_summary = rstan::summary(model$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert'))$summary |> add_annotate()
a_summary = rstan::summary(model$outputs, pars=c('a0', 'a'))$summary |> add_annotate()
b_summary = rstan::summary(model$outputs, pars=c('b_RE','b_HIV', 'b'))$summary |> add_annotate()


save(a_summary, b_summary, z_summary, p_summary, p_hiv_neg_summary, p_hiv_pos_summary,  file='export/m3_summary.Rdata')

cutpoints = c(-20,-16,-8,-4, -2,-1, log10(.5), log10(.9) ,0, log10(10/9),log10(2),1,2,4)
lab.cutpoints = c(paste0('1e', c(-20,-16,-8,-4)), .01, .1,.5, .9,1, 1.1, 2, 10 ,100, '1e4')
a_plot <-
  (td.misc::mcmc_intervals_multi(
    fits = list(a_orig),
    pars = c('a0',
             paste0('a[',c(1:5,18,6:7,11:13,15,16,14,10,19,20,8,9),']')),
    multi_point_est = FALSE,
    point_est = c('median'),
    transformations = function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))),
    point_size = 2,
    prob_outer = .95) +
     geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
     scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)], 
                          breaks = c('a', 'a_orig'),
                          labels=c('Matched', 'Original')) + 
     scale_x_continuous(name='TBM odds ratio',
                        # breaks=log(10^cutpoints),
                        minor_breaks = NULL,
                        breaks = (function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))))(log(10^cutpoints)),
                        labels=lab.cutpoints) +
     ylab('')+
     theme_bw() +
     theme(axis.text.y = ggtext::element_markdown()) +
  guides(color=guide_legend(title="Scale"))) %>%
  td.misc::change_ylabs(
    labs = labs,
    top_down = T
  ) + theme(legend.position = 'none') #+ 

b = rstan::extract(model$outputs, pars=c('b_HIV', 'b'))
b$b[,6]=-b$b[,6]
b = cbind(b$b_HIV, b$b)
colnames(b) = c('b_HIV', paste0('b[',1:(ncol(b)-1), ']'))
b_orig = b
re_b = as.numeric(model$.META$re_b)
for (i in 2:(ncol(b))) b_orig[,i] = b_orig[,i]/(scale$`scaled:scale`[[re_b[i-1]]])
b_plot = td.misc::mcmc_intervals_multi(list(b_orig),
                        pars = c('b_HIV',  paste0('b[',c(6,1,4,5,2,3), ']')),
                        multi_point_est = FALSE,
                        point_est = c('median'),
                        point_size = 2,
                        prob_outer = .95) |>
  td.misc::change_ylabs(
    'HIV +' ,
    'GCS',
    '*log<sub>2</sub>* (CSF glucose)',
    '*log<sub>2</sub>* (CSF protein)',
    '*log<sub>2</sub>* (CSF lactate)',
    '*log<sub>10</sub>* (CSF lymphocyte)',
    '*log<sub>10</sub>* (CSF white cell)',
    top_down = T
  ) +
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)], 
                       breaks=c('b', 'b_orig'),
                        labels=c('Matched', 'Original')
                       ) + 
  guides(color=guide_legend(title="Scale")) +
  ylab('') +
  scale_x_continuous(name='Standardised bacillary burden', breaks=c(-3,-2,-1,0,1,2,3))+
  theme_bw() + theme(axis.text.y = ggtext::element_markdown(), legend.position = 'none')

library(patchwork)
ab_plot = a_plot / b_plot + plot_layout(guides='collect', ncol=1,
                                        nrow = 2, heights=c(3.5,1))

z = rstan::extract(model$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert'))
z = do.call(cbind, z)
colnames(z) = paste0('z_',rep(c('Smear', 'Mgit', 'Xpert'),each=2), rep(c('[1]', '[2]'), 3))
z_plot = td.misc::mcmc_intervals_multi(list(z=z), 
                        regex_pars = '^z',
                        transformations = 'plogis',
                        multi_point_est = FALSE,
                        point_est = c('median'),
                        point_size = 2,
                        prob_outer = .95)

z_plot_logit = mcmc_intervals(z, 
                             regex_pars = '^z',
                             multi_point_est = TRUE,
                             point_est = 'median',
                             point_size = 2,
                             prob_outer = .95)

wbc = rbind(Xc[,7], Xc[,7]^2)
wbca = a[,c('a[17]', 'a[21]')]
wbcpred = wbca %*% wbc
wbcpredsum = data.frame(
  mean = apply(wbcpred,2,mean),
  l = apply(wbcpred,2,quantile,.25),
  ll = apply(wbcpred,2,quantile,.025),
  h = apply(wbcpred,2,quantile,.75),
  hh = apply(wbcpred,2,quantile,.975),
  wbc = wbc[1,]
) 

trans=function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x)))
argmax_wbc = (wbcpredsum$wbc[which.max(wbcpredsum$mean)]+scale$`scaled:center`$csfwbc) * scale$`scaled:scale`$csfwbc
wbc_plot = ggplot(wbcpredsum, aes(x=wbc)) +
  geom_ribbon(aes(ymax=h,ymin=l), alpha=.4, fill = RColorBrewer::brewer.pal(name='Set1', n=3)[2])+
  geom_ribbon(aes(ymax=hh,ymin=ll), alpha=.2, fill = RColorBrewer::brewer.pal(name='Set1', n=3)[2]) + 
  scale_x_continuous(name = 'CSF white cell count (cells/mm<sup>3</sup>)',
                     breaks = c((log10(c(1, 11, 101, 1001, 10001)) - scale$`scaled:center`$csfwbc)/scale$`scaled:scale`$csfwbc, wbcpredsum$wbc[which.max(wbcpredsum$mean)]),
                     labels = c(0, 10, 100, 1000, 10000, round(10^argmax_wbc,0)), expand=c(0,0)) +
  geom_hline(aes(yintercept=0), color=grey(.5), alpha=.5) +
  geom_vline(aes(xintercept=wbc[which.max(mean)]), color=grey(.5), alpha=.5)+
  geom_line(aes(y=mean), 
            color=RColorBrewer::brewer.pal(name='Set1', n=3)[2], size=.8) + 
  scale_y_continuous(
    name = glue::glue('TBM odds ratio relative to reference value = {10^(scale$`scaled:center`$csfwbc*scale$`scaled:scale`$csfwbc) |> formatC(format = "f", digits=0)} cells/mm<sup>3</sup>') |> unclass(),
    breaks = log(10^c(-24,-20, -16,-12, -8, -6, -4,-2, -1, 0)),
    labels = c('1e-24','1e-20','1e-16','1e-12','1e-8','1e-6','1e-4', 0.01, .1, 1)
  ) +
  coord_flip(
    xlim = (log10(c(1, 10001)) - scale$`scaled:center`$csfwbc)/scale$`scaled:scale`$csfwbc,
    ylim = c(log(10^-6), 0))+
  theme_bw() +
  theme(axis.title.x = ggtext::element_markdown(), axis.title.y = ggtext::element_markdown())
save(z_plot, z_plot_logit, a_plot, b_plot, ab_plot, wbc_plot, file = 'export/m3_plot.Rdata')
saveRDS(a_orig, 'export/a_orig.RDS')
