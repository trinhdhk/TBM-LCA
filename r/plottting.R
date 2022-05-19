library(bayesplot)
library(ggplot2)
load('data/cleaned/data_input.Rdata')
scale = purrr::transpose(scale_Xc)
model = m3t3
a = rstan::extract(model$outputs, pars=c('a0', 'a'))
# a$a[,1:ncol(Xd)] = a$a[,1:ncol(Xd)] * 2
a$a[,18] = -a$a[,18]
a = cbind(a$a0, a$a) 
colnames(a) = c('a0', paste0('a[',1:(ncol(a)-1),']'))
a_orig = a
for (i in 11:(ncol(a_orig)-2)) a_orig[,i] = a_orig[,i]/scale$`scaled:scale`[[i-10]]
# intercept_translation = unlist(scale$`scaled:center`) * a_orig[, 11:(ncol(a_orig)-2)]
# intercept_translation[, 'a[18]'] = -intercept_translation[, 'a[18]'] + 15 * a_orig[, 'a[18]']
# intercept_translation[, 'a[18]'] = intercept_translation[, 'a[18]']  (15 * (-a_orig[,'a[18]']))
# a_orig[,1] = a_orig[,1] - apply(intercept_translation,1,sum)
scale_wbc2 = sd(Xc[,7])
a_orig[,ncol(a_orig)]=a_orig[,ncol(a_orig)]/scale_wbc2
labs = c(
  'Intercept (centred)', 
  'HIV +',
  'TB-suggested symptoms',
  'Focal neurological deficit +',
  'Cranial nerve palsy +',
  'Past noticed TB contact',
  'Glasgow coma score',
  'Pulmonary TB/X-ray',
  'Miliary TB/X-Ray',
  '*log<sub>2</sub>* (Days from onset)',
  '*log<sub>2</sub>* (Paired blood glucose)',
  '*log<sub>2</sub>* (CSF Glucose)',
  '*log<sub>2</sub>* (CSF Protein)',
  '*log<sub>2</sub>* (CSF Lactate)',
  '*log<sub>10</sub>* (CSF Lymphocyte)',
  # '*log<sub>10</sub>* (CSF WBC count)',
  # '*log<sub>10</sub>* (CSF WBC count)<sup>2</sup>',
  'CSF Eosinophil > 0',
  '*log<sub>10</sub>* (CSF Eosinophil)',
  '*log<sub>10</sub>* (CSF RBC)',
  'Cryptococcus + b',
  'CSF Gram stain +'
) 

add_annotate = function(tab, text = '{round(mean,2)} [{round(`2.5%`,2)}; {round(`97.5%`,2)}]'){
  tab |>
    as.data.frame() %>%
    dplyr::mutate(
      parameter = rownames(.),
      annotate = glue::glue(text)
    )
}

p_hiv_neg_summary = {
  z = rstan::extract(model$outputs, pars=c("z_Smear", "z_Mgit", "z_Xpert"))
  z = sapply(z, plogis, simplify = FALSE)
  z = sapply(z, function(x) {apply(x, 2, function(.x) {
    m = c(mean(.x), median(.x), quantile(.x, c(.025, .975)))
    names(m) <- c('mean', '50%', '2.5%', '97.5%')
    m
    })}, simplify = FALSE)
  z = sapply(z, t, simplify = FALSE)
  z = sapply(1:3, function(i) {
    m = z[[i]]
    rownames(m) = paste0(names(z)[i], '[',1:2,']')
    m
  }, simplify = FALSE)
  do.call(rbind, z)
} |> add_annotate()

p_hiv_pos_summary = {
  z = rstan::extract(model$outputs, pars=c("z_Smear", "z_Mgit", "z_Xpert"))
  b = rstan::extract(model$outputs, pars=c("b_HIV", "b_RE"))
  k = 1
  for (i in c("z_Smear", "z_Mgit", "z_Xpert")){
    z[[i]][,2] = z[[i]][,2] + b$b_HIV * b$b_RE[,k]
    k = k+1
  }
  
  z = sapply(z, plogis, simplify = FALSE)
  z = sapply(z, function(x) {apply(x, 2, function(.x) {
    m = c(mean(.x), median(.x), quantile(.x, c(.025, .975)))
    names(m) <- c('mean', '50%', '2.5%', '97.5%')
    m
  })}, simplify = FALSE)
  z = sapply(z, t, simplify = FALSE)
  z = sapply(1:3, function(i) {
    m = z[[i]]
    rownames(m) = paste0(names(z)[i], '[',1:2,']')
    m
  }, simplify = FALSE)
  do.call(rbind, z)
} |> add_annotate()

p_summary = {
  z = rstan::extract(model$outputs, pars=c("z_Smear", "z_Mgit", "z_Xpert"))
  z.neg = sapply(z, plogis, simplify = FALSE)
  
  b = rstan::extract(model$outputs, pars=c("b_HIV", "b_RE"))
  k = 1
  for (i in c("z_Smear", "z_Mgit", "z_Xpert")){
    z[[i]][,2] = z[[i]][,2] + b$b_HIV * b$b_RE[,k]
    k = k+1
  }
  z.pos = sapply(z, plogis, simplify = FALSE)
  
  hiv = rstan::extract(model$outputs, pars='X')$X[,,1]
  p_hiv = apply(hiv,1,mean)
  z2 = sapply(1:3, function(i) z.pos[[i]][,2] * (p_hiv) + z.neg[[i]][,2] * (1-p_hiv) , simplify = FALSE)
  z3 = lapply(1:3, function(i) cbind(z.neg[[i]][,1], z2[[i]]))
  names(z3) = names(z)
  z = z3
  
  
  z = sapply(z, function(x) {apply(x, 2, function(.x) {
    m = c(mean(.x), median(.x), quantile(.x, c(.025, .975)))
    names(m) <- c('mean', '50%', '2.5%', '97.5%')
    m
  })}, simplify = FALSE)
  z = sapply(z, t, simplify = FALSE)
  z = sapply(1:3, function(i) {
    m = z[[i]]
    rownames(m) = paste0(names(z)[i], '[',1:2,']')
    m
  }, simplify = FALSE)
  do.call(rbind, z)
}

z_summary = rstan::summary(model$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert'))$summary |> add_annotate()
a_summary = rstan::summary(model$outputs, pars=c('a0', 'a'))$summary |> add_annotate()
b_summary = rstan::summary(model$outputs, pars=c('b_RE','b_HIV', 'b'))$summary |> add_annotate()


save(a_summary, b_summary, z_summary, p_summary, p_hiv_neg_summary, p_hiv_pos_summary,  file='export/m3_summary.Rdata')

cutpoints = c(-20,-16,-8,-4, -2,-1, log10(.5), log10(.9) ,0, log10(10/9),log10(2),1,2,4)
# cutpoints = c(-6,-4,-2,-1 ,0 ,1,2,4, 6)
# lab.cutpoints = c(paste0('1e', c(-6,-4)), .01, .1, 1, 10 ,100, '1e4', '1e6')
lab.cutpoints = c(paste0('1e', c(-20,-16,-8,-4)), .01, .1,.5, .9,1, 1.1, 2, 10 ,100, '1e4')
a_plot <-
  (td.misc::mcmc_intervals_multi(
    fits = list(a_orig),
    pars = c('a0',
             paste0('a[',c(1:5,18,6:7,11:13,15,16,14,10,19,20,8,9),']')),
    multi_point_est = TRUE,
    point_est = c('mean', 'median'),
    transformations = function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))),
    # point_size = 4,
    prob_outer = .95) +
     geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
     # geom_text(aes(x = 40, y=parameter, label=annotate), 
     #           data = a_summary, size = 3, hjust='right',
     #           inherit.aes = FALSE) +
     # scale_x_continuous(limits=c(-40, 40), breaks=seq(-40, 40, 10)) + 
     scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)], 
                          breaks = c('a', 'a_orig'),
                          labels=c('Matched', 'Original')) + 
     scale_x_continuous(name='TBM odd ratio', 
                        # breaks=log(10^cutpoints),
                        breaks = (function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))))(log(10^cutpoints)),
                        labels=lab.cutpoints) +
     ylab('')+
     theme_bw() +
     theme(axis.text.y = ggtext::element_markdown(face='bold')) + 
  guides(color=guide_legend(title="Scale"))) %>%
  td.misc::change_ylabs(
    labs = labs,
    top_down = T
  ) + theme(legend.position = 'none') #+ 
  # ggforce::facet_zoom(xlim=c(-8,8))
  # (\(x) x )
a_plot$data$m = bayesplot::mcmc_intervals_data(a_orig,
                                             pars = c('a0',
                                                      paste0('a[',c(1:5,18,6:7,11:13,15,16,14,10,19,20,8,9),']')),
                                             transformations = exp, prob_outer=.95, point_est='mean') %>%
  mutate(across(c(ll:hh), ~ sign(qlogis(.x)) * sqrt(abs(qlogis(.x))))) %>% .$m

b = rstan::extract(model$outputs, pars=c('b_HIV', 'b'))
# b$b = b$b * 2
b$b[,6]=-b$b[,6]
b = cbind(b$b_HIV, b$b)
colnames(b) = c('b_HIV', paste0('b[',1:(ncol(b)-1), ']'))
b_orig = b
re_b = as.numeric(model$.META$re_b)
for (i in 2:(ncol(b))) b_orig[,i] = b_orig[,i]/(scale$`scaled:scale`[[re_b[i-1]]])
# intercept_translation_b = unlist(scale$`scaled:center`[re_b]) * b_orig[,-1]
# b_orig[,1] = b_orig[,1] - apply(intercept_translation_b,1,sum)
b_plot = td.misc::mcmc_intervals_multi(list(b_orig),
                        pars = c('b_HIV',  paste0('b[',c(1,6,4,5,2,3), ']')),
                        multi_point_est = TRUE,
                        point_est = c('mean', 'median'),
                        # point_size = 2.5,
                        prob_outer = .95) |>
  td.misc::change_ylabs(
    'HIV +',
    # 'Local neurological deficit +',
    'GCS',
    '*log<sub>2</sub>* (CSF Glucose)',
    '*log<sub>2</sub>* (CSF Protein)',
    '*log<sub>2</sub>* (CSF Lactate)',
    '*log<sub>10</sub>* (CSF Lymphocyte)',
    '*log<sub>10</sub>* (CSF WBC)',
    top_down = T
  ) +
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  # scale_x_continuous(breaks = seq(-14,14,2)) + 
  scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)], 
                       # breaks = 2:1,
                       breaks=c('b', 'b_orig'),
                        labels=c('Matched', 'Original')
                       ) + 
  guides(color=guide_legend(title="Scale")) +
  ylab('') +
  scale_x_continuous(name='Standardised bacillary burden', breaks=c(-3,-2,-1,0,1,2,3))+
  theme_bw() + theme(axis.text.y = ggtext::element_markdown(family=NULL, face='bold'), legend.position = 'none')

library(patchwork)
ab_plot = a_plot / b_plot + plot_layout(guides='collect', ncol=1,
                                        nrow = 2, heights=c(3.5,1))

z = rstan::extract(model$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert'))
z = do.call(cbind, z)
colnames(z) = paste0('z_',rep(c('Smear', 'Mgit', 'Xpert'),each=2), rep(c('[1]', '[2]'), 3))
z_plot = td.misc::mcmc_intervals_multi(list(z=z), 
                        regex_pars = '^z',
                        transformations = 'plogis',
                        multi_point_est = TRUE,
                        point_est = c('mean', 'median'),
                        point_size = 2.5,
                        prob_outer = .95)

z_plot_logit = mcmc_intervals(z, 
                             regex_pars = '^z',
                             multi_point_est = TRUE,
                             point_est = 'mean',
                             point_size = 2.5,
                             prob_outer = .95)

wbc = rbind(Xc[,7], Xc[,7]^2)
wbca = a[,c('a[17]', 'a[21]')]
wbcpred = wbca %*% wbc
# wbcpred = trans(wbcpred)
wbcpredsum = data.frame(
  mean = apply(wbcpred,2,mean),
  l = apply(wbcpred,2,quantile,.25),
  ll = apply(wbcpred,2,quantile,.025),
  h = apply(wbcpred,2,quantile,.75),
  hh = apply(wbcpred,2,quantile,.975),
  # wbc_org = data_19EI$csf_wbc_corrected |> tidyr::replace_na(0),
  wbc = wbc[1,]
) 

trans=function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x)))
  
wbc_plot = ggplot(wbcpredsum, aes(x=wbc)) +
  geom_ribbon(aes(ymax=h,ymin=l), alpha=.4, fill = RColorBrewer::brewer.pal(name='Set1', n=3)[2])+
  geom_ribbon(aes(ymax=hh,ymin=ll), alpha=.2, fill = RColorBrewer::brewer.pal(name='Set1', n=3)[2]) + 
  scale_x_continuous(name = 'CSF WBC',
                     breaks = (log10(c(1, 11, 101, 1001, 10001)) - scale$`scaled:center`$csfwbc)/scale$`scaled:scale`$csfwbc,
                     labels = c(0, 10, 100, 1000, 10000), expand=c(0,0)) +
  geom_hline(aes(yintercept=0), color=grey(.5), alpha=.5) +
  geom_line(aes(y=mean), 
            color=RColorBrewer::brewer.pal(name='Set1', n=3)[2], size=.8) + 
  scale_y_continuous(
    name = 'TBM odd ratio',
    breaks = log(10^c(-24,-20, -16,-12, -8,-4,-2, -1, 0)),
    labels = c('1e-24','1e-20','1e-16','1e-12','1e-8', '1e-4', 0.01, .1, 1)
  ) +
  # scale_y_continuous(name='', 
  #                    # trans = function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))),
  #                    breaks = (function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))))(log(10^cutpoints)),
  #                    labels=lab.cutpoints) +
  coord_flip()+
  theme_bw() 
  # theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
# load('data/cleaned/data_input.Rdata')
# tb_dis = !data_19EI$other_dis_dx
# X = cbind(Xd, Xc, Xc[,8]^2)
# colnames(X) = NULL
# X = X[,c(1:5,17,6:7,9:12,14,15,13,16,20,18,19,8)] |> as.data.frame()
# X <- Hmisc::`label<-`(X, self=FALSE, 
#                       value = c(
#                         'HIV +',
#                         'TB-suggested symptoms',
#                         'Local motor palsy',
#                         'Cranial nerve palsy',
#                         'Past noticed TB contact',
#                         'Reversed GCS',
#                         'Pulmonary TB/X-ray',
#                         'Miliary TB/X-Ray',
#                         'Age',
#                         'Days from onset',
#                         'Paired blood glucose',
#                         'CSF Glucose',
#                         'CSF Protein',
#                         'CSF Lactate',
#                         'CSF Lymphocyte count',
#                         'CSF WBC count',
#                         '(CSF WBC count)^2',
#                         'CSF Eosinophil count',
#                         'CSF RBC count',
#                         'Cryptococcus + '
#                       ))
# X$TBM = tb_dis
# dd = rms::datadist(X)
# options(datadist='dd')
# pseudo.fit = rms::lrm(TBM ~ ., data=X)
# pseudo.fit$coefficients = a_plot$data$m
# nomo = rms::nomogram(pseudo.fit)
save(z_plot, z_plot_logit, a_plot, b_plot, ab_plot, wbc_plot, file = 'export/m3_plot.Rdata')
