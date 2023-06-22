####### Imputation check using method modified from Cattram et.al ##
# Author: Trinh Dong
# Email: trinhdhk@oucru.org
####################################################################

# Params to be changed
# model = load(outputs/fit_object_with_r=1&k=1)
m = model$outputs
load('data/cleaned/data_input.Rdata')

# Some misc func to cope with composite variables #######
obs = with(data_19EI, obs_smear + obs_mgit + obs_xpert > 0)
test = with(data_19EI, csf_smear + csf_mgit + csf_xpert > 0)
X = rstan::extract(m, 'X')$X
fill_x <- function(x, obs, pred_fn, model = m, full=FALSE){
  model=model
  which_obs <- which(as.logical(obs))
  x_obs <- x[which_obs]
  x_pred <- vector('list', length(x))
  which_not_obs <- setdiff(seq_along(x), which_obs)
  x_obs <- x[which_obs]
  pred <- pred_fn(model=model, which_not_obs)
  j <- 0
  x_pred <- lapply(seq_along(x), function(i){
    if (obs[[i]]) return(x[[i]])
    j <<- j + 1
    pred[,j]
  })
  if (full) return(list(obs=x, pred=do.call(cbind,x_pred)))
  x_pred = x_pred[which_obs]
  list(obs = x_obs, pred = do.call(cbind, x_pred))
}

fill_x_multi <- function(x, obs, pred_fn, model = m, full=FALSE){
  model = model
  which_obs <- lapply(seq_len(dim(obs)[2]), function(i) which(obs[,i]==1))
  which_not_obs <-  lapply(seq_len(dim(obs)[2]), function(i) which(obs[,i]==0))
  x_obs <- lapply(seq_len(dim(obs)[2]), function(i) x[obs[,i]==1, i])
  x_pred <- vector('list', length(x_obs))
  for (k in seq_along(x_pred)) x_pred[[k]] <- vector('list', dim(x)[[1]])

  for (k in seq_along(x_obs)){
    pred <- pred_fn(model=model, which_not_obs[[k]])
    j <- 0
    x_pred[[k]] <- lapply(seq_len(dim(x)[[1]]), function(i){
      if (obs[i, k]) return(x[i, k])
      j <<- j + 1
      sapply(pred, function(ph) ph[j, k])
    })
  }
  for (k in seq_along(x_pred)) x_pred[[k]] = do.call(cbind, x_pred[[k]])
  if (!full) for (k in seq_along(x_pred)) x_pred[[k]] = x_pred[[k]][,which_obs[[k]]]
  list(obs = x_obs, pred = x_pred)
}

get_par = function(par, model) {
  par = rstan::extract(model, par)[[par]]
  par
}
`%.*%` = function(x, y){
  sapply(x, "*", y)
}

cs_pred = function(model, id){
  cs_a0 = get_par('cs_a0', model)
  cs_a = get_par('cs_a', model)
  cs_p = get_par('cs_p', model)
  
  
  lapply(seq_len(dim(cs_a0)[1]), 
         function(i){
           sapply(id, function(j) {
           mu_cs =  cs_a0[i] + cs_a[1]%*%X[i,j,1] + cs_a[2]%*%X[i,j,11] + cs_a[3]%*%obs[j] + cs_a[4]%*%test[j]
           cs  = rbinom(1, size=1, prob = plogis(mu_cs))
           if (cs == 0) return(rep(0,3))
           c(rbinom(1,1,cs_p[1]), rbinom(1,1,cs_p[2]), rbinom(1,1,cs_p[3]))
          }) |> t()
         })
}

cs_pred = fill_x_multi(Td[,1:3], obs_Td[,1:3], cs_pred, full = T)

#### Check zones ----------

# Check zone -------
color = list(pred = '#2c7fb8', obs = '#f03b20')
# response propensity check
prob_missingness = function(obs, completed){
  if (is.matrix(completed))
    return(
      apply(completed,1,
            \(x){
              fit <- glm(as.logical(obs)~x,family=binomial())
              predict(fit, type='response')
            })
    )
  
  sapply(completed, 
         \(x){
           fit <- glm(as.logical(obs)~as.matrix(x),family=binomial())
           predict(fit, type='response')
         })
}

resid_missingess = function(y, e, p){
  if (!is.matrix(p)) dim(p)=dim(y)
  family = if (identical(sort(unique(as.numeric(y))), c(0,1))){
    binomial
  } else gaussian
  sapply(seq_len(dim(y)[1]), 
         \(i) {
           fit = glm(y[i,]~e[i,], family=family)
           residuals(fit, type='response')
         }) |> t() 
}

# save(list=ls(), file='.cache/impute_chk.Rdata')
library(ggplot2)
library(ggfx)
library(data.table)

hiv_miss = prob_missingness(obs_Xd[Td[,7]==1,1], 
                            lapply(1:1000,
                                   function(i) 
                                     cbind(X[i,Td[,7]==1,c(11)],
                                           data_19EI[Td[,7]==1, .(obs_smear, csf_smear, csf_mgit, csf_xpert)]
                                       ))) |> t()

hiv_plot = list()

hiv_plot$plot1 = 
  ggplot() + 
   geom_hex(aes(y=X[201:700,Td[,7]==1,1] |> c(), x=hiv_miss[201:700,]|>c(), fill = rep(as.logical(obs_Xd[Td[,7]==1,1]),each=500), alpha=..ncount..)) + #fill = rep(as.logical(obs_Xd[Td[,7]==1,1]), each=2)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=X[201:700,Td[,7]==1,1]|>c(), x=hiv_miss[201:700,]|>c(), color = rep(as.logical(obs_Xd[Td[,7]==1,1]),each=500), fill = rep(as.logical(obs_Xd[Td[,7]==1,1]), each=500)), alpha=.5, span=.5, method='loess',se=F)) + 
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,0.7), trans='log10') + guides(alpha='none', fill='none')+
  scale_y_continuous(limits=c(0,1), oob=scales::squish)+
  ylab("P(HIV)")+
  xlab('P(Response)') +
  theme_bw()

hiv_plot$plot2 = 
  ggplot()+
  geom_density(aes(x=resid_missingess(X[1:500,Td[,7]==1,1], qlogis(hiv_miss[1:500,]), rep(as.logical(obs_Xd[Td[,7]==1,1]),each=500))|>c(), fill=rep(as.logical(obs_Xd[Td[,7]==1,1]),each=500)), color='transparent', alpha=.5, adjust=.5) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

hiv_plot$plot3 = 
  ggplot()+
  geom_hex(aes(y=resid_missingess(X[601:800,Td[,7]==1,1], hiv_miss[601:800,], rep(as.logical(obs_Xd[Td[,7]==1,1]),each=200))|>c(), x = hiv_miss[601:800,]|>c(), fill=rep(as.logical(obs_Xd[Td[,7]==1,1]),each=200), alpha = ..count..)) + 
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(X[601:800,Td[,7]==1,1], hiv_miss[601:800,], rep(as.logical(obs_Xd[Td[,7]==1,1]),each=200))|>c(), x = hiv_miss[601:800,]|>c(), color=rep(as.logical(obs_Xd[Td[,7]==1,1]),each=200))), method='loess', span=.5, se=FALSE)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  xlab('P(Response)') +
  ylab('Residual') +
  theme_bw()

# Clinical symptoms
cs_miss = prob_missingness(obs_Xd[,2],
                           lapply(seq_len(dim(X)[1]),
                                  function(i) 
                                    cbind(X[i,,c(1,3,11)], #cs_pred$pred[[2]][i,], cs_pred$pred[[3]][i,],
                                          obs, 
                                          (data_19EI[, .(csf_smear, csf_mgit, csf_xpert)])
                                          ))) |> t()
cs_plot = list()

cs_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=X[1:200,,2]|>c(), x=cs_miss[1:200,]|>c(), fill = rep(as.logical(obs_Xd[,2]),each=200), alpha=..count..)) + #,fill = rep(as.logical(obs_Xd[,2]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=X[1:200,,2]|>c(), x=cs_miss[1:200,]|>c(), color = rep(as.logical(obs_Xd[,2]),each=200),fill = rep(as.logical(obs_Xd[,2]), each=200)), alpha=.5,method='loess',span=.5, se=FALSE)) + 
  scale_x_continuous(trans='logit')+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  ylab("P(TB-suggestive symptoms)")+
  xlab('P(Response)') +
  theme_bw()

cs_plot$plot2 =
  ggplot()+
  geom_density(aes(x=resid_missingess(X[1:500,,2], cs_miss[1:500,], rep(as.logical(obs_Xd[,2]),each=500)) |> c(), fill=rep(as.logical(obs_Xd[,2]),each=500)), color='transparent', alpha=.5, adjust=1) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

cs_plot$plot3 =
  ggplot()+
  geom_hex((aes(y=resid_missingess(X[1:200,,2], cs_miss[1:200,], rep(as.logical(obs_Xd[,2]),each=200)) |> c(), x = cs_miss[1:200,]|>c(), fill=rep(as.logical(obs_Xd[,2]),each=200), alpha=..count..))) +
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(X[1:200,,2], cs_miss[1:200,], rep(as.logical(obs_Xd[,2]),each=200)) |> c(), x = cs_miss[1:200,] |> c(), color=rep(as.logical(obs_Xd[,2]),each=200))), method='loess',se=FALSE)) +
  scale_x_continuous(trans='logit') +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  xlab('P(Response)') +
  ylab('Residual') +
  theme_bw()


# Motor palsy
mp_miss = prob_missingness(obs_Xd[,3], 
                           lapply(1:500,
                                  function(i) 
                                    cbind(X[i,,c(1,2)], obs,
                                          (data_19EI[, .(csf_smear, csf_mgit, csf_xpert)])
                                          ))) |> t()

mp_plot = list()


mp_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=X[1:200,,3]|>c(), x=mp_miss[1:200,]|>c(), fill = rep(as.logical(obs_Xd[,3]),each=200), alpha=..count..)) + #,fill = rep(as.logical(obs_Xd[,3]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=X[1:200,,3]|>c(), x=mp_miss[1:200,]|>c(), color = rep(as.logical(obs_Xd[,3]),each=200),fill = rep(as.logical(obs_Xd[,3]), each=200)), alpha=.7,method='loess',span=.7, se=FALSE)) + 
  # scale_x_continuous(trans='logit', )+
  scale_y_continuous(limits=c(0,1), oob=scales::squish) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  ylab("P(Focal neuro-deficit)")+
  xlab('P(Response)') +
  theme_bw()

mp_plot$plot2 = 
  ggplot() +
  geom_density(aes(x=resid_missingess(X[1:500,,3], mp_miss[1:500,], rep(as.logical(obs_Xd[,3]),each=500))|>c(), fill=rep(as.logical(obs_Xd[,3]),each=500)), color='transparent', alpha=.5, adjust=.5) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  #xlim(-1, 1)+
  theme_bw()

mp_plot$plot3 = 
  ggplot()+
  geom_hex((aes(y=resid_missingess(X[1:200,,3], mp_miss[1:200,], rep(as.logical(obs_Xd[,3]),each=200))|>c(), x = mp_miss[1:200,]|>c(), fill=rep(as.logical(obs_Xd[,3]),each=200), alpha=..count..))) +
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(X[1:200,,3], mp_miss[1:200,], rep(as.logical(obs_Xd[,3]),each=200))|>c(), x = mp_miss[1:200,]|>c(), color=rep(as.logical(obs_Xd[,3]),each=200))),method='loess', se=FALSE)) +
  # scale_x_continuous(trans='logit', breaks = c(.5, .9, .99, .999, .9999)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  xlab('P(Response)') +
  ylab('Residual') +
  theme_bw()


# Illness days
id_miss = prob_missingness(obs_Xc[,1], 
                           lapply(seq_len(dim(X)[1]),
                                  function(i) 
                                    cbind(X[i,,c(1:3,18)], obs, test
                                          ))) |> t()
id_plot = list()
id_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=X[1:200,,11] |> c(), x=id_miss[1:200,] |> c(), fill = rep(as.logical(obs_Xc[,1]),each=200), alpha=..count..)) + #, fill = rep(as.logical(obs_Xc[,1]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=X[1:200,,11] |> c(), x=id_miss[1:200,] |> c(), color = rep(as.logical(obs_Xc[,1]),each=200),fill = rep(as.logical(obs_Xc[,1]), each=200)), alpha=.5, method='loess',span=.5, se=FALSE)) + 
  scale_x_continuous(trans='logit', breaks=c(0.4, .8, .95, .99, .999))+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  ylab(latex2exp::TeX("$log_{2}$(Symptom Duration, days)")) +
  xlab('P(Response)') +
  theme_bw()

id_plot$plot2 = 
  ggplot()+
  geom_density(aes(x=resid_missingess(X[1:500,,11], id_miss[1:500,], rep(as.logical(obs_Xc[,1]),each=500)) |> c(), fill=rep(as.logical(obs_Xc[,1]),each=500)), color='transparent', alpha=.5, adjust=2) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

id_plot$plot3 = 
  ggplot()+
  geom_hex((aes(y=resid_missingess(X[1:200,,11], id_miss[1:200,], rep(as.logical(obs_Xc[,1]),each=200)) |> c(), x = id_miss[1:200,] |> c(), fill=rep(as.logical(obs_Xc[,1]),each=200), alpha=..count..))) + 
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(X[1:200,,11], id_miss[1:200,], rep(as.logical(obs_Xc[,1]),each=200)) |> c(), x = id_miss[1:200,] |> c(), color=rep(as.logical(obs_Xc[,1]),each=200))), method='loess', span=.5,se=FALSE)) +
  scale_x_continuous(trans='logit', breaks=c(0.4, .8, .95, .99, .999)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  xlab('P(Response)') +
  ylab('Residual') +
  theme_bw()

# GCS
gcs_miss = prob_missingness(obs_Xc[,8], 
                            lapply(seq_len(dim(X)[1]),
                                   function(i) 
                                     cbind(X[i,,c(1,11,12:18)], obs, test
                                           ))) |> t()

gcs_plot = list()
ft = Xc[,8] > -1
gcs_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=(X[1:200,ft,18] + rnorm(length(X[1:200,ft,18]),0,.01)) |> c(), x=gcs_miss[1:200,ft] |> c(), fill = rep(as.logical(obs_Xc[ft,8]),each=200), alpha=..count..)) + #,fill = rep(as.logical(obs_Xc[,8]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=X[1:200,ft,18] |> c(), x=gcs_miss[1:200,ft] |> c(), color = rep(as.logical(obs_Xc[ft,8]),each=200),fill = rep(as.logical(obs_Xc[ft,8]), each=200)), alpha=.5,method='loess',span=.5, se=FALSE)) + 
  scale_x_continuous(trans='logit', breaks = c(.4, .8, .95, .99, .999))+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  ylab("RGCS")+
  xlab('P(Response)') +
  theme_bw()


gcs_plot$plot2 = 
  ggplot()+
  geom_density(aes(x=resid_missingess(12-3*X[1:500,ft,18], gcs_miss[1:500,ft], rep(as.logical(obs_Xc[ft,8]),each=500)) |> c(), fill=rep(as.logical(obs_Xc[ft,8]),each=500)), color='transparent', alpha=.5, adjust=1) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

gcs_plot$plot3 = 
  ggplot()+
  geom_hex(aes(y=resid_missingess(X[1:200,ft,18], gcs_miss[1:200,ft], rep(as.logical(obs_Xc[ft,8]),each=200)) |> c(), x = gcs_miss[1:200,ft] |> c(), fill=rep(as.logical(obs_Xc[ft,8]),each=200), alpha=..count..)) +
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(X[1:200,ft,18], gcs_miss[1:200,ft], rep(as.logical(obs_Xc[ft,8]),each=200)) |> c(), x = gcs_miss[1:200,ft] |> c(), color=rep(as.logical(obs_Xc[ft,8]),each=200))),method='loess', se=FALSE)) +
  scale_x_continuous(trans='logit', breaks = c(.4, .8, .95, .99, .999)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.2,.7), trans='log10') + guides(alpha='none', fill='none')+
  xlab('P(Response)') +
  ylab('Residual') +
  theme_bw()


grobs = function(plt){
  sapply(plt, ggplot2::ggplotGrob, simplify = FALSE)  
}

hiv_plot = grobs(hiv_plot)
id_plot = grobs(id_plot)
gcs_plot = grobs(gcs_plot)
cs_plot = grobs(cs_plot)
mp_plot = grobs(mp_plot)
save(hiv_plot, id_plot, gcs_plot, cs_plot, mp_plot, file='export/impute_chk_plot_hex.Rdata')

