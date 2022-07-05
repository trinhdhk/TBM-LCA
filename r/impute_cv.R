model = m3$model
holdout = m3$folds$holdout
keptin = m3$folds$keptin
dat = m3$folds$inputs
load('data/cleaned/data_input.Rdata')

obs = with(data_19EI, obs_smear + obs_mgit + obs_xpert > 0)
test = with(data_19EI, csf_smear + csf_mgit + csf_xpert > 0)
fill_x_holdout <- function(x, obs, pred_fn, holdout, m = m3, full=FALSE){
  model <- m$model
  K <- m$n_fold
  N <- m$n_rep
  which_obs <- which(as.logical(obs))
  x_obs <- x[which_obs]
  # which_obs <- which(!is.na(x))ob
  x_pred <- vector('list', length(x))
  
  i = 1
  cli::cli_progress_bar(total = length(holdout))
  for (h in holdout){
    # cat(i, ' ')
    h = which(as.logical(h))
    # x_heldout = x[h]
    pred_heldout = pred_fn(model[[i]], h)
    j = 1
    for (hh in h) {
      x_pred[[hh]] = c(x_pred[[hh]], pred_heldout[,j])
      j = j+1
    }
    i = i+1
    cli::cli_progress_update()
  }
  # browser()
  if (full) return(list(obs=x, pred=do.call(cbind,x_pred)))
  x_pred = x_pred[which_obs]
  list(obs = x_obs, pred = do.call(cbind, x_pred))
}

fill_x_holdout_multi <- function(x, obs, pred_fn, holdout, m = m3, full=FALSE){
  model <- m$model
  K <- m$n_fold
  N <- m$n_rep
  which_obs <- lapply(seq_len(dim(obs)[2]), function(i) which(obs[,i]==1))
  # browser()
  x_obs <- lapply(seq_len(dim(obs)[2]), function(i) x[obs[,i]==1, i])
  # which_obs <- which(!is.na(x))ob
  x_pred <- vector('list', length(x_obs))
  for (k in seq_along(x_pred)) x_pred[[k]] <- vector('list', dim(x)[[1]])
  # for (k in seq_len(length(x_obs))) x_pred[[k]] = vector('numeric', length=)

  i = 1
  cli::cli_progress_bar(total = length(holdout))
  for (h in holdout){
    # cat(i, ' ')
    h = which(as.logical(h))
    # x_heldout = x[h]
    pred_heldout = pred_fn(model[[i]], h)

    for (k in seq_len(length(x_obs))){
      j = 1
      for (hh in h) {
        x_pred[[k]][[hh]] = c(x_pred[[k]][[hh]], sapply(pred_heldout, function(ph) ph[j, k]))
        j = j+1
      }
    }

    i = i+1
    cli::cli_progress_update()
  }
  # browser()
  for (k in seq_along(x_pred)) x_pred[[k]] = do.call(cbind, x_pred[[k]])
  # return(list(obs=x, pred=x_pred))
  if (!full) for (k in seq_along(x_pred)) x_pred[[k]] = x_pred[[k]][,which_obs[[k]]]
  # browser()
  list(obs = x_obs, pred = x_pred)
}

get_par = function(par, model) {
  # id = id
  par = rstan::extract(model, par)[[par]]
  par
  # abind::asub(par, id, 1)
}
`%.*%` = function(x, y){
  sapply(x, "*", y)
}

# impute_misc=new.env()
# impute_misc$calib <- function(x){
#   require(ggplot2)
#   ggplot(mapping=aes(x = apply(x$pred,2,mean)-x$obs, y=x$obs)) + 
#     geom_smooth()
# }

#hiv
hiv_pred = function(model, id){
  hiv_a0 = get_par('HIV_a0', model)
  hiv_a = get_par('HIV_a', model)
  # browser()
  sapply(seq_len(dim(hiv_a0)[1]),
        function(i){
          sapply(id, function(j){
            plogis(hiv_a0[i] + hiv_a[i,1]*obs[j] + hiv_a[i,2]*test[j])
          })
          # rbinom(length(id), 1, prob = pnorm(hiv_a0[i]))
        })  |> t()
  # matrix(mean(pnorm(hiv_a0)), ncol=length(id), nrow = dim(hiv_a0)[1])
}

obs_Xd[Td[,7]==0,1] = 1
hiv = fill_x_holdout(Xd[Td[,7]==1,1], obs_Xd[Td[,7]==1,1], hiv_pred, lapply(holdout, \(x) x[Td[,7]==1]), full=F)
# hiv2 = fill_x_holdout(Xd[Td[,7]==1,1], obs_Xd[Td[,7]==1,1], hiv_pred, lapply(holdout, \(x) x[Td[,7]==1]), full=T)
hiv3 = fill_x_holdout(Xd[,1], obs_Xd[,1], hiv_pred, holdout, full=T)
hiv2 = list(obs = hiv3$obs, pred=apply(hiv3$pred, 2, function(x) sapply(x, rbinom, n=1, size=1)))
# hiv$pred[,as.logical(obs_Xd[,1])]=Xd[as.logical(obs_Xd[,1]),1]

# gcs
gcs_pred = function(model, id){
  L_sigma_gcs = get_par('L_sigma_gcs', model)
  L_omega_gcs = get_par('L_Omega_gcs', model)
  gcs_a0 = get_par('gcs_a0', model)
  gcs_a = get_par('gcs_a', model)
  future::plan(future::multisession, workers=20)
  future.apply::future_sapply(seq_len(dim(L_omega_gcs)[1]), 
         FUN = function(i){
           L_Sigma_gcs = diag(L_sigma_gcs[i,]) %*% L_omega_gcs[i,,]
           sapply(id, function(j) {
             gcs_mu = gcs_a0[i] + gcs_a[i,,1]%.*%hiv2$pred[j] + gcs_a[i,,2]%.*%obs[j] + gcs_a[i,,3]%.*%test[j] 
             GCS = sapply(id, 
                          function(...) {
                            gcs = LaplacesDemon::rmvnc(1, mu = gcs_mu, U = t(L_Sigma_gcs))
                            while (any(gcs < 0 | gcs > 1)) gcs =  LaplacesDemon::rmvnc(1, mu = gcs_mu, U = t(L_Sigma_gcs))
                            gcs
                          }) |> t()
             round(3*GCS[,1] + 5*GCS[,2] + 4*GCS[,3] - 3)/3;
           })
         }, future.seed=TRUE) |> t()
}
# gcs = fill_x_holdout(Xc[,8], obs_Xc[,8], gcs_pred, holdout)
gcs2 = fill_x_holdout(Xc[,8], obs_Xc[,8], gcs_pred, holdout, full=T)


#id
id_pred = function(model, id){
  id_a0 = get_par('id_a0', model)
  id_a = get_par('id_a', model)
  id_sigma = get_par('id_sigma', model)
  # browser()
  sapply(seq_len(dim(id_a0)[1]), function(i) {
    sapply(id, function(j) {
      mu = id_a0[i] + id_a[i,1]*hiv2$pred[j] + id_a[i,2]*obs[j] + id_a[i,3]*test[j] #*hiv$pred[i,j]
      rnorm(1, mean = mu, sd = id_sigma[i])
      # mu
    })
  }) |> t()
} 

# id = fill_x_holdout(Xc[,1], obs_Xc[,1], id_pred, holdout)
id2 = fill_x_holdout(Xc[,1], obs_Xc[,1], id_pred, holdout, full=TRUE)

# clin_symptoms
cs_pred = function(model, id, binary = FALSE){
  cs_a0 = get_par('cs_a0', model)
  cs_a = get_par('cs_a', model)
  L_omega_cs = get_par('L_Omega_cs', model)
  # browser()
  sapply(seq_len(dim(cs_a0)[1]),
         function(i){
           mu_cs =  cs_a0[i,] + cs_a[i,,1]%.*%hiv2$pred[id,1] + cs_a[i,,2]%.*%id2$pred[i,id] + cs_a[i,,3]%.*%obs[id] + cs_a[i,,4]%.*%test[id]
           cs = apply(mu_cs, 1,
                      function(mu) mvtnorm::rmvnorm(1, mean = mu, sigma = L_omega_cs[i,,] %*% t(L_omega_cs[i,,]))) |> pnorm()
           if (binary) {
             cs = apply(cs, 2, \(x) sapply(x, rbinom, n=1, size=1))
             (cs[1,] | cs[2,] | cs[3,]) |> as.numeric()
           } else {
             1 - (1-cs[1,])*(1-cs[2,])*(1-cs[3,])
           }
         }) |> t()
}

cs_pred_binary = purrr::partial(cs_pred, binary=TRUE)
cs_pred_prob = purrr::partial(cs_pred, binary=FALSE)


# clin_symptoms = fill_x_holdout(Xd[,2], obs_Xd[,2], cs_pred_prob, holdout)
clin_symptoms2 = fill_x_holdout(Xd[,2], obs_Xd[,2], cs_pred_binary, holdout, full=T)
clin_symptoms3 = fill_x_holdout(Xd[,2], obs_Xd[,2], cs_pred_prob, holdout, full=T)

# # clin_symtomps compartments -------
# cs1_pred = function(model, id){
#   cs_a0 = get_par('cs_a0', model)
#   cs_a = get_par('cs_a', model)
#   L_omega_cs = get_par('L_Omega_cs', model)
#   # browser()
#   sapply(seq_len(dim(cs_a0)[1]),
#          function(i){
#            mu_cs =  cs_a0[i,] + cs_a[i,,1]%.*%Xd[id,1] + cs_a[i,,2]%.*%id2$pred[i,id]
#            cs = apply(mu_cs, 1,
#                       function(mu) mvtnorm::rmvnorm(1, mean = mu, sigma = L_omega_cs[i,,] %*% t(L_omega_cs[i,,])) > 0) 
#            (cs[1,]) |> as.numeric()
#          }) |> t()
# }
# 
# cs2_pred = function(model, id){
#   cs_a0 = get_par('cs_a0', model)
#   cs_a = get_par('cs_a', model)
#   L_omega_cs = get_par('L_Omega_cs', model)
#   # browser()
#   sapply(seq_len(dim(cs_a0)[1]),
#          function(i){
#            mu_cs =  cs_a0[i,] + cs_a[i,,1]%.*%Xd[id,1] + cs_a[i,,2]%.*%id2$pred[i,id]
#            cs = apply(mu_cs, 1,
#                       function(mu) mvtnorm::rmvnorm(1, mean = mu, sigma = L_omega_cs[i,,] %*% t(L_omega_cs[i,,])) > 0) 
#            (cs[2,]) |> as.numeric()
#          }) |> t()
# }
# 
# cs3_pred = function(model, id){
#   cs_a0 = get_par('cs_a0', model)
#   cs_a = get_par('cs_a', model)
#   L_omega_cs = get_par('L_Omega_cs', model)
#   # browser()
#   sapply(seq_len(dim(cs_a0)[1]),
#          function(i){.
#            mu_cs =  cs_a0[i,] + cs_a[i,,1]%.*%Xd[id,1] + cs_a[i,,2]%.*%id2$pred[i,id]
#            cs = apply(mu_cs, 1,
#                       function(mu) mvtnorm::rmvnorm(1, mean = mu, sigma = L_omega_cs[i,,] %*% t(L_omega_cs[i,,]))) |> pnorm()
#            cs = apply(cs, 2, \(x) sapply(xm, rbinom, n=1, size=1))
#            (cs[3,]) |> as.numeric()
#          }) |> t()
# }
# 
# cs1 = fill_x_holdout(Td[,1], obs_Td[,1], cs1_pred, keptin[1:20])
# cs2 = fill_x_holdout(Td[,2], obs_Td[,2], cs2_pred, keptin[1:20])
# cs3 = fill_x_holdout(Td[,3], obs_Td[,3], cs3_pred, keptin[1:20])
# #---------

# motor palsy
mp_pred = function(model, id, binary=FALSE){
  mp_a0 = get_par('mp_a0', model)
  mp_a = get_par('mp_a', model)
  L_omega_mp = get_par('L_Omega_mp', model)
  # browser()
  sapply(seq_len(dim(mp_a0)[1]),
         function(i){
           mu_mp =  mp_a0[i,] + mp_a[i,,1]%.*%hiv2$pred[id,1] + mp_a[i,,2]%.*%id2$pred[i,id] + mp_a[i,,3]%.*%obs[id] + mp_a[i,,4]%.*%test[id]
           # mu_mp =  mp_a0[i,] + mp_a[i,,1]%.*%hiv2$pred[id,1] + mp_a[i,,2]%.*%id2$pred[i,id] + mp_a[i,,3]%.*%obs + mp_a[i,,4]%*.*%test
           mp = apply(mu_mp, 1,
                      function(mu) mvtnorm::rmvnorm(1, mean = mu, sigma = L_omega_mp[i,,] %*% t(L_omega_mp[i,,]))) |> pnorm()
           
           if (binary) {
             mp = apply(mp, 2, \(x) sapply(x, rbinom, n=1, size=1))
             (mp[1,] | mp[2,] | mp[3,]) |> as.numeric()
           } else {
             1 - (1-mp[1,])*(1-mp[2,])*(1-mp[3,])
           }
           
         }) |> t()
}

mp_pred_binary = purrr::partial(mp_pred, binary=TRUE)
mp_pred_prob = purrr::partial(mp_pred, binary=FALSE)

# motor_palsy = fill_x_holdout(Xd[,3], obs_Xd[,3], mp_pred_prob, holdout)
motor_palsy2 = fill_x_holdout(Xd[,3], obs_Xd[,3], mp_pred_binary, holdout, full=T)
motor_palsy3 = fill_x_holdout(Xd[,3], obs_Xd[,3], mp_pred_prob, holdout, full=T)

# CSF
csf_pred = \(model, id){
  L_omega_csf   = get_par('L_Omega_csf', model)
  L_sigma_csf   = get_par('L_sigma_csf', model)
  csf_a0   = get_par('csf_a0', model)
  csf_a    = get_par('csf_a', model)
  # print(dim(csf_a))

  lapply(seq_len(dim(L_sigma_csf)[1]),
         function(i){
           L_Sigma_csf = diag(L_sigma_csf[i,]) %*% L_omega_csf[i,,]
           # mvtnorm::rmvnorm(N, sigma = L_Sigma_csf %*% t(L_Sigma_csf), method = 'chol')
           sapply(id,
                  function(j) {
                    csf_mu = csf_a0[i] + csf_a[i,,1]%.*%hiv2$pred[j] + csf_a[i,,2]%.*%obs[j] + csf_a[i,,3]%.*%test[j]
                    # print(csf_mu)
                    # printL_Sigma_csf)
                    LaplacesDemon::rmvnc(1, mu = csf_mu, U = t(L_Sigma_csf)) 
                    
                  }) |> t()
         })
}

csf = fill_x_holdout_multi(Xc[,2:7], obs_Xc[,2:7],csf_pred, holdout, full=T)

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

resid_missingess = function(y, e){
  sapply(seq_len(dim(y)[1]), 
         \(i) {
           fit = lm(y[i,]~e[i,])
           residuals(fit)
         }) |> t()
}

# resid_missingness.binary = function(Py, e, completed.quan, completed.qual){
#   if (is.matrix(completed.quan))
#     return(
#       sapply(seq_len(dim(completed.quan)[1]),
#             \(x){
#               pca <- PCAmixdata::PCAmix(X.quanti=completed.quan[i,], X.quali = completed.qual[i,], rename.level=TRUE)
#               dim1 <- pca$ind$coord[,1]
#               fit <- lm(dim1~Py*e)
#               predict(fit)
#             })
#     )
#   
#   sapply(seq_len(length(completed.quan)),
#          \(x){
#            pca <- PCAmixdata::PCAmix(X.quanti=completed.quan[[i]], X.quali = completed.qual[[i]], rename.level=TRUE)
#            dim1 <- pca$ind$coord[,1]
#            fit <- lm(dim1~Py*e)
#            predict(fit)
#          })
# }

save(list=ls(), file='.cache/impute_cv.Rdata')
library(ggplot2)
library(ggfx)
library(data.table)
## HIV

hiv_miss = prob_missingness(obs_Xd[Td[,7]==1,1], 
                           lapply(seq_len(dim(gcs2$pred)[1]),
                                  function(i) 
                                    cbind(gcs2$pred[i,Td[,7]==1], clin_symptoms2$pred[i,Td[,7]==1], motor_palsy2$pred[i,Td[,7]==1], id2$pred[i,Td[,7]==1],
                                          csf$pred[[1]][i,Td[,7]==1], csf$pred[[2]][i,Td[,7]==1], csf$pred[[3]][i,Td[,7]==1],
                                          csf$pred[[4]][i,Td[,7]==1], csf$pred[[5]][i,Td[,7]==1], csf$pred[[6]][i,Td[,7]==1],
                                          (data_19EI[Td[,7]==1, .(csf_smear, csf_mgit, csf_xpert)])))) |> t()

hiv_plot = list()

hiv_plot$plot0 = 
  ggplot() +
    geom_density(aes(x=hiv$pred), color=color$pred, adjust=4) + geom_vline(aes(xintercept=mean(hiv$obs)), color=color$obs) +
    theme_bw() +
    xlab('Pr (HIV+)') + ylab('Density')

hiv_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=hiv3$pred[1:100,Td[,7]==1], x=hiv_miss[1:100,], fill = rep(as.logical(obs_Xd[Td[,7]==1,1]),each=100), alpha=..ncount..)) + #fill = rep(as.logical(obs_Xd[Td[,7]==1,1]), each=2)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=hiv3$pred[1:100,Td[,7]==1], x=hiv_miss[1:100,], color = rep(as.logical(obs_Xd[Td[,7]==1,1]),each=100), fill = rep(as.logical(obs_Xd[Td[,7]==1,1]), each=100)), alpha=.5, method='loess', se=F)) + 
  scale_x_continuous(trans='logit', breaks = c(.25, .75, .95, .99, .999), limits=c(0.2, .9999))+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.1,.75), trans='log10') + guides(alpha='none', fill='none')+
  ylab("Pr (HIV +)")+
  xlab('Pr (Response)') +
  theme_bw()

hiv_plot$plot2 = 
  ggplot()+
  geom_density(aes(x=resid_missingess(hiv3$pred[,Td[,7]==1], hiv_miss)[1:2000,], fill=rep(as.logical(obs_Xd[Td[,7]==1,1]),each=2000)), color='transparent', alpha=.5, adjust=.5) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

hiv_plot$plot3 = 
  ggplot()+
  geom_hex(aes(y=resid_missingess(hiv3$pred[,Td[,7]==1], hiv_miss)[1:100,], x = hiv_miss[1:100,], fill=rep(as.logical(obs_Xd[Td[,7]==1,1]),each=100), alpha = ..count..)) + 
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(hiv3$pred[,Td[,7]==1], hiv_miss)[1:100,], x = hiv_miss[1:100,], color=rep(as.logical(obs_Xd[Td[,7]==1,1]),each=100))), method='loess', se=FALSE)) +
  scale_x_continuous(trans='logit', breaks = c(.25, .75, .95, .99, .999), limits=c(0.2, .9999)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.1,.75), trans='log10') + guides(alpha='none')+
  xlab('Pr (Response)') +
  ylab('Residual') +
  theme_bw()


# Illness days
id_miss = prob_missingness(obs_Xc[,1], 
                           lapply(seq_len(dim(hiv3$pred)[1]),
                                  function(i) 
                                    cbind(hiv2$pred[i,], gcs2$pred[i,], clin_symptoms2$pred[i,],
                                          csf$pred[[1]][i,], csf$pred[[2]][i,], csf$pred[[3]][i,],
                                          csf$pred[[4]][i,], csf$pred[[5]][i,], csf$pred[[6]][i,],
                                          (data_19EI[, .(csf_smear, csf_mgit, csf_xpert)])))) |> t()
id_plot = list()
id_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=id2$pred[1:100,], x=id_miss[1:100,], fill = rep(as.logical(obs_Xc[,1]),each=100), alpha=..count..)) + #, fill = rep(as.logical(obs_Xc[,1]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=id2$pred[1:100,], x=id_miss[1:100,], color = rep(as.logical(obs_Xc[,1]),each=100),fill = rep(as.logical(obs_Xc[,1]), each=100)), alpha=.5, method='loess', se=FALSE)) + 
  scale_x_continuous(trans='logit', breaks=c(0.4, .8, .95, .99, .999))+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.1,.75), trans='log10') + guides(alpha='none', fill='none')+
  ylab("Days from onset")+
  xlab('Pr (Response)') +
  theme_bw()

id_plot$plot2 = 
  ggplot()+
  geom_density(aes(x=resid_missingess(id2$pred, id_miss)[1:500,], fill=rep(as.logical(obs_Xc[,1]),each=500)), color='transparent', alpha=.5, adjust=.5) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

id_plot$plot3 = 
  ggplot()+
  geom_point((aes(y=resid_missingess(id2$pred, id_miss)[1:100,], x = id_miss[1:100,], color=rep(as.logical(obs_Xc[,1]),each=100))), size=2, shape=21) + 
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(id2$pred, id_miss)[1:100,], x = id_miss[1:100,], color=rep(as.logical(obs_Xc[,1]),each=100))), method='loess', se=FALSE)) +
  scale_x_continuous(trans='logit', breaks=c(0.4, .8, .95, .99, .999)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Pr (Response)') +
  ylab('Residual') +
  theme_bw()

# GCS
gcs_miss = prob_missingness(obs_Xc[,1], 
                           lapply(seq_len(dim(hiv3$pred)[1]),
                                  function(i) 
                                    cbind(hiv2$pred[i,], id2$pred[i,], clin_symptoms2$pred[i,],
                                          csf$pred[[1]][i,], csf$pred[[2]][i,], csf$pred[[3]][i,],
                                          csf$pred[[4]][i,], csf$pred[[5]][i,], csf$pred[[6]][i,],
                                          (data_19EI[, .(csf_smear, csf_mgit, csf_xpert)])))) |> t()

gcs_plot = list()

gcs_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=gcs2$pred[1:100,] + rnorm(length(gcs2$pred[1:100,]),0,.05), x=gcs_miss[1:100,], fill = rep(as.logical(obs_Xc[,8]),each=100), alpha=..count..)) + #,fill = rep(as.logical(obs_Xc[,8]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=gcs2$pred[1:100,], x=gcs_miss[1:100,], color = rep(as.logical(obs_Xc[,8]),each=100),fill = rep(as.logical(obs_Xc[,8]), each=100)), alpha=.5,method='loess', span=1, se=FALSE)) + 
  scale_x_continuous(trans='logit', breaks = c(.4, .8, .95, .99, .999))+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.1,.75), trans='log10') + guides(alpha='none', fill='none')+
  ylab("Glascow coma score")+
  xlab('Pr (Response)') +
  theme_bw()

gcs_plot$plot2 = 
  ggplot()+
  geom_density(aes(x=resid_missingess(gcs2$pred, gcs_miss)[1:500,], fill=rep(as.logical(obs_Xc[,8]),each=500)), color='transparent', alpha=.5, adjust=3.5) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

gcs_plot$plot3 = 
  ggplot()+
  geom_jitter((aes(y=resid_missingess(gcs2$pred, gcs_miss)[1:100,], x = gcs_miss[1:100,], color=rep(as.logical(obs_Xc[,8]),each=100))), size=2, shape=21) +
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(gcs2$pred, gcs_miss)[1:100,], x = gcs_miss[1:100,], color=rep(as.logical(obs_Xc[,8]),each=100))),method='loess',span=1, se=FALSE)) +
  scale_x_continuous(trans='logit', breaks = c(.4, .8, .95, .99, .999)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Pr (Response)') +
  ylab('Residual') +
  theme_bw()

# Clinical symptoms
cs_miss = prob_missingness(obs_Xd[,2], 
                           lapply(seq_len(dim(hiv3$pred)[1]),
                                  function(i) 
                                    cbind(hiv2$pred[i,], id2$pred[i,], gcs2$pred[i,],
                                          csf$pred[[1]][i,], csf$pred[[2]][i,], csf$pred[[3]][i,],
                                          csf$pred[[4]][i,], csf$pred[[5]][i,], csf$pred[[6]][i,],
                                          (data_19EI[, .(csf_smear, csf_mgit, csf_xpert)])))) |> t()
cs_plot = list()

cs_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=clin_symptoms3$pred[1:100,], x=cs_miss[1:100,], fill = rep(as.logical(obs_Xd[,2]),each=100), alpha=..count..)) + #,fill = rep(as.logical(obs_Xd[,2]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=clin_symptoms3$pred[1:100,], x=cs_miss[1:100,], color = rep(as.logical(obs_Xd[,2]),each=100),fill = rep(as.logical(obs_Xd[,2]), each=100)), alpha=.5,method='loess', se=FALSE)) + 
  scale_x_continuous(trans='logit')+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.1,.75), trans='log10') + guides(alpha='none', fill='none')+
  ylab("Pr (TB symptoms)")+
  xlab('Pr (Response)') +
  theme_bw()

cs_plot$plot2 =
  ggplot()+
  geom_density(aes(x=resid_missingess(clin_symptoms3$pred, cs_miss)[1:500,], fill=rep(as.logical(obs_Xd[,2]),each=500)), color='transparent', alpha=.5, adjust=1) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

cs_plot$plot3 =
  ggplot()+
  geom_point((aes(y=resid_missingess(clin_symptoms3$pred, gcs_miss)[1:100,], x = cs_miss[1:100,], color=rep(as.logical(obs_Xd[,2]),each=100))), size=2, shape=21) +
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(clin_symptoms3$pred, gcs_miss)[1:100,], x = cs_miss[1:100,], color=rep(as.logical(obs_Xd[,2]),each=100))), method='loess',se=FALSE)) +
  scale_x_continuous(trans='logit') +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Pr (Response)') +
  ylab('Residual') +
  theme_bw()


# Motor palsy
mp_miss = prob_missingness(obs_Xd[,3], 
                           lapply(seq_len(dim(hiv3$pred)[1]),
                                  function(i) 
                                    cbind(hiv2$pred[i,], id2$pred[i,], gcs2$pred[i,],
                                          csf$pred[[1]][i,], csf$pred[[2]][i,], csf$pred[[3]][i,],
                                          csf$pred[[4]][i,], csf$pred[[5]][i,], csf$pred[[6]][i,],
                                          (data_19EI[, .(csf_smear, csf_mgit, csf_xpert)])))) |> t()

mp_plot = list()

mp_plot$plot1 = 
  ggplot() + 
  geom_hex(aes(y=motor_palsy3$pred[1:100,], x=mp_miss[1:100,], fill = rep(as.logical(obs_Xd[,3]),each=100), alpha=..count..)) + #,fill = rep(as.logical(obs_Xd[,3]), each=100)
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth(aes(y=motor_palsy3$pred[1:100,], x=mp_miss[1:100,], color = rep(as.logical(obs_Xd[,3]),each=100),fill = rep(as.logical(obs_Xd[,3]), each=100)), alpha=.5,method='loess',span=1, se=FALSE)) + 
  scale_x_continuous(trans='logit', breaks = c(.5, .9, .99, .999, .9999))+
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  scale_alpha_continuous(range=c(0.1,.75), trans='log10') + guides(alpha='none', fill='none')+
  ylab("Pr (Focal neuro-deficit)")+
  xlab('Pr (Response)') +
  theme_bw()

mp_plot$plot2 = 
  ggplot() +
  geom_density(aes(x=resid_missingess(motor_palsy3$pred, mp_miss)[1:500,], fill=rep(as.logical(obs_Xd[,3]),each=500)), color='transparent', alpha=.5, adjust=.5) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Residual') +
  ylab('Density') +
  theme_bw()

mp_plot$plot3 = 
  ggplot()+
  geom_point((aes(y=resid_missingess(motor_palsy3$pred, mp_miss)[1:50,], x = mp_miss[1:50,], color=rep(as.logical(obs_Xd[,3]),each=50))), size=2, shape=21) +
  ggfx::with_outer_glow(colour = '#ffffff', sigma=1, expand=3, x = geom_smooth((aes(y=resid_missingess(motor_palsy3$pred, mp_miss)[1:100,], x = mp_miss[1:100,], color=rep(as.logical(obs_Xd[,3]),each=100))),method='loess', se=FALSE)) +
  scale_x_continuous(trans='logit', breaks = c(.5, .9, .99, .999, .9999)) +
  scale_discrete_manual(c('color', 'fill'), name = '', breaks = c(FALSE, TRUE), labels=c('Imputed', 'Observed'), values = unlist(color) |> setNames(c('FALSE', 'TRUE')))+
  xlab('Pr (Response)') +
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
save(hiv_plot, id_plot, gcs_plot, cs_plot, mp_plot, file='export/impute_cv_plot_hex.Rdata')

