#!/usr/bin/env Rscript

## Script to check for coefficient divergence due to imputation, modified the method by Cattram et.al ####
## Author: Trinh Dong
## Email: trinhdhk@oucru.org
## Run via bash with r/scripts/checks/pp_check.R
##########################################################################################################
misc <- new.env()
source('r/include/functions.R', local=misc)

# m3 is load
cli::cli_alert_info('Load sampler')
model <- misc$my_stan_model('stan/m3.stan')

cli::cli_alert_info('Get pretrained model')
pretrained <- readRDS('outputs/m3_t00_b345678_q7_r1_k1_3103.RDS') #change to path/to/fit_object_with_r1k1.RDS

set.seed(10011963)  
simulate_data = 
  function(stanmodel){
    stopifnot(stanmodel$.META$all_params)
    m = stanmodel$outputs
    id = sample(seq_len((stanmodel$.META$iter - stanmodel$.META$warmup) * stanmodel$.META$chain)/stanmodel$.META$thin, 1)
    get_par = function(par, id) {
      id = id
      par = rstan::extract(m, par)[[par]]
      abind::asub(par, id, 1)
    }
    
    `%.*%` = function(x, y){
      sapply(x, "*", y)
    }
    hiv_a0 = get_par('HIV_a0', id)
    hiv_a  = get_par('HIV_a', id)
    
    id_a0  = get_par('id_a0', id) 
    id_a   = get_par('id_a', id) 
    id_sigma   = get_par('id_sigma',id) 
    
    cs_a0   = get_par('cs_a0', id) 
    cs_a    = get_par('cs_a', id)
    cs_p    = get_par('cs_p', id)
    
    mp_a0   = get_par('mp_a0', id)
    mp_a    = get_par('mp_a', id)
    mp_p    = get_par('mp_p', id)
    
    L_omega_csf   = get_par('L_Omega_csf', id) 
    L_sigma_csf   = get_par('L_sigma_csf', id) 
    csf_a0        = get_par('csf_a0', id)
    csf_a.        = get_par('csf_a', id)
    
    L_sigma_gcs   = get_par('L_sigma_gcs', id) 
    L_omega_gcs   = get_par('L_Omega_gcs', id) 
    gcs_a0        = get_par('gcs_a0', id)
    gcs_a         = get_par('gcs_a', id)
    
    load('data/cleaned/data_input.Rdata')
    N = dim(Xc)[1]
  
    obs = with(data_19EI, obs_smear + obs_mgit + obs_xpert) > 0
    test = with(data_19EI, csf_smear + csf_mgit + csf_xpert) > 0
    # HIV
    Xd[Td[,7]==1,1] = sapply(hiv_a0 + hiv_a[1]*obs[Td[,7]==1] + hiv_a[2]*test[Td[,7]==1],
                             function(iii) rbinom(1, size = 1, plogis(iii)))
    Xc[,1] = sapply(id_a0 + id_a[1]*Xd[,1] + id_a[2]*obs + id_a[3]*test,
                    function(iii) rnorm(1, mean = iii, sd = id_sigma))
    
    # Clinical symptoms
    mu_cs =  cs_a0 + cs_a[1]%*%Xd[,1] + cs_a[2]%*%Xc[,1] + cs_a[3]%*%obs + cs_a[4]%*%test
    Xd[,2] = sapply(mu_cs, function(iii) rbinom(1, size=1, plogis(iii)))
    # Motor palsy
    
    mu_mp =  mp_a0 + mp_a[1]%*%Xd[,1] + mp_a[2]%*%Xc[,1] + mp_a[3]%*%obs + mp_a[4]%*%test
    Xd[,3] = sapply(mu_mp, function(iii) rbinom(1, size=1, plogis(iii)))
    
    # GCS
    L_Sigma_gcs = diag(L_sigma_gcs) %*% L_omega_gcs
    mu_gcs = gcs_a0 + gcs_a[,1] %.*% Xd[,1] + gcs_a[,2] %.*% obs + gcs_a[,3] %.*% test
    
    GCS = sapply(seq_len(N), 
                 function(n) {
                   gcs = LaplacesDemon::rmvnc(1, mu = mu_gcs[n,], U = t(L_Sigma_gcs))
                   while (any(gcs < 0 | gcs > 1)) gcs =  LaplacesDemon::rmvnc(1, mu = mu_gcs[n,], U = t(L_Sigma_gcs))
                   gcs
                 }) |> t()
    
    Xc[, 8] = ifelse(Xc[,8] == -1, -1, round(3*GCS[,1] + 5*GCS[,2] + 4*GCS[,3] - 3)/3);
    with(stanmodel$.META, {
      list(
        N_all = nrow(data_19EI),
        unsure_spc = 0,
        nB = 6,
        B = c(3,4,5,6,7,8),# as.array(7), #vector(),
        quad_RE = quad_RE,
        nFA = pca_n_FA,
        nA_pos = length(pos_a),
        nA_neg = length(neg_a),
        A_pos = as.integer(pos_a),
        A_neg = as.integer(neg_a),
        nXc = ncol(Xc),
        nXd = ncol(Xd),
        nTd = ncol(Td),
        nTc = ncol(Tc),
        nD  = (include_d) * ncol(D),
        D_all = if (include_d) D else array(dim=c(nrow(data_19EI), 0)),
        Y_Smear_all = data_19EI$csf_smear,
        Y_Mgit_all = data_19EI$csf_mgit,
        Y_Xpert_all = data_19EI$csf_xpert,
        obs_Smear_all = data_19EI$obs_smear,
        obs_Mgit_all = data_19EI$obs_mgit,
        obs_Xpert_all = data_19EI$obs_xpert,
        Xc_all = Xc,
        Xd_all = Xd,
        Td_all = Td,
        Tc_all = Tc,
        obs_Xc_all = cbind(obs_Xc),
        obs_Xd_all = obs_Xd,
        obs_Td_all = obs_Td,
        obs_Tc_all = obs_Tc,
        nQ = length(quad_Xc),
        quad_Xc_idx = as.array(as.integer(quad_Xc)),
        penalty_family = 0,
        penalty_term = c(0, 0),
        keptin = rep(1, nrow(data_19EI))
      )
    })
  }




cli::cli_alert_info('Generate replicated data')
options(future.globals.maxSize = 1610612736*3)
rep_data <- lapply(cli::cli_progress_along(1:100), function(.) simulate_data(pretrained))
saveRDS(rep_data, 'export/rep_data.RDS')
dir.create('.cache/sampling/ppc_check', showWarnings = FALSE)
cli::cli_alert_info('Refit model on replicated data')
sample <- 
  misc$stan_kfold(sampler = model,
                  list_of_datas=rep_data,
                  backend = "rstan",
                  chains = 1, cores = 15, 
                  thin = 10, 
                  merge = FALSE,
                  control = list(adapt_delta=.70, max_treedepth=12), 
                  init_r = 1, seed = 123213,
                  sample_dir = '.cache/sampling/ppc_check',
                  pars = c('p_Smear', 'p_Mgit', 'p_Xpert', 'theta', 'a0', 'a', 'b_HIV', 'b'),
                  iter=10000, warmup=2000)
saveRDS(sample, file = 'outputs/ppc_check.RDS')