LCAModel <- R6::R6Class(
  "LCAModel",
  public = list(
    model_name = character(),
    model = NULL,
    model_cmb = NULL,
    recipe = NULL,
    n_fold = integer(),
    n_rep = integer(),
    folds = integer(),
    holdout = integer(),
    p = NULL,
    initialize = function(model_name, model_dir = "outputs",
                          model_file = file.path(model_dir, paste0(model_name, ".RDS"))
    ){
      model <- readRDS(model_file)
      self$recipe <- new.env()
      load('data/cleaned/data_input.Rdata', envir = self$recipe)
      self$folds <- model$folds
      self$holdout <- self$folds$holdout
      self$n_fold <- model$fold
      self$n_rep <- model$n_rep
      self$model_name <- model_name
      self$model <- model$outputs
      self$model_cmb <- rstan::sflist2stanfit(self$model)
      source('r/include/functions.R', local=private$.misc)
    },
    print = function(){
      cat("TBM LCA model", self$model_name, "\\n")
      print(rstan::summary(self$model_cmb)$summary)
      invisible(self)
    },
    summary = function(...){
      rstan::summary(self$model_cmb, ...)
    },
    mcmc_intervals = rlang::new_function(
      args = c(rlang::fn_fmls(bayesplot::mcmc_intervals)[-1], ylabs = NULL),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args[names(args)!="ylabs"])
        pl <- do.call(bayesplot::mcmc_intervals, args)
        if (!length(args$ylabs)) return(pl)
        suppressWarnings(pl <- private$.misc$change_ylabs(pl, labs = args$ylabs))
        pl
      })
    ),
    mcmc_dens_chains = rlang::new_function(
      args = c(rlang::fn_fmls(bayesplot::mcmc_dens_chains)[-1]),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args)
        do.call(bayesplot::mcmc_dens_chains, args)
      })
    ),
    mcmc_dens = rlang::new_function(
      args = c(rlang::fn_fmls(bayesplot::mcmc_dens_chains)[-1]),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args)
        do.call(bayesplot::mcmc_dens_chains, args)
      })
    ),
    mcmc_areas = rlang::new_function(
      args = c(rlang::fn_fmls(bayesplot::mcmc_areas)[-1]),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x = self$model_cmb, args)
        do.call(bayesplot::mcmc_areas, args)
      })
    ),
    calibration_plot = function(
      which = c("Y", "C"),
      method = c("loess", "splines"),
      C = !self$recipe$data_19EI$other_dis_dx,
      span = .75,
      knots = 5,
      ...
    ){
      which <- match.arg(which); method <- match.arg(method)
      require(ggplot2)
      require(patchwork)
      
      if (is.null(self$p)){
        p <- private$.misc$extract_K_fold(self$model, self$holdout, pars = c("theta", "p_Smear", "p_Mgit", "p_Xpert"))
        p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
        p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)
        self$p <- p_summary
      } else p_summary <- self$p
      
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        if (method == "loess"){
          private$.misc$calib_curve(p_summary$p_Smear$mean,Y_Smear_all, "Smear", span=span, ...) +
          private$.misc$calib_curve(p_summary$p_Mgit$mean,Y_Mgit_all,   "Mgit",  span=span, ...) +
          private$.misc$calib_curve(p_summary$p_Xpert$mean,Y_Xpert_all, "Xpert", span=span, ...) 
        } else {
          ggplot(mapping=aes(x=p_summary$p_Smear$mean, y=Y_Smear_all)) + 
          stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red", ...) + 
          geom_line(aes(y=p_summary$p_Smear$mean)) +
          ggplot(mapping=aes(x=p_summary$p_Mgit$mean, y=Y_Mgit_all)) + 
          stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red",...) + 
          geom_line(aes(y=p_summary$p_Mgit$mean)) + 
          ggplot(mapping=aes(x=p_summary$p_Xpert$mean, y=Y_Xpert_all)) + 
          stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red", ...) + 
          geom_line(aes(y=p_summary$p_Xpert$mean))
        }
      } else {
        
        if (method == "loess")
          private$.misc$calib_curve(p_summary$theta$mean, C, "Positive TBM", span=span)
        else 
          ggplot(mapping=aes(x=p_summary$theta$mean, y=as.numeric(C))) + 
            stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red", ...) + 
            geom_line(aes(y=p_summary$theta$mean))
      }
    },
    roc_plot = function(
      which = c("Y", "C"),
      C = !self$recipe$data_19EI$other_dis_dx,
      resamps = 2000,
      force_bootstrap = NULL
    ){
      which <- match.arg(which)
      require(ggplot2)
      require(patchwork)
      
      if (is.null(self$p)){
        p <- private$.misc$extract_K_fold(self$model, self$holdout, pars = c("theta", "p_Smear", "p_Mgit", "p_Xpert"))
        p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
        p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)
        self$p <- p_summary
      } else p_summary <- self$p
      
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        classifierplots::roc_plot(Y_Smear_all, p_summary$p_Smear$mean, resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Smear") + 
        classifierplots::roc_plot(Y_Mgit_all , p_summary$p_Mgit$mean , resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Mgit") + 
        classifierplots::roc_plot(Y_Xpert_all, p_summary$p_Xpert$mean, resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Xpert")
      } else 
        classifierplots::roc_plot(C, p_summary$theta$mean, resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Positive TBM")
    }
  ),
  private = list(.misc = new.env(), .elpd = NULL, .loglik = NULL),
  active = list(
    loglik = function(){
      if (length(private$.loglik)) return(private$.loglik)
      private$.loglik <- private$.misc$extract_log_lik_K(self$model, self$holdout)
    },
    elpd =function(){
      if (length(private$.elpd)) return(private$.elpd)
      loglik <- self$loglik
      private$.elpd <- private$.misc$kfold(loglik)
    }
  )
)