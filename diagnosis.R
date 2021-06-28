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
      cat("TBM LCA model", self$model_name, "\n")
      print(rstan::summary(self$model_cmb)$summary)
      invisible(self)
    },
    summary = function(...){
      rstan::summary(self$model_cmb, ...)
    },
    mcmc_intervals = rlang::new_function(
      args = append(rlang::fn_fmls(bayesplot::mcmc_intervals)[-1], list(ylabs = NULL)),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args[names(args)!="ylabs"])
        pl <- do.call(bayesplot::mcmc_intervals, args)
        if (!length(args$ylabs)) return(pl)
        pl <- suppressWarnings(private$.misc$change_ylabs(pl, labs = args$ylabs))
        pl
      })
    ),
    mcmc_dens_chains = rlang::new_function(
      args = rlang::fn_fmls(bayesplot::mcmc_dens_chains)[-1],
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args)
        do.call(bayesplot::mcmc_dens_chains, args)
      })
    ),
    mcmc_dens = rlang::new_function(
      args = c(rlang::fn_fmls(bayesplot::mcmc_dens)[-1]),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args)
        do.call(bayesplot::mcmc_dens, args)
      })
    ),
    mcmc_dens_overlay = rlang::new_function(
      args = c(rlang::fn_fmls(bayesplot::mcmc_dens_overlay)[-1]),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args)
        do.call(bayesplot::mcmc_dens_overlay, args)
      })
    ),
    mcmc_plot = function(fn, ..., eval = TRUE, pkg = 'bayesplot'){
      fn = getFromNamespace(fn, pkg)
      if (eval) return(fn(x = self$model_cmb, ...))
      rlang::new_function(
        args = c(formals(fn)[-1]),
        body = rlang::expr({
          args <- as.list(match.call())[-1]
          args <- c(x = self$model_cmb, args)
          do.call(fn, args)
        })
      )
    },
    mcmc_areas = rlang::new_function(
      args = append(rlang::fn_fmls(bayesplot::mcmc_areas)[-1], list(ylabs = NULL)),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x = self$model_cmb, args)
        pl <- do.call(bayesplot::mcmc_areas, args)
        if (!length(args$ylabs)) return(pl)
        pl <- suppressWarnings(private$.misc$change_ylabs(pl, labs = args$ylabs))
        pl
      })
    ),
    calibration_plot = function(
      which = c("Y", "C"),
      method = c("loess", "splines"),
      C = !self$recipe$data_19EI$other_dis_dx,
      span = .75,
      knots = 5,
      est = c("mean", "median"),
      ...
    ){
      which <- match.arg(which); method <- match.arg(method)
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      
      p_summary <- self$p
      
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        private$.misc$calib_curve(p_summary$p_Smear[[est]],Y_Smear_all, "Smear", span=span, method = method, knots=knots,...) +
          private$.misc$calib_curve(p_summary$p_Mgit[[est]],Y_Mgit_all,   "Mgit",  span=span,method = method, knots=knots,...) +
          private$.misc$calib_curve(p_summary$p_Xpert[[est]],Y_Xpert_all, "Xpert", span=span,method = method, knots=knots,...) 
        # if (method == "loess"){
        #   private$.misc$calib_curve(p_summary$p_Smear[[est]],Y_Smear_all, "Smear", span=span, ...) +
        #   private$.misc$calib_curve(p_summary$p_Mgit[[est]],Y_Mgit_all,   "Mgit",  span=span, ...) +
        #   private$.misc$calib_curve(p_summary$p_Xpert[[est]],Y_Xpert_all, "Xpert", span=span, ...) 
        # } else {
        #   ggplot(mapping=aes(x=p_summary$p_Smear[[est]], y=Y_Smear_all)) + 
        #   stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red", ...) + 
        #   geom_line(aes(y=p_summary$p_Smear[[est]])) +
        #   ggplot(mapping=aes(x=p_summary$p_Mgit[[est]], y=Y_Mgit_all)) + 
        #   stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red",...) + 
        #   geom_line(aes(y=p_summary$p_Mgit[[est]])) + 
        #   ggplot(mapping=aes(x=p_summary$p_Xpert[[est]], y=Y_Xpert_all)) + 
        #   stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red", ...) + 
        #   geom_line(aes(y=p_summary$p_Xpert[[est]]))
        # }
      } else {
        private$.misc$calib_curve(p_summary$theta[[est]], C, "Positive TBM", span=span, method = method, knots=knots,...)
        # if (method == "loess")
        #   private$.misc$calib_curve(p_summary$theta[[est]], C, "Positive TBM", span=span)
        # else 
        #   ggplot(mapping=aes(x=p_summary$theta[[est]], y=as.numeric(C))) + 
        #     stat_smooth(method="glm", formula=y~splines::ns(x,knots), color = "red", ...) + 
        #     geom_line(aes(y=p_summary$theta[[est]]))
      }
    },
    density_plot = function(
      which = c("Y", "C"),
      C = !self$recipe$data_19EI$other_dis_dx,
      est = c("mean", "median")
    ){
      which <- match.arg(which)
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      p_summary <- self$p
      
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        classifierplots::density_plot(Y_Smear_all, p_summary$p_Smear[[est]]) + ggplot2::ggtitle("Smear") + 
          classifierplots::density_plot(Y_Mgit_all , p_summary$p_Mgit[[est]]) + ggplot2::ggtitle("Mgit") + 
          classifierplots::density_plot(Y_Xpert_all, p_summary$p_Xpert[[est]]) + ggplot2::ggtitle("Xpert")
      } else 
        classifierplots::density_plot(as.integer(C), p_summary$theta[[est]]) + ggplot2::ggtitle("Positive TBM")
    },
    roc_plot = function(
      which = c("Y", "C"),
      C = !self$recipe$data_19EI$other_dis_dx,
      resamps = 2000,
      force_bootstrap = NULL,
      est = c("mean", "median")
    ){
      which <- match.arg(which)
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      p_summary <- self$p
      
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        classifierplots::roc_plot(Y_Smear_all, p_summary$p_Smear[[est]], resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Smear") + 
        classifierplots::roc_plot(Y_Mgit_all , p_summary$p_Mgit[[est]] , resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Mgit") + 
        classifierplots::roc_plot(Y_Xpert_all, p_summary$p_Xpert[[est]], resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Xpert")
      } else 
        classifierplots::roc_plot(C, p_summary$theta[[est]], resamps = resamps, force_bootstrap = force_bootstrap) + ggplot2::ggtitle("Positive TBM")
    }
  ),
  private = list(.misc = new.env(), .elpd = NULL, .loglik = NULL, .p = NULL),
  active = list(
    loglik = function(){
      if (length(private$.loglik)) return(private$.loglik)
      private$.loglik <- private$.misc$extract_log_lik_K(self$model, self$holdout)
    },
    elpd =function(){
      if (length(private$.elpd)) return(private$.elpd)
      loglik <- self$loglik
      private$.elpd <- private$.misc$kfold(loglik)
    },
    p = function(){
      if (is.null(private$.p)){
        p <- private$.misc$extract_K_fold(self$model, self$holdout, pars = c("theta", "p_Smear", "p_Mgit", "p_Xpert"))
        p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
        p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)
        private$.p <- p_summary
      } 
      
      private$.p
    }
  )
)