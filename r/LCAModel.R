LCAModel <- R6::R6Class(
  "LCAModel",
  public = list(
    model_name = character(),
    model = NULL,
    model_cmb = NULL,
    recipe = NULL,
    # params = NULL,
    n_fold = integer(),
    n_rep = integer(),
    folds = integer(),
    # holdout = integer(),
    initialize = function(model_file, model_dir = "outputs",
                          model_name = model_file,
                          recipe = 'data/cleaned/data_input.Rdata'
    ){
      model_file <- file.path(model_dir, paste0(model_file, ".RDS"))
      model <- readRDS(model_file)
      self$recipe <- new.env()
      load(recipe, envir = self$recipe)
      self$folds <- model$folds
      # self$holdout <- self$folds$holdout
      self$n_fold <- model$.META$fold
      self$n_rep <- model$.META$rep
      self$model_name <- model_name
      self$model <-  model$outputs 
      self$model_cmb <- rstan::sflist2stanfit(self$model)
      # self$params <- model$.META$params
      private$.META <- model$.META
      self$update_misc()
      invisible(model)
    },
    print = function(){
      cat("TBM LCA model", self$model_name, "\n")
      print(rstan::summary(self$model_cmb)$summary)
      invisible(self)
    },
    summary = function(...){
      rstan::summary(self$model_cmb, ...)
    },
    update_misc = function(){
      source('r/include/functions.R', local=private$.misc)
    },
    extract = function(...){
      rstan::extract(self$model_cmb, ...)
    },
    extract_heldout = function(pars, bind = TRUE){
      extracted_data <- vector('list', len = self$n_rep)
      for (n in seq_len(self$n_rep)){
        m <- private$.misc$extract_K_fold(self$model[((n-1)*self$n_fold+1):(n*self$n_fold)], self$folds$holdout[((n-1)*self$n_fold+1):(n*self$n_fold)], pars = pars)
        extracted_data[[n]] <- m
      }
      extracted_data <- purrr::transpose(extracted_data)
      if (bind) {
        extracted_data <- sapply(extracted_data,
                                 function(data){
                                   do.call(rbind, data)
                                 }, simplify = FALSE, USE.NAMES = TRUE)
      }
      extracted_data
    },
    get_metrics = function(file = NULL){
      s <- list(
        model_name = self$model_name,
        p = self$p,
        p_rep = self$p_rep,
        elpd = self$elpd
      )
      if (!length(file)) return(s) 
      saveRDS(s, file)
    },
    export_to_Shiny = function(file = NULL){
      pars <- self$meta$params
      pars <- grep('^[ab]', pars, value = TRUE)
      extracted <- rstan::extract(self$model_cmb, pars=pars)
      model_name <- self$model_name
      obj <- list(coef = extracted, model_mode = 'full', model_name = model_name)
      if (length(file)) return(
        saveRDS(
          obj,
          file = file
        )
      )
      obj
    },
    export_to_simplified = function(file = NULL, .full = FALSE){
      ztheta <- rstan::extract(self$model_cmb,
                          pars = 'theta')$theta %>% qlogis()
      export <- list(
        mu_ztheta = apply(ztheta, 2, mean),
        sd_ztheta = apply(ztheta, 2, sd),
        if (.full) ztheta = ztheta,
        recipe = self$recipe,
        model_name = self$model_name
      )
      
      class(export) <- 'TBMLCA_export'
      if (length(file)) return(saveRDS(export, file))
      export
    },
    mcmc_intervals = rlang::new_function(
      args = append(rlang::fn_fmls(bayesplot::mcmc_intervals)[-1], list(ylabs = NULL)),
      body = rlang::expr({
        args <- as.list(match.call())[-1]
        args <- c(x=self$model_cmb, args[names(args)!="ylabs"])
        pl <- do.call(bayesplot::mcmc_intervals, args)
        if (!length(ylabs)) return(pl)
        pl <- suppressWarnings(private$.misc$change_ylabs(pl, labs = ylabs))
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
      method = c("loess", "splines", "gam"),
      C = self$recipe$data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
      span = .75,
      knots = 3,
      est = c("mean", "median"),
      plot_rep = FALSE, theme = ggplot2::theme_bw(),
      ...
    ){
      which <- match.arg(which);
      method <- match.arg(method)
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      
      p_summary <- self$p
      p_rep <- if (self$n_rep == 1 || !plot_rep) NULL else self$p_rep
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        private$.misc$calib_curve(p_summary$p_Smear[[est]], lapply(p_rep, function(.) .$p_Smear[[est]]), Y_Smear_all, "ZN-Smear", span=span, method = method, knots=knots,..., theme = theme) + xlab("") + 
          private$.misc$calib_curve(p_summary$p_Mgit[[est]], lapply(p_rep, function(.) .$p_Mgit[[est]]), Y_Mgit_all,   "MGIT",  span=span,method = method, knots=knots,..., theme = theme) + ylab("") +
          private$.misc$calib_curve(p_summary$p_Xpert[[est]], lapply(p_rep, function(.) .$p_Xpert[[est]]),Y_Xpert_all, "Xpert", span=span,method = method, knots=knots,..., theme = theme) + xlab("") + ylab("")
      
      } else {
        not.na <- which(!is.na(C))
        C <- na.omit(C)
        private$.misc$calib_curve(p_summary$theta[[est]][not.na], lapply(p_rep, function(.) .$theta[[est]][not.na]), C, "TBM diagnosis", span=span, method = method, knots=knots,..., theme = theme)
      }
    },
    density_plot = function(
      which = c("Y", "C"),
      C = self$recipe$data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
      est = c("mean", "median"),
      theme = ggplot2::theme_bw()
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
        
        (classifierplots::density_plot(as.integer(Y_Smear_all), p_summary$p_Smear[[est]]) + ggplot2::ggtitle("ZN-Smear") + scale_x_continuous(name="", trans='reverse') + theme) +
          (classifierplots::density_plot(as.integer(Y_Mgit_all), p_summary$p_Mgit[[est]]) + ggplot2::ggtitle("MGIT")  + scale_y_continuous(name="", expand=expansion(0))+ theme) + 
          (classifierplots::density_plot(as.integer(Y_Xpert_all), p_summary$p_Xpert[[est]]) + ggplot2::ggtitle("Xpert")  + scale_x_continuous(name="", trans='reverse') + scale_y_continuous(name="", expand=expansion(0)) + theme)+ 
          plot_layout(guides='collect') & ggplot2::theme(legend.position = 'bottom')
      } else {
        not.na <- which(!is.na(C))
        C <- na.omit(C)
        classifierplots::density_plot(as.integer(C), p_summary$theta[[est]][not.na]) + theme + ggplot2::ggtitle("TBM diagnosis")
      }
    },
    roc_plot = function(
      which = c("Y", "C", "simulated_C", "smear", "mgit", "xpert"),
      C = self$recipe$data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
      est = c("mean", "median"),
      plot_rep = FALSE,
      j_index = TRUE,
      cutoffs.at = seq(0,1,.1),
      theme = ggplot2::theme_bw(),
      ...
    ){
      which <- match.arg(which)
      est <- match.arg(est)
      require(ggplot2)
      require(patchwork)
      
      my_ggroc <- private$.misc$my_ggroc
      
      p_summary <- self$p
      p_rep <- if (self$n_rep == 1 || !plot_rep) NULL else self$p_rep
      get_pred <- function(x){
        lapply(p_rep, function(.) .[[x]][[est]])
      }
      
      # if plotting y, do faceting -------------
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        pred_reps = list(
          force(get_pred('p_Smear')), 
          force(get_pred('p_Mgit')),
          force(get_pred('p_Xpert'))
        )
        patch <-
          my_ggroc(
            Y_Smear_all, 
            p_summary$p_Smear[[est]], 
            pred_reps[[1]],
            j_index,
            cutoffs.at, ...
          ) + theme + ggplot2::ggtitle("ZN-Smear") + ggplot2::theme(axis.title = element_blank()) +
          my_ggroc(
            Y_Mgit_all, 
            p_summary$p_Mgit[[est]], 
            pred_reps[[2]],
            j_index,
            cutoffs.at, ...
          ) + theme + ggplot2::ggtitle("MGIT") + ggplot2::theme(axis.title = element_blank()) +
          my_ggroc(
            Y_Xpert_all, 
            p_summary$p_Xpert[[est]], 
            pred_reps[[3]],
            j_index,
            cutoffs.at, ...
          ) + theme + ggplot2::ggtitle("Xpert") + ggplot2::theme(axis.title = element_blank())
        
        gt <- patchwork::patchworkGrob(patch)
        x <- list(gt, left = "True Positive Rate (%)    (Sensitivity)", bottom = "False Positive Rate (%)    (1-Specificity)")
        class(x) <- 'multggplotGrob'
        return(x)
      }
      
      if (which == "smear") {
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        pred_rep = get_pred('p_Smear')
        x <- my_ggroc(
          Y_Smear_all, 
          p_summary$p_Smear[[est]], 
          pred_rep,
          j_index,
          cutoffs.at, ...
        ) + theme + ggplot2::ggtitle("ZN-Smear")
        return(x)
      }
      
      if (which == "mgit") {
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        pred_rep = get_pred('p_Mgit')
        x <- my_ggroc(
          Y_Mgit_all, 
          p_summary$p_Mgit[[est]], 
          pred_rep,
          j_index,
          cutoffs.at, ...
        ) + theme + ggplot2::ggtitle("MGIT")  
        return(x)
      } 
      
      if (which == "xpert") {
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        pred_rep = get_pred('p_Xpert')
        x <- my_ggroc(
          Y_Xpert_all, 
          p_summary$p_Xpert[[est]], 
          pred_rep,
          j_index,
          cutoffs.at, ...
        ) + theme + ggplot2::ggtitle("Xpert")
        
        return(x)
      } 
      
      # if C ----------
      not.na <- which(!is.na(C))
      C <- na.omit(C)
      pred_rep = get_pred('theta') |> lapply(function(x) x[not.na])
      plt <- my_ggroc(
        C, 
        p_summary$theta[[est]][not.na], 
        pred_rep,
        j_index,
        cutoffs.at, ...) + theme + ggplot2::ggtitle("TBM diagnosis")
      
      plt
    },
    roc_plot_legacy = function(
      which = c("Y", "C", "simulated_C", "smear", "mgit", "xpert"),
      C = self$recipe$data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
      resamps = 2000,
      force_bootstrap = NULL,
      est = c("mean", "median"),
      plot_rep = FALSE,
      theme = ggplot2::theme_bw(),
      ...
    ){
      which <- match.arg(which)
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      
      p_summary <- self$p
      p_rep <- if (self$n_rep == 1 || !plot_rep) NULL else self$p_rep
      
      
      
      if (which == "Y"){
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        
        suppressMessages(
          patch <-
            private$.misc$my_roc_plot(Y_Smear_all, p_summary$p_Smear[[est]], lapply(p_rep, function(.) .$p_Smear[[est]]), resamps = resamps, force_bootstrap = force_bootstrap, ...) + theme + ggplot2::ggtitle("ZN-Smear") + ggplot2::theme(axis.title = element_blank()) + # + scale_x_continuous(name="") +
            private$.misc$my_roc_plot(Y_Mgit_all , p_summary$p_Mgit[[est]] ,lapply(p_rep, function(.) .$p_Mgit[[est]]), resamps = resamps, force_bootstrap = force_bootstrap, ...) + theme + ggplot2::ggtitle("MGIT") + ggplot2::theme(axis.title = element_blank()) + #+ scale_y_continuous(name="") +  
            private$.misc$my_roc_plot(Y_Xpert_all, p_summary$p_Xpert[[est]],lapply(p_rep, function(.) .$p_Xpert[[est]]), resamps = resamps, force_bootstrap = force_bootstrap, ...) + theme + ggplot2::ggtitle("Xpert") + ggplot2::theme(axis.title = element_blank())#+ scale_x_continuous(name="") + scale_y_continuous(name="")
        )
        
        gt <- patchwork::patchworkGrob(patch)
        x <- list(gt, left = "True Positive Rate (%)    (Sensitivity)", bottom = "False Positive Rate (%)    (1-Specificity)")
        class(x) <- 'multggplotGrob'
        x
      } else if (which == "smear") {
        Y_Smear_all <- self$folds$inputs[[1]]$Y_Smear_all
        private$.misc$my_roc_plot(Y_Smear_all, p_summary$p_Smear[[est]], lapply(p_rep, function(.) .$p_Smear[[est]]), resamps = resamps, force_bootstrap = force_bootstrap, ...) + theme + ggplot2::ggtitle("ZN-Smear")
      } else if (which == "mgit") {
        Y_Mgit_all  <- self$folds$inputs[[1]]$Y_Mgit_all
        private$.misc$my_roc_plot(Y_Mgit_all , p_summary$p_Mgit[[est]] ,lapply(p_rep, function(.) .$p_Mgit[[est]]), resamps = resamps, force_bootstrap = force_bootstrap,...) + theme + ggplot2::ggtitle("MGIT")  
      } else if (which == "xpert") {
        Y_Xpert_all <- self$folds$inputs[[1]]$Y_Xpert_all
        private$.misc$my_roc_plot(Y_Xpert_all, p_summary$p_Xpert[[est]],lapply(p_rep, function(.) .$p_Xpert[[est]]), resamps = resamps, force_bootstrap = force_bootstrap, ...) + theme + ggplot2::ggtitle("Xpert")
      } else {
        not.na <- which(!is.na(C))
        C <- na.omit(C)
        plt <- private$.misc$my_roc_plot(C, p_summary$theta[[est]][not.na], lapply(p_rep, function(.) .$theta[[est]][not.na]), resamps = resamps, force_bootstrap = force_bootstrap, ...) + theme + ggplot2::ggtitle("TBM diagnosis")
        # ggplot2::ggplotGrob(plt)
      }  
    }
  ),
  private = list(.misc = new.env(), .elpd = NULL, .loglik = NULL, .p = NULL, .p_rep = NULL, .META = list()),
  active = list(
    loglik = function(){
      if (length(private$.loglik)) return(private$.loglik)
      private$.loglik <- private$.misc$extract_log_lik_K(self$model, self$folds$holdout)
    },
    elpd =function(){
      if (length(private$.elpd)) return(private$.elpd)
      loglik <- self$loglik
      private$.elpd <- private$.misc$kfold(loglik)
    },
    p = function(){
      if (is.null(private$.p)){
        p <- private$.misc$extract_K_fold(self$model, self$folds$holdout, pars = c("z_theta", "theta", "p_Smear", "p_Mgit", "p_Xpert"))
        p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
        p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)
        private$.p <- p_summary
      } 
      
      private$.p
    },
    p_rep = function(){
      if (self$n_rep == 1) return(self$p)
      if (is.null(private$.p_rep)){
        private$.p_rep <- vector("list", self$n_rep)
        for (n in seq_len(self$n_rep)){
          p <- private$.misc$extract_K_fold(self$model[((n-1)*self$n_fold+1):(n*self$n_fold)], self$folds$holdout[((n-1)*self$n_fold+1):(n*self$n_fold)], pars = c("z_theta","theta", "p_Smear", "p_Mgit", "p_Xpert"))
          p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
          p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)
          private$.p_rep[[n]] <- p_summary
        }
      }
      private$.p_rep
    },
    meta = function(){
      private$.META
    }
  )
)

print.multggplotGrob <- plot.multggplotGrob <- function(x, ...){
  do.call(gridExtra::grid.arrange, x)
}

print.ggplotGrob <- plot.ggplotGrob <- function(x, ...){
  plot(x, ...)
}
