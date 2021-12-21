sLCAModel <- R6::R6Class(
  "sLCAModel",
  inherit = LCAModel,
  public = list(
    model_name = character(),
    model = NULL,
    model_cmb = NULL,
    recipe = NULL,
    n_fold = integer(),
    n_rep = integer(),
    folds = integer(),
    print = function(){
      cat("TBM sLCA model", self$model_name, "\n")
      print(rstan::summary(self$model_cmb)$summary)
      invisible(self)
    },
    calibration_plot = function(
      which = c("ztheta", "theta", "C"),
      method = c("loess", "splines"),
      C = !self$recipe$data_19EI$other_dis_dx,
      span = .75,
      knots = 3,
      est = c("mean", "median"),
      plot_rep = FALSE, theme = ggplot2::theme_bw,
      ...
    ){
      which <- match.arg(which); method <- match.arg(method)
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      
      p_summary <- self$p
      p_rep <- if (self$n_rep == 1 || !plot_rep) NULL else self$p_rep
      if (which == "ztheta"){
        Y <- self$folds$inputs[[1]]$mu_ztheta
        private$.misc$calib_curve(p_summary$lambda[[est]], lapply(p_rep, function(.) .$lambda[[est]]), Y, "TBM log-odd", type='linear', span=span, method = method, knots=knots,..., theme = theme)
      } else if (which == "theta"){
        Y <- self$folds$inputs[[1]]$mu_ztheta |> plogis()
        private$.misc$calib_curve(p_summary$lambda[[est]] |> plogis(), lapply(p_rep, function(.) plogis(.$lambda[[est]])), Y, "TBM probability", type='linear', span=span, method = method, knots=knots,..., theme = theme)
      } else {
        private$.misc$calib_curve(p_summary$lambda[[est]] |> plogis(), lapply(p_rep, function(.) plogis(.$lambda[[est]])), C, "Positive TBM", span=span, method = method, knots=knots,..., theme = theme)
      }
    },
    density_plot = function(
      C = !self$recipe$data_19EI$other_dis_dx,
      est = c("mean", "median"),
      theme = ggplot2::theme_bw 
    ){
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      p_summary <- self$p
      classifierplots::density_plot(as.integer(C), plogis(p_summary$lambda[[est]])) + theme() + ggplot2::ggtitle("Positive TBM")
    },
    roc_plot = function(
      C = !self$recipe$data_19EI$other_dis_dx,
      resamps = 2000,
      force_bootstrap = NULL,
      est = c("mean", "median"),
      plot_rep = FALSE,
      theme = ggplot2::theme_bw
    ){
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      
      p_summary <- self$p
      p_rep <- if (self$n_rep == 1 || !plot_rep) NULL else self$p_rep
      
      private$.misc$my_roc_plot(C, plogis(p_summary$lambda[[est]]), lapply(p_rep, function(.) plogis(.$lambda[[est]])), resamps = resamps, force_bootstrap = force_bootstrap) + theme() + ggplot2::ggtitle("Positive TBM")
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
        p <- private$.misc$extract_K_fold(self$model, self$folds$holdout, pars = "lambda")
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
          p <- private$.misc$extract_K_fold(self$model[((n-1)*self$n_fold+1):(n*self$n_fold)], self$folds$holdout[((n-1)*self$n_fold+1):(n*self$n_fold)], pars = "lambda")
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