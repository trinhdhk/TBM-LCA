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
    initialize = function(...){
      model = super$initialize(...)
      private$.theta = model$theta
      invisible(model)
    },
    print = function(){
      cat("TBM sLCA model", self$model_name, "\n")
      print(rstan::summary(self$model_cmb)$summary)
      invisible(self)
    },
    calibration_plot = function(
      which = c("theta", "lambda", "C", "simulated_C"),
      method = c("loess", "splines", "gam"),
      C = self$recipe$data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
      span = .75,
      knots = 3,
      est = c("mean", "median"),
      plot_rep = FALSE,
      theme = ggplot2::theme_bw(),
      ...
    ){
      which <- match.arg(which); method <- match.arg(method)
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      
      p_summary <- self$p
      p_rep <- if (self$n_rep == 1 || !plot_rep) NULL else self$p_rep
      if (which == "theta"){
        Y <- apply(private$.theta, 2, mean)
        private$.misc$calib_curve(p_summary$theta[[est]], lapply(p_rep, function(.) (.$theta[[est]])), Y, "TBM probability", type='linear', span=span, method = method, knots=knots,..., theme = theme)
      } else if (which == "lambda"){
        Y <- apply(qlogis(private$.theta), 2, mean)
        private$.misc$calib_curve(p_summary$lambda[[est]], lapply(p_rep, function(.) .$lambda[[est]]), Y, "TBM log_odd", type='linear', span=span, method = method, knots=knots,..., theme = theme)
      }
      else if (which == 'C'){
        not.na <- (!is.na(C))
        C <- na.omit(C)
        private$.misc$calib_curve(p_summary$theta[[est]][not.na], lapply(p_rep, function(.) (.$theta[[est]][not.na])), C, "Positive TBM", span=span, method = method, knots=knots,..., theme = theme)
      }
      else if (which == 'simulated_C'){
        C = self$sample_C
        private$.misc$calib_curve(p_summary$theta[[est]], lapply(p_rep, function(.) (.$theta[[est]])), C, "Positive TBM", span=span, method = method, knots=knots,..., theme = theme)
      }
    },
    density_plot = function(
      C = self$recipe$data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
      est = c("mean", "median"),
      theme = ggplot2::theme_bw() 
    ){
      require(ggplot2)
      require(patchwork)
      est <- match.arg(est)
      p_summary <- self$p
      not.na <- which(!is.na(C))
      C <- na.omit(C)
      classifierplots::density_plot(as.integer(C), (p_summary$theta[[est]][not.na])) + theme() + ggplot2::ggtitle("Positive TBM")
    },
    roc_plot = function(
      which = c("C", "simulated_C"),
      C = self$recipe$data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
      resamps = 2000,
      force_bootstrap = NULL,
      est = c("mean", "median"),
      plot_rep = FALSE,
      theme = ggplot2::theme_bw()
    ){
      require(ggplot2)
      require(patchwork)
      require(data.table)
      est <- match.arg(est)
      which <- match.arg(which)
      
      p_summary <- self$p
      p_rep <- if (self$n_rep == 1 || !plot_rep) NULL else self$p_rep
      
      if (which=="C"){
        not.na <- (!is.na(C))
        C <- na.omit(C)
        
        private$.misc$my_roc_plot(obs=C, pred=p_summary$theta[[est]][not.na], lapply(p_rep, function(.) (.$theta[[est]][not.na])), resamps = resamps, force_bootstrap = force_bootstrap) + theme + ggplot2::ggtitle("Positive TBM")
      } else {
        plt <- ggplot()
        auc <- numeric(resamps)
        sen <- spc <- threshold <- matrix(nrow = 201, ncol = resamps)
        bounds <- c(NA, NA)
        for (i in seq_len(resamps)){
          C <- self$sample_C
          j <- sample(seq_along(self$p_rep), 1)
          plt <- private$.misc$._add_roc_plot(plt, C, self$p_rep[[j]]$theta[[est]], resamps, force_bootstrap, .rep = TRUE, .has_rep = TRUE, bounds=bounds)
          roc <- pROC::roc(C ~ self$p_rep[[j]]$theta[[est]])
          auc[i] <- roc$auc
          # threshold[[i]] <- roc$threshold
          # sen[[i]] <- roc$sensitivities
          # spc[[i]] <- roc$specificities
          coord.mat <- pROC::coords(roc, seq(0, 1, .005))
          sen[,i] <- coord.mat[,'sensitivity']
          spc[,i] <- coord.mat[,'specificity']
          threshold[,i] <- coord.mat[,'threshold']
        }
        sen.mean <- apply(sen, 1, mean)
        spc.mean <- apply(spc, 1, mean)
        threshold.mean <- apply(threshold, 1, mean)
        cutoff <- which.max(sen.mean + spc.mean)
        auc.mean <- mean(auc)
        auc.ci <- quantile(auc, c(0.25, .975))
        break.x <- c(0, 25, 50, 75, 100)
        far.break.x <- abs((round(100 - spc.mean[cutoff]*100,1)) - break.x) > 5
        break.x <- c(break.x[far.break.x], (round(100 - spc.mean[cutoff]*100,1)))
        break.y <- c(0, 25, 50, 75, 100)
        far.break.y <- abs(round(100 - spc.mean[cutoff]*100,1) - break.y) > 5
        break.y <- c(break.y[far.break.y], round(100 - spc.mean[cutoff]*100,1))
    
        plt <- plt + 
        annotate("text", x = 62.5, y = 22.5, label = paste0("AUC ", formatC(auc.mean*100, digits = 1, format="f"), "%"),
                 parse = F, size = 5, colour = classifierplots:::fontgrey_str) + 
          annotate("text", x = 62.5, y = 16.5, label = paste0("(95% CI ", 
                                                              formatC(auc.ci[1]*100, digits = 1, format="f"), '% - ', formatC(auc.ci[2]*100, digits = 1, format="f"), '%)'),
                   parse = F, size = 3.5, colour = classifierplots:::fontgrey_str) + 
          scale_x_continuous(name = "False Positive Rate (%)    (1-Specificity)", 
                             limits = c(0, 100), expand = c(0.05, 0.05),
                             breaks = break.x
          ) + 
          scale_y_continuous(name = "True Positive Rate (%)    (Sensitivity)", 
                             limits = c(0, 100), expand = c(0.05, 0.05),
                             breaks = break.y) + 
          geom_hline(aes(yintercept=sen.mean[cutoff]*100), colour = classifierplots:::fontgrey_str, linetype = "dashed") +
          geom_vline(aes(xintercept=100-spc.mean[cutoff]*100), colour = classifierplots:::fontgrey_str, linetype = "dashed") +
          annotate("text", y=sen.mean[cutoff]*100-3, x=100-spc.mean[cutoff]*100+3,
                   hjust='left', vjust='top',
                   label = glue::glue("Pr(Y) = {round(threshold.mean[cutoff]*100,1)}"),
                   size=3.8)+
          geom_line(mapping=aes(x = 100 * (1-spc.mean), 
                                y = 100 * sen.mean),
                    color =  classifierplots:::green_str, 
                    size = 1,
                    alpha =  1) 
        plt + theme
      }
    }
  ),
  private = list(.misc = new.env(), .elpd = NULL, .loglik = NULL, .p = NULL, .p_rep = NULL, .META = list(), .theta = numeric(), .mean_theta = numeric()),
  active = list(
    theta = function(){
      private$.theta
    },
    mean_theta = function(){
      mean_theta = private$.mean_theta
      if (length(mean_theta)) return(mean_theta)
      private$.mean_theta = apply(private$.theta, 2, mean)
    },
    sample_C = function(){
      mean_theta = self$mean_theta
      sapply(mean_theta, rbinom, n=1, size=1)
    },
    loglik = function(){
      if (length(private$.loglik)) return(private$.loglik)
      private$.loglik <- private$.misc$extract_log_lik_K(self$model, self$folds$holdout)
    },
    elpd = function(){
      if (length(private$.elpd)) return(private$.elpd)
      loglik <- self$loglik
      private$.elpd <- private$.misc$kfold(loglik)
    },
    p = function(){
      if (is.null(private$.p)){
        p <- private$.misc$extract_K_fold(self$model, self$folds$holdout, pars = c("lambda", "theta"))
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
          p <- private$.misc$extract_K_fold(self$model[((n-1)*self$n_fold+1):(n*self$n_fold)], self$folds$holdout[((n-1)*self$n_fold+1):(n*self$n_fold)], pars = c("lambda", "theta"))
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