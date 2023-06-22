StanKfold <- R6::R6Class("StanKfold",
                         public = list(
                           mc.cores = 1,
                           compile_sampler = function(stan_file = NULL, stan_obj,
                                                      compiler_args = list(), force = FALSE){
                             backend <- private$.backend
                             if (force) private$.sampler <- private$.stan_file <- NULL
                             else if (is.null(private$.sampler)){
                               if (missing(stan_obj)){
                                 compiler <- private$.sampler_compiler(stan_file, backend)   
                                 sampler  <- do.call(compiler, compiler_args)
                               } else {
                                 sampler <- stan_obj
                                 backend2 <- switch(class(sampler), 
                                                    stanmodel = "rstan", CmdStanModel = "cmdstanr")
                                 if (is.na(backend2)) stop("Unrecognised Stan object.")
                                 if (backend != backend2){
                                   cat("Backend set to ", backend2)
                                   private$.backend <- backend2
                                 }
                               }
                               
                               private$.stan_file <- stan_file
                               private$.sampler <- sampler
                             }
                             
                             invisible(self)
                           },
                           set_backend = function(backend = c("rstan", "cmdstanr")){
                             backend <- match.arg(backend)
                             private$.backend <- backend
                             invisible(self)
                           },
                           set_data = function(data){
                             if (!length(data$N)) stop("data must have an elements called N!")
                             
                             private$.raw_data <- data
                             private$.N <- data$N
                             
                             holdout_obj <- private$.create_holdout()
                             private$.keptin <- holdout_obj$keptin
                             private$.holdout <- holdout_obj$holdout
                             private$.inputs <- holdout_obj$inputs
                             cat("Succesfully created k_fold inputs.\n")
                             invisible(self)
                           },
                           sampling = function(chains, mc.cores = self$mc.cores,...){
                             sampler <- private$.sampler
                             if (!length(sampler)) stop("Sampler has not been set!")
                             backend <- private$.backend
                             seed <- private$.seed
                             K <- private$.K
                             list_of_datas <- private$.inputs
                             future::plan(future::multisession, workers=mc.cores)
                             wd <- getwd()
                             progressr::with_progress({
                               p <- progressr::progressor(steps = K*chains)
                               
                               sflist <- furrr::future_map(seq_len(K*chains), function(i){
                                 setwd(wd)
                                 k <- ceiling(i / chains)
                                 if (backend == "rstan")
                                   sf <- rstan::sampling(sampler, data = list_of_datas[[k]], chains=1, seed = seed+i, 
                                                         chain_id = i, ...)
                                 else {
                                   m <- sampler$clone()
                                   s <- m$sample(data = list_of_datas[[k]],
                                                 chains = 1, seed = seed + i, chain_ids = i,...)
                                   sf <- rstan::read_stan_csv(s$output_files())
                                 }
                                 p()
                                 sf
                               })
                             })
                             stanfit <- list()
                             for(k in 1:K){
                               inchains <- (chains*k - (chains - 1)):(chains*k)
                               #  Merge `chains` of each fold
                               stanfit[[k]] <- rstan::sflist2stanfit(sflist[inchains])
                             }  
                             StanKfold_MCMC$new(K = K, N_Rep = private$.N_Rep, inputs = private$.inputs,
                                                keptin = private$.keptin, private$.holdout,
                                                stanfits = stanfit, mc.cores = self$mc.cores)
                           },
                           
                           initialize = function(stan_file, stan_obj,
                                                 K = 10, N_Rep = 1, seed = 123,
                                                 backend,
                                                 complier_args = list(),
                                                 mc.cores = 1){
                             
                             if (!missing(backend)) self$set_backend(backend)
                             if (!missing(stan_obj)) self$compile_sampler(stan_obj = stan_obj, compiler_args = complier_args)
                             else if (!missing(stan_file)) self$compile_sampler(stan_file = stan_file, compiler_args = complier_args)
                             
                             private$.K <- K
                             private$.N_Rep <- N_Rep
                             private$.seed <- seed
                             self$mc.cores <- mc.cores
                             
                             invisible(self)
                           },
                           print = function(){
                             binded <- if (length(private$.data)) "Binded" else "Unbinded"
                             writeLines(glue::glue("{binded} Stan Kfold object with {private$.K} folds and {private$.N_Rep} repetitions.\n"))
                             writeLines(glue::glue("Using {private$.backend} as backend and {private$.seed} as random seed.\n"))
                           }
                         ),
                         private = list(
                           .backend = "rstan", .stan_file = NULL, .sampler = NULL,
                           .K = NULL, .N_Rep = NULL, .seed = NULL, 
                           .raw_data = NULL, .N = NULL, .inputs = NULL,
                           .keptin = NULL, .holdout = NULL,
                           .create_holdout = function(K = private$.K, N_rep = private$.N_Rep, N_obs = private$.N, seed = private$.seed, cores = self$mc.cores){
                             holdout <- vector("list", N_rep * K)
                             keptin <- vector("list", N_rep * K)
                             for(n in seq_len(N_rep)){
                               set.seed(seed + n - 1)
                               hh <- loo::kfold_split_random(K = K, N = N_obs)
                               foldkept <- matrix(1, nrow = N_obs, ncol = K)
                               for(i in 1:N_obs) foldkept[i, hh[i]] <- 0
                               foldhold  <- 1- foldkept
                               keptin[((n-1)*K+1):(n*K)] <- split(foldkept,rep(1:ncol(foldkept),each=nrow(foldkept)))
                               holdout[((n-1)*K+1):(n*K)]  <- split(foldhold,rep(1:ncol(foldhold),each=nrow(foldhold)))
                             }
                             inputs <- pbmcapply::pbmclapply(seq_len(K * N_rep),
                                                             mc.cores=self$mc.cores,
                                                             function(k) modifyList(private$.raw_data, list(keptin=keptin[[k]])))
                             list(keptin=keptin, holdout=holdout, inputs = inputs)
                           },
                           .sampler_compiler = function(stan_file, backend){
                             switch(backend,
                                    cmdstanr = purrr::partial(cmdstanr::cmdstan_model, stan_file=stan_file),
                                    rstan    = purrr::partial(rstan::stan_model, file = stan_file))
                           }
                         ),
                         active = list(
                           seed = function(){
                             private$.seed
                           },
                           backend = function(backend){
                             if (missing(backend)) return(private$.backend)
                             self$set_backend(backend)
                           },
                           data = function(data){
                             if (missing(data)) return(private$.raw_data)
                             self$set_data(data)
                           },
                           holdout = function(){
                             return(simplify2array(private$.holdout))
                           },
                           keptin = function(){
                             return(simplify2array(private$.keptin))
                           },
                           inputs = function(){
                             return(private$.inputs)
                           }
                         ))