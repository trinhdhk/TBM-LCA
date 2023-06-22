# An Latent Class Diagnostic Model for *Tuberculous meningitis* (#TBM-LCA)

## Introduction

This is the repository for the TBM-LCA project. The main objectives of the project are to (1) provide a statistically calibrated diagnostic model for TBM patients and (2) provide an estimation of current myco-bacterial tests and guide for future ones.

## Shiny app

An Shiny app is available at: ![https://trinhdong.shinyapps.io/TBMLCA](https://trinhdong.shinyapps.io/TBMLCA). Source code for the Shiny app is available at: ![github.com/trinhdhk/Shiny.TBMLCA](github.com/trinhdhk/Shiny.TBMLCA) (to be updated).

## Reproduction

The main program to fit the LCA model is `sample.R`. The program providing the simplication is `simplify.R`.

Prerequisites:

-   R ( \>= 4.1)

-   Packages used: 

    ```
    abind_1.4-5             askpass_1.1             assertthat_0.2.1       
    backports_1.4.1         base64enc_0.1.3         BayesFactor_0.9.12.4.3 
    bayesplot_1.9.0         bayestestR_0.12.1       BH_1.78.0.0            
    bit_4.0.4               bit64_4.0.5             bitops_1.0.7           
    brio_1.1.3              broom_0.8.0             bslib_0.3.1            
    callr_3.7.0             caret_6.0-92            caTools_1.18.2         
    cellranger_1.1.0        checkmate_2.1.0         class_7.3-20           
    classifierplots_1.4.0   cli_3.0.1               cmdstanr_0.5.2         
    coda_0.19-4             codetools_0.2-18        colorspace_2.0-3       
    compiler_4.2.0          contfrac_1.1.12         coro_1.0.2             
    correlation_0.8.1       cpp11_0.4.2             crayon_1.5.1           
    curl_4.3.2              data.table_1.14.2       datasets_4.2.0         
    datawizard_0.4.1        DBI_1.1.2               dcurves_0.3.0          
    desc_1.4.1              deSolve_1.32            diffobj_0.3.5          
    digest_0.6.29           distributional_0.3.0    dplyr_1.0.9            
    e1071_1.7.9             effectsize_0.6.0.1      ellipsis_0.3.2         
    elliptic_1.4.0          emmeans_1.7.4-1         estimability_1.3       
    evaluate_0.15           fansi_1.0.3             farver_2.1.0           
    fastmap_1.1.0           foreach_1.5.2           fs_1.5.2               
    future_1.25.0           future.apply_1.9.0      generics_0.1.2         
    ggfx_1.0.0              gghighlight_0.3.2       ggplot2_3.3.6          
    ggrepel_0.9.1           ggridges_0.5.3          ggsignif_0.6.3         
    ggstatsplot_0.9.3       ggtext_0.1.1            globals_0.15.0         
    glue_1.6.2              gower_1.0.0             gplots_3.1.3           
    graphics_4.2.0          grDevices_4.2.0         grid_4.2.0             
    gridExtra_2.3           gridtext_0.1.4          gtable_0.3.0           
    gtools_3.9.2.1          hardhat_0.2.0           highr_0.9              
    hms_1.1.1               htmltools_0.5.2         hypergeo_1.2.13        
    import_1.2.0            inline_0.3.19           insight_0.17.1         
    ipred_0.9-12            isoband_0.2.5           iterators_1.0.14       
    jpeg_0.1.9              jquerylib_0.1.4         jsonlite_1.8.0         
    KernSmooth_2.23.20      knitr_1.39              labeling_0.4.2         
    LaplacesDemon_16.1.6    lattice_0.20-45         lava_1.6.10            
    lifecycle_1.0.1         listenv_0.8.0           loo_2.5.1              
    lubridate_1.8.0         magick_2.7.3            magrittr_2.0.3         
    markdown_1.1            MASS_7.3-57             Matrix_1.3-4           
    MatrixModels_0.5.0      matrixStats_0.62.0      mc2d_0.1.21            
    methods_4.2.0           mgcv_1.8.40             mime_0.12              
    ModelMetrics_1.2.2.2    multcomp_1.4-19         munsell_0.5.0          
    mvtnorm_1.1-3           nlme_3.1-157            nnet_7.3-17            
    numDeriv_2016.8.1.1     openssl_2.0.2           packrat_0.7.0          
    paletteer_1.4.0         pander_0.6.5            parallel_4.2.0         
    parallelly_1.31.1       parameters_0.18.1       patchwork_1.1.1        
    pbapply_1.5.0           pbmcapply_1.5.1         performance_0.9.0      
    pillar_1.7.0            pkgbuild_1.3.1          pkgconfig_2.0.3        
    pkgload_1.2.4           plyr_1.8.7              png_0.1-7              
    posterior_1.2.1         praise_1.0.0            prettyunits_1.1.1      
    prismatic_1.1.0         pROC_1.18.0             processx_3.5.3         
    prodlim_2019.11.13      progress_1.2.2          progressr_0.10.0       
    proxy_0.4.26            ps_1.7.0                purrr_0.3.4            
    R6_2.5.1                ragg_1.2.2              rappdirs_0.3.3         
    RColorBrewer_1.1-3      Rcpp_1.0.8.3            RcppEigen_0.3.3.9.2    
    RcppParallel_5.1.5      RCurl_1.98.1.6          readxl_1.4.0           
    recipes_0.2.0           rematch_1.0.1           rematch2_2.1.2         
    reshape_0.8.9           reshape2_1.4.4          rlang_1.0.2            
    rmarkdown_2.14          rmda_1.6                ROCR_1.0-11            
    rpart_4.1.16            rprojroot_2.0.3         rsconnect_0.8.25       
    rstan_2.29.2.9000       rstudioapi_0.13         sandwich_3.0-1         
    sass_0.4.1              scales_1.2.0            splines_4.2.0          
    SQUAREM_2021.1          StanHeaders_2.29.2.9000 stats_4.2.0            
    stats4_4.2.0            statsExpressions_1.3.2  stringi_1.7.6          
    stringr_1.4.0           survival_3.2-12         sys_3.4                
    systemfonts_1.0.4       tensorA_0.36.2         
    testthat_3.1.4          textshaping_0.3.6       TH.data_1.1-1          
    tibble_3.1.7            tidyr_1.2.0             tidyselect_1.1.2       
    timeDate_3043.102       tinytex_0.39            tools_4.2.0            
    torch_0.8.0             tzdb_0.3.0              utf8_1.2.2             
    utils_4.2.0             V8_4.2.0                vctrs_0.4.1            
    viridisLite_0.4.0       vroom_1.5.7             waldo_0.4.0            
    withr_2.5.0             WRS2_1.1.3              xfun_0.31              
    xml2_1.3.3              xtable_1.8-4            yaml_2.3.5             
    zeallot_0.1.0           zoo_1.8-10             
    ```

Run the program:

-   Unix(-like) operating systems:

    ``` bash
    cd /path/to/TBM-LCA
    ./sampler.R --config-file sampler-config.yml
    ```

-   Windows:

    ``` cmd
    cd /path/to/TBM-LCA
    Rscript sampler.R --config-file sampler-config.yml
    ```
    
To seek help use flag `--help`.

&copy; Trinh Huu-Khanh Dong \@ Biostatistics Group - Oxford University Clinical Research Unit - Viet Nam 2022
