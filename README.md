# An Latent Class Diagnostic Model for *Tuberculous meningitis* (#TBM-LCA)

## Introduction

This is the repository for the TBM-LCA project. The main objectives of the project are to (1) provide a statistically calibrated diagnostic model for TBM patients and (2) provide an estimation of current myco-bacterial tests and guide for future ones.

## Reproduction

The main program to fit the LCA model is `sample.R`. The program providing the simplication is `simplify.R`.

Prerequisites:

-   R ( \>= 4.1)

-   Packages: `argparser, data.table, dplyr, Stan (>=2.27), abind, ggplot2, bayesplot, rmda, ggfx, R6, pROC, classifierplots, future, furrr`

Run the program:

-   Unix and Unix-like operating systems:

    ``` bash
    cd TBM-LCA
    ./sampler.R --config-file sample-config.yml
    ```

-   Windows:

    ``` cmd
    cd TBM-LCA
    Rscript sampler.R --config-file sample-config.yml
    ```

@copy; Trinh Huu-Khanh Dong \@ Biostatistics Group - Oxford University Clinical Research Unit - Viet Nam 2022
