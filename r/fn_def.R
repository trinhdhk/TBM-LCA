# Definitions of necessary functions of use
# Author: trinhdhk
# Vers 0.1.2004

#Automatically import excel files
import.Excel <- function(dataPath){
  requireNamespace('readxl')
  sheets <- readxl::excel_sheets(dataPath)
  tables <- sapply(sheets, 
                   function(sheet) {
                     tryCatch(
                       readxl::read_excel(dataPath, sheet = sheet, guess_max = 3000),
                       warning = function(w) {
                         print(sheet)
                         readxl::read_excel(dataPath, sheet = sheet, guess_max = 3000)
                       })
                   }, 
                   simplify = FALSE,
                   USE.NAMES = TRUE)
  return(list(tables = tables, tableNames = sheets))
}

#auto loo
quickloo <- function(model, cores=getOption('mc.cores'), ...){
  log_lik_ <- loo::extract_log_lik(model, merge_chains = FALSE)
  r_eff_ <- loo::relative_eff(exp(log_lik_), cores=cores)
  loo_ <- loo::loo(log_lik_, r_eff = r_eff_, cores=cores, ...)
  loo_
}

# Get the nearest date that have tests
get_nearest <- function(x, deltaDate){
  assertthat::assert_that(length(x) == length(deltaDate),msg = 'Length mismatched!')
  require(data.table)
  dt <- data.table(x, deltaDate)
  sort(dt[!is.na(x) & deltaDate==min(deltaDate)]$x,decreasing = TRUE)[1]
}

# Load model
load_model <- function(id, data, chain=getOption('mc.cores'),..., save_file = file.path('outputs', paste0('results_',id,'.RDS'))){
  model <- stan_model(file.path('stan', paste0('model_', id, '.stan')))
  result <- sampling(model, data=data, chain=chain, ...)
  if (length(save_file)) try(saveRDS(result, file=save_file))
  result
}
