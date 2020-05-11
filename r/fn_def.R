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
  if (length(save_file)) saveRDS(result, file=save_file)
  result
}
