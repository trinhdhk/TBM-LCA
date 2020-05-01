# Definitions of necessary functions of use
# Author: trinhdhk
# Vers 0.1.2004

# Get the nearest date that have tests
get_nearest <- function(x, deltaDate){
  assertthat::assert_that(length(x) == length(deltaDate),msg = 'Length mismatched!')
  require(data.table)
  dt <- data.table(x, deltaDate)
  sort(dt[!is.na(x) & deltaDate==min(deltaDate)]$x,decreasing = TRUE)[1]
}