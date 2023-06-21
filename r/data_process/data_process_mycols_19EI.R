# Extract from 19EI myco data
library(data.table)
# myco19 <- readxl::read_excel('data/raw/19EI/mycols 010421.xlsx')
# myco19add <- readxl::read_excel('data/raw/19EI/req2004chithu2.xlsx')
setDT(myco19); setDT(myco19add)

myco19sumAdd = myco19add[,`:=`(
  GeneXpert = fifelse(GeneXpert == "0" | GeneXpert ==  "NOT DONE", NA_character_, GeneXpert)
)][, .(
  USUBJID = `19EICODE`,
  DateSample = as.Date(NA),
  ZN = ZielNeelsen == "Positive",
  MycoResult = (MGITculture == "Positive"),
  XpertResult = stringr::str_detect(GeneXpert, "MTB DETECTED"),
  XpertLevel = tolower(fifelse(stringr::str_detect(GeneXpert, "MTB DETECTED"), gsub("MTB DETECTED ", "", GeneXpert), NA_character_))
)]

myco19sum = myco19[, `:=`(
  GeneXpert = fifelse(GeneXpert != "0" & GeneXpert !=  "NOT DONE", GeneXpert, NA_character_),
  Xpert_Ultra = fifelse(Xpert_Ultra != "0" & Xpert_Ultra !=  "NOT DONE", Xpert_Ultra, NA_character_))][!USUBJID %in% myco19sumAdd$USUBJID, .(
  USUBJID,
  DateSample = as.Date(DateSample),
  ZN = (ZielNeelsen == "Positive"),
  MycoResult = (MGITculture == "Positive"),
  XpertResult = fifelse(is.na(Xpert_Ultra), stringr::str_detect(GeneXpert, "MTB DETECTED"), stringr::str_detect(Xpert_Ultra, "MTB DETECTED")),
  XpertLevel = tolower(fifelse(stringr::str_detect(GeneXpert, "MTB DETECTED") %in% TRUE, gsub("MTB DETECTED ", "", GeneXpert), 
                      fifelse(stringr::str_detect(Xpert_Ultra, "MTB DETECTED") %in% TRUE, gsub("MTB DETECTED ", "", Xpert_Ultra), NA_character_)))
)]

myco19sum = rbind(myco19sum, myco19sumAdd)

saveRDS(myco19sum, "data/cleaned//mycol19.RDS")
