# Load libraries ----
library(magrittr)
library(data.table)
library(dplyr)

# Prepare input ----
data_19EI <- readRDS('data/cleaned/data_19EI_3.RDS')
data_19EI[, `:=`(
  csf_smear = ifelse(is.na(csf_smear), 0, csf_smear),
  csf_mgit =  ifelse(is.na(csf_mgit),  0, csf_mgit),
  csf_xpert = ifelse(is.na(csf_xpert), 0, csf_xpert)
)]
data_19EI[is.na(GCSV) & GLASCOW == 15, GCSV := 5]
data_19EI[, clin_contact_tb := clin_contact_tb %in% TRUE]
data_19EI <- data_19EI[!USUBJID %in% c("003-306","003-335", "003-102")] # death?
## Add some known test result from other sources and remove unknown one (for now)
data_19EI <- data_19EI[!USUBJID %in% c('003-407')]
data_19EI[USUBJID == "003-038", c("csf_smear", "csf_mgit", "csf_xpert", "XpertLevel") := list(1,1,1, "very low")] 
data_19EI[USUBJID == "003-100", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,1,0)] 
data_19EI[USUBJID == "003-132", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,0,0)] 
data_19EI[USUBJID == "003-178", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,0,0)] 
data_19EI[USUBJID == "003-512", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,1,0)] 
data_19EI[USUBJID == "003-625", c("csf_smear", "csf_mgit", "csf_xpert", "XpertLevel") := list(1,1,1, "very low")] #mgit is not done but get 1
data_19EI[USUBJID == "003-551", c("csf_smear", "csf_mgit", "csf_xpert", "XpertLevel") := list(1,1,1, "very low")]
data_19EI[USUBJID == "003-557", c("csf_smear", "csf_mgit", "csf_xpert") := list(1,0,0)]
#create dummy variable of XpertLevel
data_19EI[, c("csf_bac_verylow", "csf_bac_low", "csf_bac_med") := 
    #list((XpertLevel %in% c("very low", "low", "medium")) %in% TRUE, 
     # (XpertLevel %in% c("low", "medium")) %in% TRUE,
     # (XpertLevel %in% c("medium")) %in% TRUE)]
    list((XpertLevel == "very low") %in% TRUE,
      (XpertLevel == "low") %in% TRUE,
      (XpertLevel == "medium") %in% TRUE)]
#remove patient with age=0 bc outlier
data_19EI = data_19EI[age!=0]
#impute gcsv using: https://pubmed.ncbi.nlm.nih.gov/32898843/
#if missing all, single imputation using mice
data_19EI[, GCSV:=fifelse(is.na(GCSV), 
                         fcase(GCSE+GCSM<=6, 1,
                               GCSE+GCSM==7, 2,
                               GCSE+GCSM >7, 4),
                               #GCSE+GCSM<=9, 4,  #NO
                               #GCSE+GCSM==10,5), #NO BECAUSE ISREDUCED==TRUE
                         GCSV)]
dummi_19EI <- mice::mice(data_19EI[,.(hiv_stat, GCSE, GCSM, GCSV, 
  age, WHITE, LYMP, NEUTRO, EOSI, PLATE)], m=1, maxit=0, seed = 219)
pred_mat <- dummi_19EI$predictorMatrix
pred_mat[2:4, 6:10] <- 0
pred_mat[6:10, 2:4] <- 0
pred_mat[1,] <- 0

mi_19EI <- mice::mice(data_19EI[,.(hiv_stat, GCSE, GCSM, GCSV, 
  age, WHITE, LYMP, NEUTRO, EOSI, PLATE)], m=1, maxit=1500, seed = 219, predictorMatrix = pred_mat)
imp_19EI <- mice::complete(mi_19EI)[c("GCSE", "GCSM", "GCSV", "WHITE", "LYMP", "NEUTRO", "EOSI", "PLATE")]
data_19EI <- cbind(data_19EI[,
  c("GCSE", "GCSM", "GCSV", "WHITE", "LYMP", "NEUTRO", "EOSI", "PLATE") := NULL], imp_19EI)
# impgcs_19EI <- lapply(impgcs_19EI, function(x) 
#   cbind(row = row.names(x), x)[row.names(x) %in% which(is.na(data_19EI$clin_gcs)),])
# impgcs_19EI <- cbind(impgcs_19EI[[1]], impgcs_19EI[[2]][,2], impgcs_19EI[[3]][,2])
# data_19EI[is.na(GCSE)&is.na(GCSM)&is.na(GCSV)&is.na(clin_gcs),
#           c("GCSE", "GCSM", "GCSV"):=impgcs_19EI[,2:4]]
# data_19EI[,clin_gcs:=fifelse(is.na(clin_gcs), GCSE+GCSM+GCSV, clin_gcs)]

#Red blood cell lumbar puncture might need correction
#There is interaction but let's loosen it b/c we are on log scale
data_19EI[,csf_rbc:=ifelse(is.na(REDCELL),0,REDCELL)]

Xd <- data_19EI %$% cbind(
  hiv_stat,                             #1
  clin_symptoms,                        #2
  clin_motor_palsy,                     #3
  clin_nerve_palsy,                     #4
  clin_contact_tb,                      #5
  xray_pul_tb,                          #6
  xray_miliary_tb,                      #7
  csf_ink | csf_crypto
)

Xc <- data_19EI %$% cbind(
  #age=(age-mean(age, na.rm=TRUE))/10,   #1    #8
  age=log2(age)-mean(log2(age)),                      #1    #8 
  id=log2(clin_illness_day),            #2    #9
  glu=log10(BLDGLU),                     #3    #10 
  csfglu=log10(1+csf_glucose),           #4    #11
  csflym=log10(csf_lympho+1),            #5    #12   
  csfpro=log10(csf_protein),             #6    #13   
  csflac=log10(csf_lactate),             #7    #14   
  gcs=15-clin_gcs,                      #8    #15
  csfeos=log10(csf_eos+1),               #9    #16
  csfred=log10(REDCELL+1)                #10   #17
)

Td <- data_19EI %$% cbind(
  ISWEIGHT,                             #1
  ISNSWEAT,                             #2
  ISCOUGH,                              #3
  HEMIPLEGIA,                           #4
  PARAPLEGIA,                           #5
  TETRAPLEGIA,                          #6
  ISDIABETE                             #7
)


Tc <- data_19EI %$% cbind(
  4-GCSE, 
  6-GCSM, 
  5-GCSV, 
  log2(WHITE), 
  log2(LYMP),
  csf_bac_verylow, csf_bac_low, csf_bac_med
)

is.not.na <- Negate(is.na)
obs_Xd <- apply(Xd, 2, is.not.na) %>% apply(2, as.numeric)
obs_Xc <- apply(Xc, 2, is.not.na) %>% apply(2, as.numeric)
obs_Td <- apply(Td, 2, is.not.na) %>% apply(2, as.numeric)
obs_Tc <- apply(Tc, 2, is.not.na) %>% apply(2, as.numeric)

my.replace_na <- function(data, replace=0){
  apply(data, 2, tidyr::replace_na, replace=replace)
}

Xd <- my.replace_na(Xd)
Xc <- my.replace_na(Xc)
Td <- my.replace_na(Td)
Tc <- my.replace_na(Tc)

# Prepare the input ----
save(data_19EI, Xc, Xd, Td, Tc, obs_Xc, obs_Xd, obs_Td, obs_Tc, file='data/cleaned/data_input.Rdata')

