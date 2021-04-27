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
data_19EI <- data_19EI[!USUBJID %in% c("003-335", "003-102")] # death?

## Add some known test result from other sources and remove unknown one (for now)
data_19EI <- data_19EI[!USUBJID %in% c("003-102",
                                       '003-625', 
                                       '003-407', '003-512')]

#data_19EI[USUBJID=="003-144", `:=`(csf_smear = 1, csf_mgit = 1, csf_xpert = 0)]
#data_19EI[USUBJID=="003-425", `:=`(csf_smear = 0, csf_mgit = 1, csf_xpert = 0)]
#data_19EI[USUBJID=="003-450", `:=`(csf_smear = 1, csf_mgit = 1, csf_xpert = 1)]
#data_19EI[USUBJID=="003-456", `:=`(csf_smear = 0, csf_mgit = 0, csf_xpert = 0)]
#data_19EI[USUBJID=="003-469", `:=`(csf_smear = 1, csf_mgit = 1, csf_xpert = 1)]
#data_19EI[USUBJID=="003-501", `:=`(csf_smear = 1, csf_mgit = 1, csf_xpert = 1)]
#data_19EI[USUBJID=="003-516", `:=`(csf_smear = 1, csf_mgit = 0, csf_xpert = 1)]
#data_19EI[USUBJID=="003-519", `:=`(csf_smear = 1, csf_mgit = 1, csf_xpert = 0)]
#data_19EI[USUBJID=="003-520", `:=`(csf_smear = 0, csf_mgit = 0, csf_xpert = 0)]
#data_19EI[USUBJID=="003-537", `:=`(csf_smear = 0, csf_mgit = 0, csf_xpert = 0)]
#data_19EI[USUBJID=="003-542", `:=`(csf_smear = 1, csf_mgit = 1, csf_xpert = 0)]

Xd <- data_19EI %$% cbind(
  hiv_stat,                             #1
  clin_symptoms,                        #2
  clin_motor_palsy,                     #3
  clin_nerve_palsy,                     #4
  clin_contact_tb,                      #5
  xray_pul_tb,                          #6
  xray_miliary_tb                       #7
)

Xc <- data_19EI %$% cbind(
  age=(age-mean(age, na.rm=TRUE))/10,   #1    #8
  id=log2(clin_illness_day),            #2    #9
  glu=sqrt(BLDGLU),                     #3    #10 
  csfglu=sqrt(csf_glucose),             #4    #11
  csflym=log2(csf_lympho+1),            #5    #12   
  csfpro=log2(csf_protein),             #6    #13   
  csflac=log2(csf_lactate),             #7    #14   
  gcs=15-clin_gcs                       #8    #16             
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
  log(WHITE), 
  log(LYMP))

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

