# Load libraries ----
library(magrittr)
library(data.table)
library(dplyr)

# Prepare input ----
data_19EI <- readRDS('data/cleaned/data_19EI_3.RDS')

data_19EI[is.na(GCSV) & GLASCOW == 15, GCSV := 5]
data_19EI[is.na(GCSE) & GLASCOW == 15, GCSE := 4]
data_19EI[is.na(GCSM) & GLASCOW == 15, GCSM := 6]
data_19EI[, clin_contact_tb := clin_contact_tb %in% TRUE]
# data_19EI <- data_19EI[!USUBJID %in% c("003-102")] # death?
# data_19EI <- data_19EI[!USUBJID %in% c("003-306","003-335", "003-102")] # death?
## Add some known test result from other sources and remove unknown one (for now)

# data_19EI <- data_19EI[!USUBJID %in% c('003-407')] # no information on TBM group side
# data_19EI[USUBJID == "003-038", c("csf_smear", "csf_mgit", "csf_xpert", "XpertLevel") := list(1,1,1, "very low")] 
# data_19EI[USUBJID == "003-100", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,1,0)] 
# data_19EI[USUBJID == "003-132", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,0,0)] 
# data_19EI[USUBJID == "003-178", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,0,0)] 
# data_19EI[USUBJID == "003-512", c("csf_smear", "csf_mgit", "csf_xpert") := list(0,1,0)] 
# data_19EI[USUBJID == "003-625", c("csf_smear", "csf_mgit", "csf_xpert", "XpertLevel") := list(1,1,1, "very low")] #mgit is not done but get 1
# data_19EI[USUBJID == "003-551", c("csf_smear", "csf_mgit", "csf_xpert", "XpertLevel") := list(1,1,1, "very low")]
# data_19EI[USUBJID == "003-557", c("csf_smear", "csf_mgit", "csf_xpert") := list(1,0,0)]
# #create dummy variable of XpertLevel
# data_19EI[, c("csf_bac_verylow", "csf_bac_low", "csf_bac_med") := 
#             list((XpertLevel == "very low") %in% TRUE,
#                  (XpertLevel == "low") %in% TRUE,
#                  (XpertLevel == "medium") %in% TRUE)]

# Load cleaned data export from TB Group
test_data = readRDS('data/cleaned/full_tb_test_extracted_good.RDS')
# test_data = readRDS('data/cleaned/full_tb_test_extracted_good.RDS')
setDT(test_data)
#remove age=0 bc outlier, replace it with the one in the test data
new_age = test_data[PatientCode==data_19EI[age==0, USUBJID],Age]
data_19EI[, age := ifelse(age==0, new_age, age)]
# data_19EI = merge(data_19EI, test_data[,.(USUBJID=PatientCode, Contaminated)], all.x=TRUE)
#Red blood cell lumbar puncture might need correction
#There is interaction but let's loosen it b/c we are on log scale
data_19EI[,csf_rbc:=ifelse(is.na(REDCELL),0,REDCELL)]

#join the test results & remove conataminated
data_19EI = dplyr::left_join(
  data_19EI |> dplyr::select(-csf_xpert, -csf_mgit, -csf_smear, XpertLevel_EI=XpertLevel),
  test_data |> dplyr::select(USUBJID=PatientCode,
                             csf_smear=ZielNeelsen,
                             csf_mgit=MGITculture,
                             csf_xpert=GeneXpert,
                             csf_mgit_contaminated = Contaminated,
                             Volume, diffday, wrong_name, GrowthUnit, TimeToPositive, XpertLevel),
  by='USUBJID') |> filter(USUBJID != '003-335') #not really in the study

#csf smear should not be missing if xpert or mgit is available. they were forgotten to input
data_19EI[, csf_smear := ifelse(is.na(csf_smear)&(!is.na(csf_mgit)|!is.na(csf_xpert)), FALSE, csf_smear)]
data_19EI[USUBJID == '003-357', age := 30]

saveRDS(data_19EI, 'export/data_dirty.RDS')

# filter out only those with uncontaminated CSF, Volume >= 3, diffday <= 7
data_19EI = data_19EI |>
  filter(is.na(csf_mgit_contaminated) | !csf_mgit_contaminated) |>
  filter(is.na(Volume) | Volume >= 3) |>
  filter(is.na(diffday) | diffday < 7) |>
  filter(!wrong_name)
  # filter(!USUBJID%in% c('003-102', '003-407'))  #death

data_19EI[, `:=`(
  Volume    = fifelse(is.na(Volume), 6, Volume),
  obs_smear = fifelse(is.na(csf_smear), 0, 1        ),
  obs_mgit  = fifelse(is.na(csf_mgit ), 0, 1        ),
  obs_xpert = fifelse(is.na(csf_xpert), 0, 1        ),
  csf_smear = fifelse(is.na(csf_smear), F, csf_smear),
  csf_mgit  = fifelse(is.na(csf_mgit ), F, csf_mgit ),
  csf_xpert = fifelse(is.na(csf_xpert), F, csf_xpert)
)]

# Correction for traumatic lumbar puncture
# If csf_rbc >= 10
debloody = function(csf_cell, bld_cell, csf_rbc, bld_rbc){
  ifelse(csf_rbc >= 10, sapply(csf_cell - bld_cell * csf_rbc / bld_rbc, max, 0), csf_cell)
}
data_19EI[, `:=`(
  csf_neutro = csf_wbc - csf_lympho - csf_eos)][, `:=`(
  csf_wbc_corrected = debloody(csf_wbc, WHITE*1000, csf_rbc, HEMO*1e6),
  csf_lympho_corrected = debloody(csf_lympho, LYMP/100*WHITE*1000, csf_rbc, HEMO*1e6),
  csf_neutro_corrected = debloody(csf_neutro, NEUTRO/100*WHITE*1000, csf_rbc, HEMO*1e6),
  csf_eos_corrected = debloody(csf_eos, EOSI/100*WHITE*1000, csf_rbc, HEMO*1e6),
  csf_protein_corrected = csf_protein - 0.011 * csf_rbc / 1000
)]
################ Not used
#impute gcsv using: https://pubmed.ncbi.nlm.nih.gov/32898843/
#if missing all, single imputation using mice
# data_19EI[, GCSV:=fifelse(is.na(GCSV), 
#                           fcase(GCSE+GCSM<=6, 1,
#                                 GCSE+GCSM==7, 2,
#                                 GCSE+GCSM >7, 4),
#                           #GCSE+GCSM<=9, 4,  #NO
#                           #GCSE+GCSM==10,5), #NO BECAUSE ISREDUCED==TRUE
#                           GCSV)]
# dummi_19EI <- mice::mice(data_19EI[,.(hiv_stat, GCSE, GCSM, GCSV, 
#                                       age, WHITE, LYMP, NEUTRO, EOSI, PLATE)], m=1, maxit=0, seed = 219)
# pred_mat <- dummi_19EI$predictorMatrix
# pred_mat[2:4, 6:10] <- 0
# pred_mat[6:10, 2:4] <- 0
# pred_mat[1,-c(6,7)] <- 0
# 
# mi_19EI <- mice::mice(data_19EI[,.(hiv_stat, GCSE, GCSM, GCSV, 
#                                    age, WHITE, LYMP, NEUTRO, EOSI, PLATE,csf_wbc)], m=1, maxit=500, seed = 219, predictorMatrix = pred_mat)
# imp_19EI <- mice::complete(mi_19EI)[c("GCSE", "GCSM", "GCSV", "WHITE", "LYMP", "NEUTRO", "EOSI", "PLATE")]
# data_19EI <- cbind(data_19EI[,
#                              c("GCSE", "GCSM", "GCSV", "WHITE", "LYMP", "NEUTRO", "EOSI", "PLATE") := NULL], imp_19EI) # "csf_wbc"
##################### End not used

Xd <- data_19EI %$% cbind(
  hiv_stat,                                                                 #1
  clin_symptoms,                                                            #2
  clin_motor_palsy,                                                         #3
  clin_nerve_palsy,                                                         #4
  clin_contact_tb,                                                          #5
  xray_pul_tb,                                                              #6
  xray_miliary_tb,                                                          #7
  csf_ink | csf_crypto ,                                                    #8
  GRAM = csf_gram,                                                          #9
  eos = csf_eos_corrected > 0                                               #10
)

scale.nonzero = function(x, ...) {
  x[x==0] = NA;
  scale(x, ...)
}
Xc <- data_19EI %$% tibble(
  id=scale(log2(clin_illness_day),scale=T),                          #1    #11
  glu=scale(log2(BLDGLU), scale=T),                                  #2    #12 
  csfglu=scale(log2(.1+csf_glucose), scale=T),                       #3    #13
  csflym=scale(log10(csf_lympho_corrected+1), scale=T),              #4    #14   
  csfpro=scale(log2(csf_protein_corrected), scale=T),                #5    #15   
  csflac=scale(log2(csf_lactate), scale=T),                          #6    #16   
  csfwbc=scale(log10(csf_wbc_corrected+1), scale=T),                 #7    #17
  gcs=structure((15-clin_gcs - 3)/3, dim=c(nrow(data_19EI),1), `scaled:center`=3, `scaled:scale`=3),  #8    #18
  csfeos=scale.nonzero(log10(csf_eos_corrected+1), scale=T),                 #9   #19
  csfred=scale(log10(csf_rbc+1), scale=T)                            #10   #20
)

scale_Xc <- summarise(Xc, across(everything(), attributes)) |> as.list() 

Td <- data_19EI %$% cbind(
  ISWEIGHT,                                                           #1
  ISNSWEAT,                                                           #2
  ISCOUGH,                                                            #3
  HEMIPLEGIA,                                                         #4
  PARAPLEGIA,                                                         #5
  TETRAPLEGIA,                                                        #6
  suspectedHIV = DISDIA %in% c('DIA1', 'DIA10', 'DIA12', 'DIA3') | ISTBTREAT == 'C49488' | !is.na(TBTREATDATE),     #7
  sex = sex                                                           #8
)
 
#Td[Td[,6]==1, 4:5] <- 1 #everyone has tetra must have semi and para
#Td[Td[,4]==0&Td[,5]==0, 6] <- 0 #''
Xd[,3] = apply(Td[,4:6],1, any)
Xd[,2] = apply(Td[,1:3],1, any)


Tc <- data_19EI %$% cbind(
  (4-GCSE)/3,                                                         #1
  (6-GCSM)/5,                                                         #2
  (5-GCSV)/4,                                                         #3
  scale(log2(WHITE), scale=T),                                        #4
  scale(log2(LYMP), scale=T)                                          #5
  # csf_bac_verylow, csf_bac_low, csf_bac_med 
)

D = data_19EI %$% cbind(
  scale(Volume/2, scale=T)
)

# Record the missingness and Set all missing to 0
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
D  <- my.replace_na(D, 6)

# Prepare the input ----
save(data_19EI, Xc, Xd, Td, Tc, obs_Xc, obs_Xd, obs_Td, obs_Tc, D, scale_Xc,
     file='data/cleaned/data_input.Rdata')
