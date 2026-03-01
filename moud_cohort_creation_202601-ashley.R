# data management only
setwd("/Volumes/HPLWork1/p50_oud/Ashley/Ahra2") 

library(Hmisc)
library(tangram)
library(readr)
library(haven)
library(dplyr)
library(tidyr)
library(purrr)
# library(ggplot2)
# library(ggh4x)
# library(plotly)
# library(data.table)
library(reshape2)
library(janitor)
library(kableExtra)
# library(UpSetR)
# library(survival)
# library(rms)
# library(mice)
# library(PSW)
library(readxl)

load("/Volumes/HPLWork1/p50_oud/Ashley/Ahra2/G5_data_202509.Rdata")
ndc <- read.csv("/Volumes/HPLWork1/p50_oud/Ashley/Ahra2/NDCs_2020.csv", colClasses = c(ndcnum = "character"))

# exclude from cohort if they have these meds: H0020, G2067, and G2078. 

af_extract <- as.data.frame(af_extract)
extract <- af_extract[which(af_extract$ndc %in% c(ndc$ndcnum, 'Q9991',
                                                  'Q9992',
                                                  'J0571',
                                                  'J0572',
                                                  'J0573',
                                                  'J0574',
                                                  'J0575',
                                                  'G2068',
                                                  'G2069',
                                                  'G2070',
                                                  'G2072',
                                                  'G2079',
                                                  'J0575',
                                                  'G2068',
                                                  'G2069',
                                                  'G2070',
                                                  'G2079',
                                                  'G2072',
                                                  'J2315',
                                                  'G2073',
                                                  'G2057',
                                                  'G2067',
                                                  'G2078',
                                                  'S0109')), ]
methadonereceipt <- af_extract[af_extract$ndc %in% c("H0020", "G2067", "G2078"), ]
extract <- extract[!extract$StudyID %in% methadonereceipt$StudyID, ]

# AL: 
cat("Records with G2067 or G2078 in extract:",
    sum(extract$ndc %in% c("G2067", "G2078")), "\n")
# Returns 0 matches in current data -- TennCare did not cover OTP methadone until June 2020

extract <- merge(extract, ndc, by.x="ndc", by.y= "ndcnum", all.x = TRUE)
extract <- extract[extract$BTrack==90, ]
# extract$gennme <- ifelse(is.na(extract$gennme) & extract$ndc =="H0020", "Methadone", extract$gennme)
# extract$daysup <- ifelse(is.na(extract$daysup) & extract$ndc =="H0020", 1, extract$daysup)

# subset people who are in the exposure group
# basic <- af_master_new[af_master_new$StudyID %in% extract$StudyID, ]
# basic <- af_master[af_master$Expose1==1, ] # 998421 to 31199

# AL CHANGES
# Starting cohort: OUD evidence via diagnosis (Expose1==1) OR MOUD medication in pharmacy file
# Union ensures MOUD users not captured by Expose1 DX flag are included
basic <- af_master[af_master$BTrack==90 &
                     (af_master$Expose1==1 | af_master$StudyID %in% extract$StudyID), ]
n_step1 <- nrow(basic) #this union added 1,983 pregnancies beyond the original 15,951 = 17,934 total (Box 1 in flow diagram)

basic$filledrx <- ifelse(basic$StudyID %in% extract$StudyID, 1, 0)
table(basic$filledrx) #9,818 have filled MOUD prescription captured in extract; 8,116 do not have an MOUD prescription in extract-but they could still have OUD evidence via Expose1==1 (diagnosis or medication captured through a different route)
  
# exclude people if they have an overdose or death from -90 to 41 days
# not a problem in the new data; everyone has outcome of 42+ days
# af_outcome$dead_before_follow <- ifelse(!is.na(af_outcome$death_day) & af_outcome$death_day <=41, 1, 
#                              ifelse(!is.na(af_outcome$overdose_day) & af_outcome$overdose_day <= 41, 1, 0))
# dead_before_follow <- af_outcome[af_outcome$dead_before_follow==1, ]
# af_outcome <- af_outcome[af_outcome$dead_before_follow==0 , ]

# want 90 day only since 180 is for sensitivity

af_master$unique_id <- sub(".*-", "", af_master$StudyID)
basic$unique_id <- sub(".*-", "", basic$StudyID)

# merge covariates
af_cov$unique_id <- sub(".*-", "", af_cov$StudyID)
# only include one of the pregnancies if twins/triplets 
basic <- basic[order(basic$unique_id, basic$StudyID), ]
# separate those with multiple pregnancies
# create indicators: if next row covariates are the same as previous row then they are twins/triplets
basic <- basic %>% mutate(rownum= row_number()) %>% data.frame()
# possibly twins? 
# these are the records where the study id is duplicated (either twins or multiple children over years)
twins <- basic %>% group_by(unique_id) %>% 
  mutate(n = n()) %>% # count number of same study id
  filter(n> 1) %>%
  mutate(mom_dob_prevrow = lag(MomDOB_day), 
         gweeks20_prevrow = lag(GWeeks20_dt_day))
twins$multbirths <- ifelse(twins$mom_dob_prevrow == twins$MomDOB_day & twins$gweeks20_prevrow == twins$GWeeks20_dt_day, 1, 0)
twins$multbirths <- ifelse(is.na(twins$multbirths), 0, twins$multbirths)  
not_twins <- twins[twins$multbirths==0 ,]
twins2 <- twins[twins$multbirths==1,] # these are twins/triplets

# exclude twins/triplets
basic <- basic[! (basic$rownum %in% twins2$rownum),  ] # n=15328
# AL ADDED:
n_step2      <- nrow(basic)
n_excl_twins <- n_step1 - n_step2
cat("Unique pregnancies after excluding twins/triplets:", formatC(n_step2, format = "d", big.mark = ","),
    "| Excluded:", n_excl_twins, "\n")
# Result from above code: Unique pregnancies after excluding twins/triplets: 17,233 | Excluded: 701
# Note: the twins/triplets exclusion removes duplicate rows within a single pregnancy 
# --> where the same mother, same DOB, same gestational week date appears more than once
# because she delivered multiples--we want to keep one record per pregnancy***

# AL ADDED:
# Exclude individuals with any overdose during the exposure/baseline window (-90 to +41 days)                      
# Applied early as an eligibility criterion, using af_outcome directly (af_out not yet created)                    
# Note: basic still has uppercase StudyID at this point                                                            
od_exposure_ids <- unique(                                                                                         
  af_outcome$StudyID[                                                                                              
    af_outcome$BTrack == 90 &
      af_outcome$StudyID %in% basic$StudyID &
      (af_outcome$Overdose_Opioid_exposure == 1 | af_outcome$Overdose_NonOpioid_exposure == 1)
  ]
)
basic              <- basic[!(basic$StudyID %in% od_exposure_ids), ]
n_step3            <- nrow(basic)
n_excl_od_exposure <- n_step2 - n_step3
cat("After excluding overdose in exposure window:", formatC(n_step3, format = "d", big.mark = ","),
    "| Excluded:", n_excl_od_exposure, "\n")
##Results from above: After excluding overdose in exposure window: 17,221 | Excluded: 12 


# some covariates are in Btrack 180 only and not 90
af_cov <- as.data.frame(af_cov)
# af_cov2 <- af_cov[af_cov$BTrack==90, ]
af_cov2 <- af_cov[af_cov$StudyID %in% basic$StudyID, ]

af_cov_u <- af_cov2[ , !(colnames(af_cov2) %in% c("BTrack", "condition", "condition_from_day", "CodeID"))]
af_cov_u <- af_cov_u[!duplicated(af_cov_u), ]

basic <- merge(basic, af_cov_u, by= "StudyID", all.x = TRUE)

# df <- basic #creating a copy to work off of just in case
names(basic) <- tolower(names(basic))

# duplicated ids 
# dups <- basic[which(duplicated(basic$unique_id)), "unique_id"]
# dup <- basic[which(duplicated(basic$unique_id)), ]
# dupdf <- basic[which(basic$unique_id %in% dups), ]
# include only the first pregnancy for now 

# Format covariates
label(basic$momage) <- "Age at delivery"
basic$momrace <- ifelse(basic$momrace=="", NA, basic$momrace)
basic$momrace_desc <- ifelse(basic$momrace_desc=="", NA, basic$momrace_desc)
basic$momethnicity_desc <- ifelse(basic$momethnicity_desc=="", NA, basic$momethnicity_desc)
basic$momrace.new <- ifelse(basic$momrace_desc %in% c("White", "Bridged White"), "White",
                            ifelse(basic$momrace_desc %in% c("Black or African-American", "Bridged Black"), "Black of African American",
                                   ifelse(basic$momrace_desc %in% c("American Indian or Alaskan Native", "Bridged American Indian or Alaskan Native"), "American Indian or Alaska Native",
                                          ifelse(basic$momrace_desc %in% c("Asian Indian", "Bridged Asian or Pacific Isander", "Chinese", "Filipino", "Korean", "Native Hawaiian", "Other Asian", "Vietnamese"), "Asian or Pacific Islander",
                                                 ifelse(basic$momrace_desc %in% c("Other Race", "Refused or Unknown"), "Other or Unknown", NA)))))

label(basic$momrace.new) <- "Race"

basic$momhisp.new <- ifelse(basic$momethnicity_desc %in% c("\tNot Hispanic", "\tOther", "\tUnknown"), "Not Hispanic", 
                            ifelse(is.na(basic$momethnicity_desc), NA, "Hispanic"))

basic$momraceeth.new <- ifelse(basic$momhisp.new=="Hispanic", "Hispanic",
                               ifelse((basic$momhisp.new!="Hispanic"| is.na(basic$momhisp.new)) & basic$momrace.new=="Asian or Pacific Islander", "Non-Hispanic Asian",
                                      ifelse((basic$momhisp.new!="Hispanic"| is.na(basic$momhisp.new)) & basic$momrace.new=="Black of African American", "Non-Hispanic Black",
                                             ifelse((basic$momhisp.new!="Hispanic"| is.na(basic$momhisp.new)) & basic$momrace.new=="White", "Non-Hispanic White", NA))))
basic$momraceeth.new <- ifelse(basic$momhisp.new=="Not Hispanic" & is.na(basic$momrace.new) & is.na(basic$momraceeth.new), "Non-Hispanic Other", 
                               ifelse(is.na(basic$momhisp.new) & basic$momrace.new=="White" & is.na(basic$momraceeth.new), "Non-Hispanic White",
                                      ifelse(is.na(basic$momhisp.new) & basic$momrace.new %in% c("American Indian or Alaska Native", "Other or Unknown"), "Non-Hispanic Other", 
                                             ifelse(basic$momhisp.new=="Not Hispanic" & basic$momrace.new %in% c("American Indian or Alaska Native", "Other or Unknown") & is.na(basic$momraceeth.new), "Non-Hispanic Other", basic$momraceeth.new))))
label(basic$momraceeth.new) <- "Race/Ethnicity combined"

label(basic$rpl_themes) <- "CDC SVI overall percentile ranking of all 4 themes"

#  information related to each delivery/pregnancy
basic$routedeliv.new <- ifelse(basic$deliverycesarean=="Y", "Cesarean", ifelse(basic$deliveryvaginal=="Y", "Vaginal", ifelse(basic$deliveryforcep=="Y", "Vaginal", ifelse(basic$deliveryvacuum=="Y", "Vaginal", NA))))
basic$routedeliv.new <- as.factor(basic$routedeliv.new)
label(basic$routedeliv.new) <- "Route of delivery"

# AL CHECK
key_covs <- c("routedeliv.new", "deliverycesarean", "deliveryvaginal",                                             
              "momrace_desc", "momethnicity_desc", "momage",                                                       
              "rpl_themes", "previouspregnanciestotaln")                                                           

# NAs per key covariate
na_by_cov <- sapply(key_covs, function(v) {
  if (v %in% names(basic)) sum(is.na(basic[[v]])) else NA_integer_
})
print(na_by_cov)

# Which IDs are missing at least one key covariate?
key_cols_present <- intersect(key_covs, names(basic))
missing_any_key <- basic$studyid[
  apply(basic[, key_cols_present], 1, function(x) any(is.na(x)))
]
cat("StudyIDs missing at least one key covariate:", length(missing_any_key), "\n")


# AL ADDITION: Check if any missingness in the delivery category
# Check af_cov directly - filter to current basic StudyIDs, BTrack==90                                       
af_cov_check <- af_cov[af_cov$StudyID %in% basic$studyid & af_cov$BTrack == 90, ]                                  

# Deduplicate since af_cov has multiple rows per StudyID (one per condition)                                       
af_cov_check_u <- af_cov_check[!duplicated(af_cov_check$StudyID), ]                                                

# How many have all four delivery fields not "Y"?                                                                  
sum(!af_cov_check_u$DeliveryCesarean %in% "Y" &
      !af_cov_check_u$DeliveryVaginal  %in% "Y" &
      !af_cov_check_u$DeliveryForcep   %in% "Y" &
      !af_cov_check_u$DeliveryVacuum   %in% "Y")

# How many StudyIDs in basic have no record in af_cov at all?
sum(!basic$studyid %in% af_cov$StudyID)
# Results above: 0 & 0, meaning everyone is in the af_cov and everyone has at least one delivery route marked "Y"


# AL COMMENTED OUT; NO NEED TO EXCLUDE HERE, DO AT ANALYSIS STAGE SINCE OTHER COVARIATES ARE ALSO MISSING, NOT JUST 63 that are missing for all 
# These were excluded previously if ALL covariates were missing for individuals, but in reality, some people have some but others 
# that are missing/etc. and those are handled at the analysis stage; makes sense to handle all at the analysis stage. 
# exclude if route of delivery everything is marked as "NO"
# excluded one person
#basic <- basic[!is.na(basic$routedeliv.new), ] #15328 to 15272
# AL ADDITIONS: 
#n_step4      <- nrow(basic)                                                                                        
#n_excl_route <- n_step3 - n_step4                                                                                  
#cat("After excluding missing route of delivery:", formatC(n_step4, format = "d", big.mark = ","),                  
#    "| Excluded:", n_excl_route, "\n")                                                                           


label(basic$previouspregnanciestotaln) <- "Number of previous pregnancies"
label(basic$cat4_seg1) <- "Enrollment up to 365 days post-partum"

basic$enroll365.new.f <- factor(basic$cat4_seg1, levels = c(0,1), labels=c("No", "Yes"))
label(basic$enroll365.new.f) <- "Enrollment up to 365 days post-partum"

label(basic$momhisp.new) <- "Hispanic"

# Death
# want to capture deaths from day 42 and onwards but before loss of enrollment; some deaths captured after they disenroll (for now, but will go back and capture them later)
basic$death.new <- ifelse(!is.na(basic$momdeath_day) & basic$momdeath_day <= 323 & basic$momdeath_day <= basic$fdate_seg1_day, 1, 0)
basic$death.new.f <- factor(basic$death.new, levels=c(0, 1), labels=c("No", "Yes"))
label(basic$death.new.f) <- "Death during follow-up only (not after disenroll)"

basic$death.new.2 <- ifelse(!is.na(basic$momdeath_day) & basic$momdeath_day <= 323 , 1, 0)
basic$death.new.2.f <- factor(basic$death.new.2, levels=c(0, 1), labels=c("No", "Yes"))
label(basic$death.new.f) <- "Death within 1 year (even after disenroll)"

# to check deaths
# dead323 <- basic[which(basic$cat2_seg1==1 & basic$momdeath_day<= 323), ]

# comorbidities
af_cov2 <- af_cov2[af_cov2$condition!="", ]
af_cov2 <- af_cov2[af_cov2$BTrack==90, ]

# transform conditions from long to wide
af_cov3 <- af_cov2[ , c("StudyID", "condition_from_day", "condition")]

# comorb <- af_cov3 %>% group_by(StudyID) %>%
#   pivot_wider(id_cols = "StudyID", names_from = "condition", values_from = "condition") %>%
#   mutate_all(~replace(., . == "NULL", 0)) %>%
#   mutate_all(~replace(., . != "0", 1)) %>% 
#   mutate_if(is.list, as.numeric) %>%
#   ungroup() %>% as.data.frame()

comorb <- af_cov3 %>% group_by(StudyID) %>%
  pivot_wider(id_cols = "StudyID", names_from = "condition", values_from = "condition") %>%
  mutate(across(everything(), ~replace(., . == "NULL", 0))) %>%
  mutate(across(everything(), ~replace(., . != "0", 1))) %>% 
  mutate_if(is.list, as.numeric) %>%
  ungroup() %>% as.data.frame()

basic <- merge(basic, comorb, by.x = "studyid", by.y = "StudyID", all.x = TRUE)
## AL: This merge code makes sense; what didn't make sense was the subsequent code after this in the original
# script: "basic <- basic[basic$studyid %in% comorb$StudyID, ]", which I commented out. 
# See check below. 

# remove those with missing comorbidities
# basic <- basic[basic$studyid %in% comorb$StudyID, ] # original code line 202; AL commented this line out; 

# AL check: 
# Those removed by the original line 285--though Ahra meant to just remove the 63 without any of the birth
# certificate variables. 
# AL: handle missingness at the analysis stage; We want all people for unadjusted analysis

# removed_by_285 <- basic$studyid[!basic$studyid %in% comorb$StudyID]
# cat("Total that original line 285 would remove:", length(removed_by_285), "\n")

# Truly missing from af_cov vs. in af_cov but no conditions
# absent_from_afcov <- removed_by_285[!toupper(removed_by_285) %in% toupper(af_cov$StudyID)]
# in_afcov_no_cond  <- removed_by_285[toupper(removed_by_285) %in% toupper(af_cov$StudyID)]

# cat("  Absent from af_cov entirely (truly missing):", length(absent_from_afcov), "\n")
# cat("  In af_cov but zero comorbidities (healthy, wrongly removed):", length(in_afcov_no_cond), "\n")

# meant to exclude 63 missing birth certificate data but excluded an additional 527 who had no comorbidities
# (depression/anxiety/other substance use/etc.)

comorbidities <- names(comorb)[-1] #AL changed to -1, removing manual index update, will capture all
comorblabs <- gsub("_", " ", comorbidities)
comorblabs <- capitalize(comorblabs)
comorblabs <- gsub("nsaids", "NSAIDS", comorblabs)
comorblabs <- gsub("substance", "substance use", comorblabs)
comorblabs <- gsub("Bipolard", "Bipolar disorder", comorblabs)
comorblabs <- gsub("psych", "psychotic disorders", comorblabs)
comorblabs <- gsub("benzo", "benzodiazepine use", comorblabs)
comorblabs <- gsub("usedis", "use disorder", comorblabs)
comorblabs <- gsub("anticonv", "anticonvulsant use", comorblabs)
comorblabs <- gsub("Personality dis", "Personality disorder", comorblabs)
comorblabs <- gsub("Alc", "Alcohol", comorblabs)
comorblabs <- gsub("Mat", "Maternal", comorblabs)
comorblabs <- gsub("relax", "relaxant use", comorblabs)
comorblabs <- gsub("Med barbit", "Med barbiturates", comorblabs)

# barbiturates are new in data as of January 2026 after removing Exposure1 criteria

# create factor variables
for(i in 1:length(comorbidities)) {
  basic[[paste0(comorbidities[i], ".f")]] <- basic[[comorbidities[i]]]
  basic[[paste0(comorbidities[i], ".f")]] <- ifelse(is.na(basic[[paste0(comorbidities[i], ".f")]]), 0, basic[[paste0(comorbidities[i], ".f")]])
  basic[[paste0(comorbidities[i], ".f")]] <- factor(basic[[paste0(comorbidities[i], ".f")]], levels = c(0:1), labels=c('No', 'Yes'))
  label(basic[[paste0(comorbidities[i], ".f")]]) <- comorblabs[i]
}

label(basic$gweeks) <- "Number of Gestation weeks"

# Combine variables for analyses
basic$pain.combine <- ifelse(basic$back_pain.f=="Yes" | basic$dental_pain.f=="Yes" | basic$fibromyalgia.f=="Yes" | basic$med_nsaids.f=="Yes" | basic$med_muscle_relax.f=="Yes", 1, 0)
basic$pain.combine.f <- factor(basic$pain.combine, levels = 0:1, labels=c("No", "Yes"))
label(basic$pain.combine.f) <- "Combined pain"

basic$personality.combine <- ifelse(basic$personality_dis.f=="Yes"| basic$schizophrenia_psych.f=="Yes", 1, 0)
basic$personality.combine.f <- factor(basic$personality.combine, levels = 0:1, labels=c("No", "Yes"))
label(basic$personality.combine.f) <- "Schizophrenia or personality disorder"

basic$alcuse.othersub <- ifelse(basic$alc_usedis.f=="Yes" | basic$other_substance.f=="Yes", 1, 0)
basic$alcuse.othersub.f <- factor(basic$alcuse.othersub, levels = 0:1, labels=c("No", "Yes"))
label(basic$alcuse.othersub.f) <- "Alcohol or other substance use"

basic$benzo.anticonv <- ifelse(basic$med_benzo.f=="Yes" | basic$med_anticonv.f == "Yes", 1, 0)
basic$benzo.anticonv.f <- factor(basic$benzo.anticonv, levels=0:1, labels=c("No", "Yes"))
label(basic$benzo.anticonv.f) <- "Benzodiazepine or anticonvulsant use"

basic$anxiety <- ifelse(is.na(basic$anxiety), 0, basic$anxiety)
basic$depression <- ifelse(is.na(basic$depression), 0, basic$depression)
basic$mat_morbidity <- ifelse(is.na(basic$mat_morbidity), 0, basic$mat_morbidity)

# overdose
af_out <- af_outcome[af_outcome$BTrack==90 & af_outcome$StudyID %in% basic$studyid, ]
names(af_out) <- tolower(names(af_out))
# change to count outcomes that occur after disenrollment too
af_out$fdate_seg1_day <- basic[match(af_out$studyid, basic$studyid), "fdate_seg1_day"]
# check to see which overdose or death occurs after disenrollment
# x <- af_out[which(af_out$overdose_day - 42 > af_out$fdate_seg1_day), ]
# subset events that occured before disenrollment (before)
af_out$momdeath_day <- basic[match(af_out$studyid, basic$studyid), "momdeath_day"]
# af_out <- af_out[(!is.na(af_out$overdose_day) & (af_out$overdose_day <= af_out$fdate_seg1_day)) | (!is.na(af_out$momdeath_day) & (af_out$momdeath_day<= af_out$fdate_seg1_day)), ]
af_out$momdeath_day <- ifelse(af_out$momdeath_day > 323, NA, af_out$momdeath_day)

# outcomes file contains multiple overdose records (previous ODs can show up on next visits)
# summarize information on total number of overdoses
af_out2 <- af_out %>% group_by(studyid) %>%
  filter(overdose_opioid_outcome ==1  | overdose_nonopioid_outcome ==1) %>%
  mutate(n= n()) %>% data.frame()

af_out2$unique_id <- sub(".*-", "", af_out2$studyid)

# AL Note for below code: exposure==1 conditions below are largely redundant since
# individuals with any exposure-window overdose were already excluded from "basic" earlier.

# exclude if they had OD in both exposure and outcome
af_out2 <- af_out2[!(af_out2$overdose_nonopioid_exposure==1 & af_out2$overdose_nonopioid_outcome==1), ]
af_out2 <- af_out2[!(af_out2$overdose_opioid_exposure==1 & af_out2$overdose_opioid_outcome==1), ]
af_out2 <- af_out2[!(af_out2$overdose_nonopioid_exposure==1 & af_out2$overdose_opioid_outcome==1), ]
af_out2 <- af_out2[!(af_out2$overdose_opioid_exposure==1 & af_out2$overdose_nonopioid_outcome==1), ]

af_out2 <- af_out2[af_out2$overdose_nonopiod_day >= 0 & af_out2$overdose_nonopiod_day <= 323 | 
                     af_out2$overdose_opiod_day >= 0 & af_out2$overdose_opiod_day <= 323,  ]


# If someone who had multiple births had overdoses for both births, then both information will be included (but this did not happen in the data)


# af_out2 <- af_out %>% group_by(studyid) %>%
#   filter(!is.na(overdose_day)) %>%
# #   mutate(n= n()) %>% data.frame()
# 
# basic$num.od <- af_out2[match(basic$studyid, af_out2$studyid), "n"]
# label(basic$num.od) <- "Number of any overdoses recorded during follow-up"

# Filter out opioid related overdose only
# af_out_op <- af_out %>% group_by(studyid) %>%
#   filter(overdose_opioid_outcome==1) %>%
#   mutate(n= n()) %>% data.frame()
# 
# basic$num.opioid.od <- af_out_op[match(basic$studyid, af_out_op$studyid), "n"]
# label(basic$num.opioid.od) <- "Number of opioid overdoses recorded during follow-up"  

# ever had od during follow-up
# basic$outcome.od <- ifelse(!is.na(basic$num.od), 1, 0)
# label(basic$outcome.od) <- "Had any overdose during follow-up"

# ever had opioid od during follow-up
# basic$outcome.opioid.od <- ifelse(!is.na(basic$num.opioid.od), 1, 0)
# label(basic$outcome.opioid.od) <- "Had any opioid overdose during follow-up"
#   

# More deaths from outcome table; merge into large table
af_out_death <- af_out %>% group_by(studyid) %>% 
  filter(!is.na(momdeath_day) & momdeath_day <= 323) %>% data.frame()

af_out_death <- af_out_death[ , c("studyid", grep("dcause", names(af_out_death), value=TRUE), "momdeath_day")]
af_out_death <- unique(af_out_death)
# basic <- merge(basic, af_out_death, by="studyid", all.x = TRUE)
basic$death_day <- basic$momdeath_day
basic$death_day <- ifelse(basic$death_day > 323, NA, basic$death_day)


# get dosage information
# redbook <- read_sas("/Volumes/HPDrdw/MarketScan/RawData/CCAE/Random_1Pct/SAS Data 2022/redbook_yr22.sas7bdat", NULL)
# names(redbook) <- tolower(names(redbook))
# redbook$gennme2 <- tolower(redbook$gennme)
# redbook$methadone <- ifelse(grepl("methadone", redbook$gennme2), 1, 0)
# redbook_methadone <- redbook[redbook$methadone==1, ]
# write.csv(redbook_methadone, file="redbook_methadone.csv", row.names = FALSE)
redbook_methadone <- read.csv("/Volumes/HPLWork1/p50_oud/Ashley/Ahra2/redbook_methadone.csv")
redbook_methadone$ndcnum  <-sprintf("%0*s", 11, as.character(redbook_methadone$ndcnum)) 

# redbook2 <- redbook[redbook$ndcnum %in% extract$ndc, ]
# write.csv(redbook2, file="redbook2.csv", row.names = FALSE)
redbook2 <- read.csv("/Volumes/HPLWork1/p50_oud/Ashley/Ahra2/redbook2.csv")
redbook2$ndcnum <-sprintf("%0*s", 11, as.character(redbook2$ndcnum)) 
extract <- merge(extract, redbook2, by.x="ndc", by.y="ndcnum", all.x = TRUE)
names(extract)[names(extract)=="gennme.x"] <- "gennme"
names(extract)[names(extract)=="prodnme.x"] <- "prodnme"


# make sure the individual is excluded if they have used methadone for OUD treatment                                                                         
# AL: restricted to exposure window (BTrack==90, Rx_B90_DOB==1 or Rx_DOB_P41==1)                                                                             
# AL: HCPCS codes only — redbook NDCs removed as these capture methadone for pain, not OUD
methadone <- af_extract[af_extract$BTrack==90 &
                          (af_extract$Rx_B90_DOB==1 | af_extract$Rx_DOB_P41==1) &
                          af_extract$ndc %in% c("H0020", "G2067", "G2078", "S0109", "J1230"), ]
# extract <- extract[!(extract$ndc=="H0020" | extract$ndc %in% redbook_methadone$ndcnum), ]
basic <- basic[!(basic$studyid %in% methadone$StudyID), ]

# AL ADDED:
n_step4          <- nrow(basic)
n_excl_methadone <- n_step3 - n_step4
cat("After excluding methadone users:", formatC(n_step4, format = "d", big.mark = ","),
    "| Excluded:", n_excl_methadone, "\n")
## Excluded = 0
## This was an error in the original code; people were dropped who shouldn't have been. 



# Unique pregnancies with BTrack==90 methadone
methadone_90_ids <- unique(methadone$StudyID[methadone$BTrack==90])
cat("Pregnancies with BTrack==90 methadone:", length(methadone_90_ids), "\n")

# Unique pregnancies with BTrack==180 methadone ONLY (no BTrack==90 record)
methadone_180_only_ids <- unique(methadone$StudyID[methadone$BTrack==180 &
                                                     !methadone$StudyID %in% methadone_90_ids])
cat("Pregnancies with BTrack==180 methadone only:", length(methadone_180_only_ids), "\n")




# -132 to -43 days = 90 days
# -42 to -1 days = 42 days

# get days of OUD medications
# 90 days before delivery
# extract90 <- extract[extract$Rx_B90_DOB==1 & extract$Expose1==1, ] # no longer subset on expose1 (per Andrew, January 2026)
extract90 <- extract[extract$Rx_B90_DOB==1, ]
names(extract90) <- tolower(names(extract90))
extract90 <- extract90[extract90$studyid %in% basic$studyid, ]
extract90 <- extract90[order(extract90$studyid, extract90$from), ]

# extract90 <- extract90 %>% group_by(studyid) %>% 
#                   mutate(earliest.rx = min(from),
#                         last.rx = max(from),
#                             next.rx.exp = as.Date(from) + daysup,
#                                 next.rx.actual = lead(from),
#                                   prev.rx = lag(from),
#                                     diff.days.rx = from - prev.rx) 


# if expected next prescription date (date when med will run out is the same in two rows, then delete the second row)
# extract90 <- extract90 %>% group_by(studyid) %>% 
#                   mutate(keep = case_when(prev.rx + lag(daysup) > from + diff.days.rx ~ 0))
# 
# extract90 <- extract90 %>% group_by(studyid) %>% 
#                   mutate(keep = case_when(next.rx.exp == lag(next.rx.exp) ~ 0,
#                                           from == earliest.rx ~ 1, 
#                                           from == last.rx ~ 1))

# extract90$keep <- ifelse(is.na(extract90$keep) & )
# IDs to check 
# extract90[which(extract90$studyid=="002009032034-04544465-03494573"), c("studyid", "daysup", "from", "total.rx.days")]
# "002008060134-05537053-05774374"
# "002008032250-05512799-02853395"

# summarize the number of days supplied cap at 90+42 days = 132 days (Method 1)
extract90 <- extract90 %>% group_by(studyid) %>% 
  mutate(total.rx.days = cumsum(daysup)) %>% data.frame()

#extract90[which(extract90$studyid=="002009032034-04544465-03494573"), c("studyid", "daysup", "from", "total.rx.days")]
extract90$total.rx.days <- ifelse(extract90$total.rx.days > 90, 90, extract90$total.rx.days)

extract90s <- extract90 %>% group_by(studyid) %>% slice_max(total.rx.days, with_ties = FALSE) %>% data.frame()


# 41 days after delivery
# extract41 <- extract[extract$Rx_DOB_P41==1 & extract$Expose1==1, ] (no longer needed per Andrew, January 2026)
extract41 <- extract[extract$Rx_DOB_P41==1, ]
names(extract41) <- tolower(names(extract41))
extract41 <- extract41[extract41$studyid %in% basic$studyid, ]
extract41 <- extract41[order(extract41$studyid, extract41$from), ]

extract41 <- extract41 %>% group_by(studyid) %>% 
  mutate(total.rx.days = cumsum(daysup)) %>% data.frame()
extract41$total.rx.days <- ifelse(extract41$total.rx.days > 42, 42, extract41$total.rx.days)

extract41s <- extract41 %>% group_by(studyid) %>% slice_max(total.rx.days, with_ties = FALSE) %>% data.frame()


basic$total.rx.days.b90 <- extract90s[match(basic$studyid, extract90s$studyid), "total.rx.days"]
basic$total.rx.days.p41 <- extract41s[match(basic$studyid, extract41s$studyid), "total.rx.days"]

basic$total.rx.days.b90 <- ifelse(is.na(basic$total.rx.days.b90), 0, basic$total.rx.days.b90)
basic$total.rx.days.b90 <- ifelse(basic$total.rx.days.b90 > 90, 90, basic$total.rx.days.b90)
basic$total.rx.days.p41 <- ifelse(is.na(basic$total.rx.days.p41), 0, basic$total.rx.days.p41)



# Method 2
# Add up total days supplied, if > 90 then add that number to days after delivery instead of capping at 90
extract90_2 <- extract90 %>% group_by(studyid) %>% 
  mutate(total.rx.days = cumsum(daysup),
         max.rx.days = max(total.rx.days)) %>% data.frame()

extract90_2$total.rx.minus90 <- ifelse(extract90_2$max.rx.days > 90, extract90_2$max.rx.days - 90, 0)


extract41_2 <- extract41 %>% group_by(studyid) %>% 
  mutate(total.rx.days = cumsum(daysup),
         max.rx.days = max(total.rx.days)) %>% data.frame()

basic$total.rx.days.b90_2 <- extract90_2[match(basic$studyid, extract90_2$studyid), "max.rx.days"]
basic$total.rx.days.b90_2 <- ifelse(basic$total.rx.days.b90_2 > 90, 90, basic$total.rx.days.b90_2)
basic$total.rx.days.b90_2 <- ifelse(is.na(basic$total.rx.days.b90_2), 0, basic$total.rx.days.b90_2)
basic$total.rx.days.b90_minus90 <- extract90_2[match(basic$studyid, extract90_2$studyid), "total.rx.minus90"]
basic$total.rx.days.b90_minus90 <- ifelse(is.na(basic$total.rx.days.b90_minus90), 0, basic$total.rx.days.b90_minus90)

basic$total.rx.days.p41_2max <- extract41_2[match(basic$studyid, extract41_2$studyid), "max.rx.days"]
basic$total.rx.days.p41_2max <- ifelse(basic$total.rx.days.p41_2max > 42, 42, basic$total.rx.days.p41_2max)
basic$total.rx.days.p41_2max <- ifelse(is.na(basic$total.rx.days.p41_2max), 0, basic$total.rx.days.p41_2max)

basic$total.rx.days.p41_2 <- rowSums(basic[, c("total.rx.days.p41_2max", "total.rx.days.b90_minus90")])
basic$total.rx.days.p41_2 <- ifelse(basic$total.rx.days.p41_2 > 42, 42, basic$total.rx.days.p41_2)

# Method 3
# counting for overlapping prescription dates
# extract90_3_old <- extract90 %>% group_by(studyid) %>%
#     arrange(studyid, from_day) %>%
#     mutate(indx = c(0, cumsum(as.numeric(lead(from_day)) >
#                      cummax(as.numeric(from_day + daysup)))[-n()])) %>%
#     group_by(studyid, indx) %>%
#     summarise(start = min(from_day), end = max(from_day + daysup))

extract90_3 <- extract90 
extract90_3$strength <- as.numeric(gsub( " .*$", "", extract90_3$strngth))
extract41_3 <- extract41
extract41_3$strength <- as.numeric(gsub( " .*$", "", extract41_3$strngth))

newdf <- bind_rows(extract90_3, extract41_3)
# exclude outlier where qty is 56000
newdf <- newdf[!newdf$qty==56000, ]
# exclude record where days supply is 0
newdf <- newdf[!newdf$daysup==0, ]

# powder buprenorphine ndc: 51927101200; exclude in days supply but not dosage calculation 
# powder bup was taken by someone with Expose1==0, so they are not included 
# exclude sublocade in dosage calculation, but they can be included in the days supply calculation
# same with vivitrol; calculate days supply but not the dosage

# who got injectable (solution) versus not
# who got buprenorphine versus not

newdf_dose <- newdf[!(newdf$prodnme %in% "sublocade 300mg"), ] #only 10 observations differ and 1 study id ("002020059663-06267662-03931614")
newdf_dose <- newdf_dose[newdf_dose$prodnme != "vivitrol", ] #differ by 78 study ids
newdf_dose <- newdf_dose[newdf_dose$prodnme != "sublocade 300mg", ] # this did not reduce anymore people


calc_coverage_and_dose <- function(df, pre_cap = 90, post_cap = 42) {
  
  # Expand ranges to unique day integers
  get_unique_days <- function(start_vec, end_vec) {
    sort(unique(unlist(purrr::map2(start_vec, end_vec, seq))))
  }
  
  # Safely coerce possible NULL/NA to integer(0) for set operations
  safe_int <- function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) integer(0) else as.integer(x)
  }
  
  df <- df %>%
    mutate(
      fill_start = from_day,
      fill_end   = from_day + daysup - 1
    )
  
  # ---------------- PRE-PERIOD ----------------
  pre_df <- df %>% filter(from_day >= -132, from_day <= -43) %>%
    mutate(
      fill_start = pmax(fill_start, -132),
      fill_end   = fill_start + daysup - 1
    )
  
  pre_summary <- pre_df %>%
    group_by(studyid) %>%
    summarise(
      all_covered_days = list(get_unique_days(fill_start, fill_end)),
      .groups = "drop"
    ) %>%
    mutate(
      pre_period_days   = map(all_covered_days, ~ .x[.x >= -132 & .x <= -43]),
      capped_days       = map(pre_period_days, ~ head(.x, pre_cap)),
      pre_days_covered  = map_int(capped_days, length),
      
      # days that didn't fit in capped pre; then restrict to post window for carryover
      carryover_day_list = map2(all_covered_days, capped_days, ~ setdiff(.x, .y)),
      carryover_days_vec = map(carryover_day_list, ~ .x[.x >= -42 & .x <= -1]),
      carryover_days     = map_int(carryover_days_vec, length)
    )
  
  capped_pre_df <- pre_df %>%
    inner_join(pre_summary, by = "studyid") %>%
    rowwise() %>%
    mutate(
      actual_days        = list(seq(fill_start, fill_end)),
      overlap_days       = length(intersect(actual_days, capped_days)),
      proportion         = overlap_days / daysup,
      adjusted_quantity  = qty * proportion,
      strength_quantity  = qty * strength * proportion
    ) %>%
    ungroup()
  
  pre_dose_summary <- capped_pre_df %>%
    group_by(studyid) %>%
    summarise(
      total_adjusted_quantity_pre = sum(adjusted_quantity, na.rm = TRUE),
      total_strength_quantity_pre = sum(strength_quantity, na.rm = TRUE),
      total_daysup_pre            = sum(daysup * proportion, na.rm = TRUE),
      avg_dose_pre = if_else(total_daysup_pre > 0,
                             total_strength_quantity_pre / total_daysup_pre, 0),
      .groups = "drop"
    )
  
  pre_summary <- pre_summary %>%
    left_join(pre_dose_summary, by = "studyid")
  
  # ---------------- POST-PERIOD ----------------
  post_df <- df %>% filter(from_day >= -42, from_day <= -1) %>%
    mutate(
      fill_start = pmax(fill_start, -42),
      fill_end   = fill_start + daysup - 1
    )
  
  post_summary <- post_df %>%
    group_by(studyid) %>%
    summarise(
      all_covered_days = list(get_unique_days(fill_start, fill_end)),
      .groups = "drop"
    ) %>%
    mutate(
      post_period_days  = map(all_covered_days, ~ .x[.x >= -42 & .x <= -1]),
      post_capped_days  = map(post_period_days, ~ head(.x, post_cap)),
      post_days_covered = map_int(post_capped_days, length),
      effective_days    = post_days_covered
    )
  
  capped_post_df <- post_df %>%
    inner_join(post_summary, by = "studyid") %>%
    rowwise() %>%
    mutate(
      actual_days        = list(seq(fill_start, fill_end)),
      overlap_days       = length(intersect(actual_days, post_capped_days)),
      proportion         = overlap_days / daysup,
      adjusted_quantity  = qty * proportion,
      strength_quantity  = qty * strength * proportion
    ) %>%
    ungroup()
  
  post_dose_summary <- capped_post_df %>%
    group_by(studyid) %>%
    summarise(
      total_adjusted_quantity_post = sum(adjusted_quantity, na.rm = TRUE),
      total_strength_quantity_post = sum(strength_quantity, na.rm = TRUE),
      total_daysup_post            = sum(daysup * proportion, na.rm = TRUE),
      avg_dose_post = if_else(total_daysup_post > 0,
                              total_strength_quantity_post / total_daysup_post, 0),
      .groups = "drop"
    )
  
  post_summary <- post_summary %>%
    left_join(post_dose_summary, by = "studyid")
  
  # ---------------- COMBINE PRE + POST ----------------
  full_summary <- pre_summary %>%
    select(studyid, capped_days, carryover_days_vec, pre_days_covered, carryover_days,
           total_adjusted_quantity_pre, total_strength_quantity_pre,
           total_daysup_pre, avg_dose_pre) %>%
    full_join(
      post_summary %>%
        select(studyid, post_capped_days, post_days_covered, effective_days,
               total_adjusted_quantity_post, total_strength_quantity_post,
               total_daysup_post, avg_dose_post),
      by = "studyid"
    ) %>%
    mutate(across(
      c(pre_days_covered, carryover_days, total_adjusted_quantity_pre,
        total_strength_quantity_pre, total_daysup_pre, avg_dose_pre,
        post_days_covered, effective_days, total_adjusted_quantity_post,
        total_strength_quantity_post, total_daysup_post, avg_dose_post),
      ~ replace_na(.x, 0)
    )) %>%
    mutate(
      total_days_covered           = pre_days_covered + post_days_covered,
      total_adjusted_quantity_all  = total_adjusted_quantity_pre + total_adjusted_quantity_post,
      total_strength_quantity_all  = total_strength_quantity_pre + total_strength_quantity_post,
      total_daysup_all             = total_daysup_pre + total_daysup_post,
      avg_dose_entire_period       = if_else(total_daysup_all > 0,
                                             total_strength_quantity_all / total_daysup_all, 0),
      
      # NEW: unique post-period days including carryover spillover
      total_unique_post_days_incl_carryover = map2_int(
        carryover_days_vec, post_capped_days,
        ~ length(unique(c(safe_int(.x), safe_int(.y))))
      ),
      
      # NEW: true total unique days supplied across periods (pre capped + carryover + post capped)
      total_unique_days_supplied = pmap_int(
        list(capped_days, carryover_days_vec, post_capped_days),
        function(a, b, c) length(unique(c(safe_int(a), safe_int(b), safe_int(c))))
      )
    ) %>%
    # drop helper list-columns from final output (keep if you want to audit)
    select(-capped_days, -carryover_days_vec, -post_capped_days)
  
  return(full_summary)
}

# manuscript example: id 002009013002-04514929-02081695 (old); now use: (002007035729-02903621-04746761)
# total_daysup_pre adds the overlapping/duplicate rx dates
# use pre_days_covered and post_days_covered (adjusted)
# need to take into account for carryover_days that do not overlap in the post-period and add that to post_days_covered

extract_method3 <- calc_coverage_and_dose(newdf)
extract_method3 <- as.data.frame(extract_method3)
extract <- extract[order(extract$StudyID, extract$from_day),]

# to describe dose use newdf_dose
extract_method3_dose <- calc_coverage_and_dose(newdf_dose)
extract_method3_dose <- as.data.frame(extract_method3_dose)


# check

# extract[which(extract$StudyID=="002012002626-01447011-05711550"), c("StudyID", "from_day", "daysup", "qty", "strngth")]
# extract[which(extract$StudyID=="002012081578-01529551-05299669"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]
# extract[which(extract$StudyID=="002007047985-02911169-05614828"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]
# extract[which(extract$StudyID=="002007035729-02903621-04746761"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]
# extract[which(extract$StudyID=="002007050401-02350074-04750282"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]
# extract[which(extract$StudyID=="002009013002-04514929-02081695"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]

# merge total days and dose
extract_method3 <- as.data.frame(extract_method3)
# pre-days covered
basic$total.rx.days.b90_3 <- extract_method3[match(basic$studyid, extract_method3$studyid), "pre_days_covered"]
# post-days covered
basic$total.rx.days.p41_3 <- extract_method3[match(basic$studyid, extract_method3$studyid), "total_unique_post_days_incl_carryover"]
# total days covered
# basic$total.rx.days.total_3 <- extract_method3[match(basic$studyid, extract_method3$studyid), "total_unique_days_supplied"]

#total days shifted
calc_coverage_and_dose_shifted <- function(df, pre_cap = 90, post_cap = 42) {
  
  library(dplyr)
  library(purrr)
  library(tidyr)
  
  # ---------------- HELPERS ----------------
  get_unique_days <- function(start_vec, end_vec) {
    sort(unique(unlist(map2(start_vec, end_vec, seq))))
  }
  
  safe_int <- function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) integer(0) else as.integer(x)
  }
  
  # NEW: shifted-days helper (sequential, no overlap)
  calc_shifted_days <- function(from_day, daysup,
                                window_start = -132,
                                window_end   = -1) {
    
    ord <- order(from_day)
    from_day <- from_day[ord]
    daysup   <- daysup[ord]
    
    current_day <- -Inf
    total_days  <- 0
    
    for (i in seq_along(from_day)) {
      start_day <- max(from_day[i], current_day)
      end_day   <- start_day + daysup[i] - 1
      
      usable_end <- min(end_day, window_end)
      usable_days <- max(0, usable_end - start_day + 1)
      
      total_days <- total_days + usable_days
      current_day <- end_day + 1
    }
    
    total_days
  }
  
  # ---------------- BASE DF ----------------
  df <- df %>%
    mutate(
      fill_start = from_day,
      fill_end   = from_day + daysup - 1
    )
  
  # ---------------- PRE-PERIOD ----------------
  pre_df <- df %>%
    filter(from_day >= -132, from_day <= -43) %>%
    mutate(
      fill_start = pmax(fill_start, -132),
      fill_end   = fill_start + daysup - 1
    )
  
  pre_summary <- pre_df %>%
    group_by(studyid) %>%
    summarise(
      all_covered_days = list(get_unique_days(fill_start, fill_end)),
      .groups = "drop"
    ) %>%
    mutate(
      pre_period_days   = map(all_covered_days, ~ .x[.x >= -132 & .x <= -43]),
      capped_days       = map(pre_period_days, ~ head(.x, pre_cap)),
      pre_days_covered  = map_int(capped_days, length),
      
      carryover_day_list = map2(all_covered_days, capped_days, ~ setdiff(.x, .y)),
      carryover_days_vec = map(carryover_day_list, ~ .x[.x >= -42 & .x <= -1]),
      carryover_days     = map_int(carryover_days_vec, length)
    )
  
  capped_pre_df <- pre_df %>%
    inner_join(pre_summary, by = "studyid") %>%
    rowwise() %>%
    mutate(
      actual_days       = list(seq(fill_start, fill_end)),
      overlap_days      = length(intersect(actual_days, capped_days)),
      proportion        = overlap_days / daysup,
      adjusted_quantity = qty * proportion,
      strength_quantity = qty * strength * proportion
    ) %>%
    ungroup()
  
  pre_dose_summary <- capped_pre_df %>%
    group_by(studyid) %>%
    summarise(
      total_adjusted_quantity_pre = sum(adjusted_quantity, na.rm = TRUE),
      total_strength_quantity_pre = sum(strength_quantity, na.rm = TRUE),
      total_daysup_pre            = sum(daysup * proportion, na.rm = TRUE),
      avg_dose_pre = if_else(total_daysup_pre > 0,
                             total_strength_quantity_pre / total_daysup_pre, 0),
      .groups = "drop"
    )
  
  pre_summary <- pre_summary %>%
    left_join(pre_dose_summary, by = "studyid")
  
  # ---------------- POST-PERIOD ----------------
  post_df <- df %>%
    filter(from_day >= -42, from_day <= -1) %>%
    mutate(
      fill_start = pmax(fill_start, -42),
      fill_end   = fill_start + daysup - 1
    )
  
  post_summary <- post_df %>%
    group_by(studyid) %>%
    summarise(
      all_covered_days = list(get_unique_days(fill_start, fill_end)),
      .groups = "drop"
    ) %>%
    mutate(
      post_period_days  = map(all_covered_days, ~ .x[.x >= -42 & .x <= -1]),
      post_capped_days  = map(post_period_days, ~ head(.x, post_cap)),
      post_days_covered = map_int(post_capped_days, length),
      effective_days    = post_days_covered
    )
  
  capped_post_df <- post_df %>%
    inner_join(post_summary, by = "studyid") %>%
    rowwise() %>%
    mutate(
      actual_days       = list(seq(fill_start, fill_end)),
      overlap_days      = length(intersect(actual_days, post_capped_days)),
      proportion        = overlap_days / daysup,
      adjusted_quantity = qty * proportion,
      strength_quantity = qty * strength * proportion
    ) %>%
    ungroup()
  
  post_dose_summary <- capped_post_df %>%
    group_by(studyid) %>%
    summarise(
      total_adjusted_quantity_post = sum(adjusted_quantity, na.rm = TRUE),
      total_strength_quantity_post = sum(strength_quantity, na.rm = TRUE),
      total_daysup_post            = sum(daysup * proportion, na.rm = TRUE),
      avg_dose_post = if_else(total_daysup_post > 0,
                              total_strength_quantity_post / total_daysup_post, 0),
      .groups = "drop"
    )
  
  post_summary <- post_summary %>%
    left_join(post_dose_summary, by = "studyid")
  
  # ---------------- SHIFTED DAYS (NEW) ----------------
  shifted_days <- df %>%
    group_by(studyid) %>%
    summarise(
      total_shifted_days_supplied =
        calc_shifted_days(from_day, daysup),
      .groups = "drop"
    )
  
  # ---------------- COMBINE ----------------
  full_summary <- pre_summary %>%
    select(studyid, capped_days, carryover_days_vec,
           pre_days_covered, carryover_days,
           total_adjusted_quantity_pre,
           total_strength_quantity_pre,
           total_daysup_pre, avg_dose_pre) %>%
    full_join(
      post_summary %>%
        select(studyid, post_capped_days,
               post_days_covered, effective_days,
               total_adjusted_quantity_post,
               total_strength_quantity_post,
               total_daysup_post, avg_dose_post),
      by = "studyid"
    ) %>%
    left_join(shifted_days, by = "studyid") %>%
    mutate(across(
      c(pre_days_covered, carryover_days,
        total_adjusted_quantity_pre,
        total_strength_quantity_pre,
        total_daysup_pre, avg_dose_pre,
        post_days_covered, effective_days,
        total_adjusted_quantity_post,
        total_strength_quantity_post,
        total_daysup_post, avg_dose_post,
        total_shifted_days_supplied),
      ~ replace_na(.x, 0)
    )) %>%
    mutate(
      total_days_covered          = pre_days_covered + post_days_covered,
      total_adjusted_quantity_all =
        total_adjusted_quantity_pre + total_adjusted_quantity_post,
      total_strength_quantity_all =
        total_strength_quantity_pre + total_strength_quantity_post,
      total_daysup_all =
        total_daysup_pre + total_daysup_post,
      avg_dose_entire_period =
        if_else(total_daysup_all > 0,
                total_strength_quantity_all / total_daysup_all, 0),
      
      total_unique_days_supplied = pmap_int(
        list(capped_days, carryover_days_vec, post_capped_days),
        function(a, b, c)
          length(unique(c(safe_int(a), safe_int(b), safe_int(c))))
      )
    ) %>%
    select(-capped_days, -carryover_days_vec, -post_capped_days)
  
  return(full_summary)
}

extract_method3_shifted <- calc_coverage_and_dose_shifted(newdf)


extract_method3_shifted <- as.data.frame(extract_method3_shifted)
basic$total.rx.days.shifted <- extract_method3_shifted[match(basic$studyid, extract_method3_shifted$studyid), "total_shifted_days_supplied"]
label(basic$total.rx.days.shifted) <- "MOUD days exposure shifted approach"

# labels
label(basic$total.rx.days.b90) <- "Total days supplied MOUD before delivery"
label(basic$total.rx.days.p41) <- "Total days supplied MOUD after delivery"

label(basic$total.rx.days.b90_2) <- "Total days supplied MOUD before delivery (method 2)"
label(basic$total.rx.days.p41_2) <- "Total days supplied MOUD after delivery (method 2)"

label(basic$total.rx.days.b90_3) <- "Total days supplied MOUD before delivery (method 3)"
label(basic$total.rx.days.p41_3) <- "Total days supplied MOUD after delivery (method 3)"


basic$avg.dose.pre <- extract_method3_dose[match(basic$studyid, extract_method3_dose$studyid), "avg_dose_pre"]
basic$avg.dose.post <- extract_method3_dose[match(basic$studyid, extract_method3_dose$studyid), "avg_dose_post"]
basic$avg.dose.total <- extract_method3_dose[match(basic$studyid, extract_method3_dose$studyid), "avg_dose_entire_period"]

label(basic$avg.dose.pre) <- "Avg dose pre-delivery (mg)"
label(basic$avg.dose.post) <- "Avg dose post-delivery (mg)"
label(basic$avg.dose.total) <- "Avg dose exposure (mg)"

# # summarize very last medication received per person
# extract_last <- extract %>% group_by(StudyID) %>% slice_max(from_day) %>% data.frame
# names(extract_last) <- tolower(names(extract_last))
# ## summarize the type of medication and number of prescriptions 
# extract_last <- extract_last[ , c("studyid", "from_day", "gennme")]
# extract_last <- unique(extract_last)
# extract_last <- reshape2::dcast(extract_last, studyid ~ gennme)
# extract_last$last_rx_type <- apply(extract_last[, 2:7], 1, function(x) toString(na.omit(x)))
# names(extract_last) <- make.names(names(extract_last))

## summarize the type of medication and number of prescriptions 
extract_u <- extract
names(extract_u) <- tolower(names(extract_u))
extract_u <- extract_u[ , c("studyid", "from_day", "gennme")]
extract_u <- unique(extract_u)
extract_u <- reshape2::dcast(extract_u, studyid ~ gennme)
names(extract_u) <- make.names(names(extract_u))
extract_u <- extract_u %>%  mutate_if(is.numeric, ~1 * (. > 0))
extract_u$unique.rx <- rowSums(extract_u[ , 2:6], na.rm = TRUE)
# filter out people in our dataset only 
extract_u <- extract_u[which(extract_u$studyid %in% basic$studyid),]
# excluded outlier with qty of 56000 so exclude this person 
extract_u <- extract_u[which(extract_u$studyid != "002020030193-06235149-04443003"), ]
# summarize how many took which medications
extract_u$bup <- ifelse(extract_u$Buprenorphine==1 | extract_u$Buprenorphine.Hydrochloride==1 | extract_u$Buprenorphine.Naloxone==1, 1, 0)
extract_u$nal <- ifelse(extract_u$Naltrexone==1 | extract_u$Naltrexone.Hydrochloride==1, 1, 0)
extract_u$rxtype <- ifelse(extract_u$bup==1 & extract_u$nal==0, "Buprenorphine only", 
                           ifelse(extract_u$nal==1 & extract_u$bup==0, "Naltrexone only",
                                  ifelse(extract_u$bup==1 & extract_u$nal==1, "Both", NA)))


# basic$last_rx_type <- extract_last[match(basic$studyid, extract_last$studyid), "last_rx_type"]
basic$unique.rx <- extract_u[match(basic$studyid, extract_u$studyid), "unique.rx"]
basic$rxtype <- extract_u[match(basic$studyid, extract_u$studyid), "rxtype"]
label(basic$rxtype) <- "Type of MOUD"

# 
# basic$last_rx_type <- as.factor(basic$last_rx_type)
# label(basic$last_rx_type) <- "Last medication type"

label(basic$unique.rx) <- "Unique number of prescriptions"


# dosage and exposure days information without naltrexone and methadone
newdf2 <- newdf
newdf2 <- newdf2[!(newdf2$gennme %in% c("Naltrexone", "Naltrexone Hydrochloride")), ]

newdf_dose2 <- newdf_dose
newdf_dose2 <- newdf_dose2[!(newdf_dose2$gennme %in% c("Naltrexone", "Naltrexone Hydrochloride")), ]

extract_method3_bup <- calc_coverage_and_dose(newdf2)
extract_method3_bup <- as.data.frame(extract_method3_bup)

extract_method3_bup_dose <- calc_coverage_and_dose(newdf_dose2)
extract_method3_bup_dose <- calc_coverage_and_dose(newdf_dose2)

# check

# extract[which(extract$StudyID=="002012002626-01447011-05711550"), c("StudyID", "from_day", "daysup", "qty", "strngth")]
# extract[which(extract$StudyID=="002012081578-01529551-05299669"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]
# extract[which(extract$StudyID=="002007047985-02911169-05614828"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]
# extract[which(extract$StudyID=="002018049083-01868012-03765766"), c("StudyID", "from_day", "daysup",  "qty", "strngth")]

# merge total days and dose
extract_method3_bup <- as.data.frame(extract_method3_bup)
# pre-days covered
basic$total.rx.days.b90_3bup <- extract_method3_bup[match(basic$studyid, extract_method3_bup$studyid), "pre_days_covered"]
# post-days covered
basic$total.rx.days.p41_3bup <- extract_method3_bup[match(basic$studyid, extract_method3_bup$studyid), "total_unique_post_days_incl_carryover"]
# total days covered
basic$total.rx.days.total_3bup <- extract_method3_bup[match(basic$studyid, extract_method3_bup$studyid), "total_unique_days_supplied"]

label(basic$total.rx.days.b90_3bup) <- "Total days buprenorphine before delivery (method 3)"
label(basic$total.rx.days.p41_3bup) <- "Total days buprenorphine after delivery (method 3)"

extract_method3_bup_dose <- as.data.frame(extract_method3_bup_dose)
basic$avg.dose.pre.bup <- extract_method3_bup_dose[match(basic$studyid, extract_method3_bup_dose$studyid), "avg_dose_pre"]
basic$avg.dose.post.bup <- extract_method3_bup_dose[match(basic$studyid, extract_method3_bup_dose$studyid), "avg_dose_post"]
basic$avg.dose.total.bup <- extract_method3_bup_dose[match(basic$studyid, extract_method3_bup_dose$studyid), "avg_dose_entire_period"]

label(basic$avg.dose.pre.bup) <- "Avg dose BUP pre-delivery (mg)"
label(basic$avg.dose.post.bup) <- "Avg dose BUP post-delivery (mg)"
label(basic$avg.dose.total.bup) <- "Avg dose BUP exposure (mg)"

# follow-up time if they had an event (counting events not captured in claims using all hospital payer data)
# if death or overdose happens after 365 days, then use overdose day or last observed day
af_out365 <- af_out[(af_out$momdeath_day >= 0 & af_out$momdeath_day <= 323) | (af_out$overdose_opiod_day >=0 & af_out$overdose_opiod_day <= 323 ) |  (af_out$fdate_seg1_day >=0 & af_out$fdate_seg1_day <= 323) | (af_out$overdose_nonopiod_day >=0 & af_out$overdose_nonopiod_day <= 323), ]

# exclude if they had OD in both exposure and outcome
af_out365 <- af_out365[!(af_out365$overdose_nonopioid_exposure==1 & af_out365$overdose_nonopioid_outcome==1), ]
af_out365 <- af_out365[!(af_out365$overdose_opioid_exposure==1 & af_out365$overdose_opioid_outcome==1), ]
af_out365 <- af_out365[!(af_out365$overdose_nonopioid_exposure==1 & af_out365$overdose_opioid_outcome==1), ]
af_out365 <- af_out365[!(af_out365$overdose_opioid_exposure==1 & af_out365$overdose_nonopioid_outcome==1), ]


af_out365$overdose_day <- pmin(af_out365$overdose_nonopiod_day, af_out365$overdose_opiod_day, na.rm = TRUE)
# af_out365 <- af_out365 %>% group_by(studyid) %>% slice_min(overdose_day, with_ties = FALSE) %>% data.frame()

af_out365$days.to.event <- pmin(af_out365$momdeath_day, af_out365$overdose_day, na.rm = TRUE)
af_out365$outcome <- ifelse(af_out365$days.to.event== af_out365$momdeath_day , "Death", 
                            ifelse(af_out365$days.to.event == af_out365$overdose_day, "Overdose", NA))
af_out365$outcome <- ifelse(is.na(af_out365$outcome) & is.na(af_out365$momdeath_day) & !is.na(af_out365$overdose_day), "Overdose", af_out365$outcome)

af_out365 <- as.data.frame(af_out365)
af_out <- as.data.frame(af_out)
af_out365$unique_id <- sub(".*-", "", af_out365$studyid)
af_out$unique_id <- sub(".*-", "", af_out$studyid)

# match on study ID to count overdoses that occurred in subsequent pregnancies
af_out$days.to.event <- af_out365[match(af_out$studyid, af_out365$studyid), "days.to.event"]
basic$days.to.event <- af_out365[match(basic$studyid, af_out365$studyid), "days.to.event"]
basic$outcome <- af_out365[match(basic$studyid, af_out365$studyid), "outcome"]

# af_out365[which(!is.na(af_out365$momdeath_day) & !is.na(af_out365$overdose_day)), c("studyid", "overdose_day", "fdate_seg1_day", "momdeath_day", "days.to.event", "outcome")]

# basic$days.to.event <- ifelse(is.na(basic$days.to.event), basic$fdate_seg1_day, basic$days.to.event)
# change so that if loss of enrollment occurs, we assume that no event happened if they are not found in death certificate or HDDS (change follow-up to 323 days or 1 year postpartum)
basic$days.to.event <- ifelse(is.na(basic$days.to.event), 323, basic$days.to.event)
basic$days.to.event <- ifelse(basic$days.to.event > 323, 323, basic$days.to.event)

# to check
# basic[which(basic$fdate_seg1_day < basic$days.to.event), c("studyid", "death.or.od",  "days.to.event", "fdate_seg1_day")]
# these people had fdate_seg1_day before death or overdose date; ???

# All dates are relative to Day 42, so everything is shifted to day-42
basic$death.day.from.deliv = basic$death_day + 42
basic$days.to.event.from.deliv = basic$days.to.event + 42
basic$fdate.from.deliv = basic$fdate_seg1_day + 42

# only want outcomes up to 1 year post delivery, which is day 323 from 42 days of delivery
basic <- basic[basic$days.to.event.from.deliv <= 365, ]

# Loss of enrollment within 1 year of delivery
basic$loss.enroll <- ifelse(basic$cat1_seg1==1 & basic$fdate.from.deliv <=365, 1, 0)
label(basic$loss.enroll) <- "Loss of enrollment within 1-year of delivery"

basic$outcome.num <- ifelse(is.na(basic$outcome), 0, 1)


basic$final.outcome <- ifelse(basic$outcome =="Death", 1,
                              ifelse(basic$outcome=="Overdose", 2, 
                                     ifelse(basic$loss.enroll==1 & is.na(basic$outcome), 3, NA)))
basic$final.outcome <- ifelse(is.na(basic$final.outcome) & basic$loss.enroll==1, 3, basic$final.outcome)
basic$final.outcome <- ifelse(is.na(basic$final.outcome) & basic$loss.enroll!=1, 4, basic$final.outcome)



# check: basic[which(basic$death.new==1), c("studyid", "days.to.event", "follow.up.days", "loss.enroll",  "death_day", "fdate_seg1_day", "death.or.od", "death.new", "final.outcome") ]
# some people have fdate that are greater than death_day (ex: ID "002007049737-02925139-05210809") as well as those with fdate before death date (maybe death recorded later on after disenrollment?)

basic$final.outcome.f <- factor(basic$final.outcome, levels=1:4, labels=c("Death", "Overdose", "Loss of enrollment", "Censored"))
label(basic$final.outcome.f) <- "Final outcome"


# follow-up times calculated without accounting for OD after disenrollment 
af_out365$days.to.event2 <- pmin(af_out365$overdose_day, af_out365$fdate_seg1_day, af_out365$momdeath_day)
af_out365$days.to.event2 <- ifelse(is.na(af_out365$overdose_day), pmin(af_out365$fdate_seg1_day, af_out365$momdeath_day), af_out365$days.to.event2)
af_out365$days.to.event2 <- ifelse(is.na(af_out365$momdeath_day), pmin(af_out365$fdate_seg1_day, af_out365$overdose_day), af_out365$days.to.event2)

af_out365$outcome2 <- ifelse(af_out365$days.to.event2== af_out365$momdeath_day , "Death", 
                             ifelse(af_out365$days.to.event2 == af_out365$overdose_day, "Overdose", NA))
af_out365$outcome2 <- ifelse(is.na(af_out365$outcome2) & is.na(af_out365$momdeath_day) & !is.na(af_out365$overdose_day), "Overdose", af_out365$outcome2)
af_out365$outcome2 <- ifelse(af_out365$days.to.event2 == af_out365$fdate_seg1_day, "Loss of enrollment", af_out365$outcome2)
# check 
# af_out365[which(af_out365$days.to.event != af_out365$days.to.event2), ]

basic$outcome2 <- af_out365[match(basic$studyid, af_out365$studyid), "outcome2"]
basic$days.to.event2 <- af_out365[match(basic$studyid, af_out365$studyid), "days.to.event2"]
basic$days.to.event2 <- ifelse(is.na(basic$days.to.event2), 365, basic$days.to.event2) # if missing then no event; change to max follow-up
basic$days.to.event2 <- ifelse(basic$days.to.event2 > 323, 323, basic$days.to.event) #follow up anchored to day 41 postpartum so event should be 323 days post partum

basic$final.outcome2 <- ifelse(basic$outcome2 =="Death", 1,
                               ifelse(basic$outcome2=="Overdose", 2, 
                                      ifelse(basic$loss.enroll==1 & is.na(basic$outcome2), 3, NA)))
basic$final.outcome2 <- ifelse(is.na(basic$final.outcome2) & basic$loss.enroll==1, 3, basic$final.outcome2)
basic$final.outcome2 <- ifelse(is.na(basic$final.outcome2) & basic$loss.enroll!=1, 4, basic$final.outcome2)

basic$final.outcome2.f <- factor(basic$final.outcome2, levels=1:4, labels=c("Death", "Overdose", "Loss of enrollment", "Censored"))
label(basic$final.outcome2.f) <- "Final outcome (w/o all payer data)"

# to check
# af_out365[which(!is.na(af_out365$momdeath_day) & !is.na(af_out365$overdose_day)), c("studyid", "overdose_day", "fdate_seg1_day", "momdeath_day", "days.to.event", "outcome", "outcome2", "days.to.event2")]

# Additional covariates to be adjusted for
basic$delivery_year <- as.integer(basic$delivery_year)
label(basic$delivery_year) <- "Year of delivery"

basic$delivery.year.cat <- as.factor(basic$delivery_year)
label(basic$delivery.year.cat) <- "Year of delivery as categorical"

basic$momraceeth_comb.new <- as.character(basic$momraceeth.new)
basic$momraceeth_comb.new <- ifelse(basic$momraceeth_comb.new %nin% c("Non-Hispanic White"), "All Other", basic$momraceeth_comb.new)
basic$momraceeth_comb.new <- as.factor(basic$momraceeth_comb.new)
label(basic$momraceeth_comb.new) <- "Race/Ethnicity binary"

label(basic$total.rx.days.b90_2) <- "Total days supplied MOUD before delivery (method 2)"
label(basic$total.rx.days.p41_2) <- "Total days supplied MOUD after delivery (method 2)"

label(basic$total.rx.days.b90_3) <- "Total days supplied MOUD before delivery (method 3)"
label(basic$total.rx.days.p41_3) <- "Total days supplied MOUD after delivery (method 3)"

basic$total.rx.days.exp <- rowSums(basic[ , c("total.rx.days.b90", "total.rx.days.p41")]) 
basic$total.rx.days.exp2 <- rowSums(basic[ , c("total.rx.days.b90_2", "total.rx.days.p41_2")]) 
basic$total.rx.days.exp3 <- extract_method3[match(basic$studyid, extract_method3$studyid), "total_unique_days_supplied"]

basic$moud.pre <- ifelse(basic$total.rx.days.b90 >0, 1, 0)
basic$moud.pre2 <- ifelse(basic$total.rx.days.b90_2 >0, 1, 0)
basic$moud.pre3 <- ifelse(basic$total.rx.days.b90_3 >0, 1, 0)

basic$moud.post <- ifelse(basic$total.rx.days.p41 >0, 1, 0)
basic$moud.post2 <- ifelse(basic$total.rx.days.p41_2 >0, 1, 0)
basic$moud.post3 <- ifelse(basic$total.rx.days.p41_3 >0, 1, 0)

basic$moud.both <- ifelse(basic$moud.pre==1 & basic$moud.post==1, 1, NA)
basic$moud.both <- ifelse(is.na(basic$moud.both) & (basic$moud.pre==1 | basic$moud.post==1), 2, basic$moud.both)
basic$moud.both <- ifelse(is.na(basic$moud.both), 0, basic$moud.both)

basic$moud.both2 <- ifelse(basic$moud.pre2==1 & basic$moud.post2==1, 1, NA)
basic$moud.both2 <- ifelse(is.na(basic$moud.both2) & (basic$moud.pre2==1 | basic$moud.post2==1), 2, basic$moud.both2)
basic$moud.both2 <- ifelse(is.na(basic$moud.both2), 0, basic$moud.both2)

basic$moud.both3 <- ifelse(basic$moud.pre3==1 & basic$moud.post3==1, 1, NA)
basic$moud.both3 <- ifelse(is.na(basic$moud.both3) & (basic$moud.pre3==1 | basic$moud.post3==1), 2, basic$moud.both3)
basic$moud.both3 <- ifelse(is.na(basic$moud.both3), 0, basic$moud.both3)

basic$moud.onlypre3 <- ifelse(basic$moud.pre3==1 & basic$moud.post3==0, 1, 0)
basic$moud.onlypost3 <- ifelse(basic$moud.post3==1 & basic$moud.pre3==0, 1, 0)

basic$moud.both3div <- ifelse(basic$moud.onlypre3==1, 1, 
                              ifelse(basic$moud.onlypost3==1, 2, 
                                     ifelse(basic$moud.both3==1, 3, 3)))
basic$moud.both3div <- ifelse(is.na(basic$moud.both3div), 0, basic$moud.both3div)
basic$moud.yes <- ifelse(basic$moud.both3div != 0, 1, 0)
basic$moud.yes.f <- factor(basic$moud.yes, levels=0:1, labels=c("No MOUD", "MOUD"))
label(basic$moud.yes.f) <- "MOUD Receipt"
basic$moud.both3div <- factor(basic$moud.both3div, levels=0:3, labels=c("No MOUD", "Only pre", "Only post", "Both pre and post"))

label(basic$moud.both3div) <- "MOUD timing (method 3)"

basic$moud.pre <- factor(basic$moud.pre, levels=0:1, labels=c("No MOUD", "MOUD"))
label(basic$moud.pre) <- "MOUD pre-delivery (method 1)"

basic$moud.pre2 <- factor(basic$moud.pre2, levels=0:1, labels=c("No MOUD", "MOUD"))
label(basic$moud.pre2) <- "MOUD pre-delivery (method 2)"

basic$moud.pre3 <- factor(basic$moud.pre3, levels=0:1, labels=c("No MOUD", "MOUD"))
label(basic$moud.pre3) <- "MOUD pre-delivery (method 3)"

basic$moud.post <- factor(basic$moud.post, levels=0:1, labels=c("No MOUD", "MOUD"))
label(basic$moud.post) <- "MOUD post-delivery (method 1)"

basic$moud.post2 <- factor(basic$moud.post2, levels=0:1, labels=c("No MOUD", "MOUD"))
label(basic$moud.post2) <- "MOUD post-delivery (method 2)"

basic$moud.post3 <- factor(basic$moud.post3, levels=0:1, labels=c("No MOUD", "MOUD"))
label(basic$moud.post3) <- "MOUD post-delivery (method 3)"

# 
basic$moud.both <- factor(basic$moud.both, levels=0:2, labels=c("No MOUD", "Both pre and post", "Either pre or post"))
label(basic$moud.both) <- "MOUD pre and post (method 1)"

basic$moud.both2 <- factor(basic$moud.both2, levels=0:2, labels=c("No MOUD", "Both pre and post", "Either pre or post"))
label(basic$moud.both2) <- "MOUD pre and post (method 2)"

basic$moud.both3 <- factor(basic$moud.both3, levels=0:2, labels=c("No MOUD", "Both pre and post", "Either pre or post"))
label(basic$moud.both3) <- "MOUD pre and post (method 3)"


label(basic$total.rx.days.exp2) <- "Total days supplied MOUD (method 2)"
label(basic$total.rx.days.exp3) <- "Total days supplied MOUD (method 3)"
# 

basic$anx.dep.comb <- ifelse(basic$anxiety==1 | basic$depression== 1, 1, 0)
basic$anx.dep.comb <- ifelse(is.na(basic$anx.dep.comb), 0, basic$anx.dep.comb)
basic$anx.dep.comb.f <- factor(basic$anx.dep.comb, levels=0:1, labels=c("No", "Yes"))
label(basic$anx.dep.comb.f) <- "Anxiety and depression combined"

basic$outcome.bin <- ifelse(basic$final.outcome %in% c(1, 2), 1, 0)

label(basic$total.rx.days.exp) <- "Total days supplied MOUD"
label(basic$total.rx.days.total_3bup) <- "Total days supplied buprenorphine"
moud <- basic[basic$moud.yes==1, ]

n_step6       <- nrow(moud)                                                                                                                                  
n_excl_nomoud <- n_step4 - n_step6     
cat("Final analytic sample (MOUD receipt):", formatC(n_step6, format = "d", big.mark = ","),
    "| Excluded (no MOUD receipt):", n_excl_nomoud, "\n")
# Final analytic sample (MOUD receipt): 9,419 | Excluded (no MOUD receipt): 7802 

# AL: 
table(moud$rxtype, useNA = "always")
#BUP only: 9084
#Naltrexone only: 232
#Both: 103




# AL: Create OPIOID PAIN FLAG; for PS; final MOUD cohort only 
# Any opioid for pain (not MOUD) in exposure window: -90 to +41 days                                                                                                         

library(readxl)                                                                                                                                                              

# CDC opioid list: all non-buprenorphine opioids only                                                                                                                        
cdc_opioids   <- read_excel("CDC_opioid NDC_oral MME conversion_update_2020.xlsx", sheet = "Opioids")
cdc_non_bup   <- cdc_opioids[cdc_opioids$Drug != "Buprenorphine", ]
cdc_pain_ndcs <- unique(cdc_non_bup$NDC)

# Redbook methadone: all routes (oral + injectable + IV)
redbook_meth        <- read.csv("redbook_methadone.csv")
redbook_meth$ndcpad <- suppressWarnings(sprintf("%011.0f", as.numeric(redbook_meth$ndcnum)))
redbook_meth_ndcs   <- unique(redbook_meth$ndcpad[!is.na(suppressWarnings(as.numeric(redbook_meth$ndcnum)))])

# Buprenorphine pain NDC overlap check
# Includes all buprenorphine pain files:
# CDC (patches + Belbuca), Buprenex/Belbuca/Butrans csv, generic transdermal, HCl powder/solution
cdc_bup_pain  <- cdc_opioids[cdc_opioids$Drug == "Buprenorphine" &
                               (cdc_opioids$Master_Form == "Patch, Extended Release" |
                                  cdc_opioids$PRODNME == "BELBUCA" |
                                  cdc_opioids$PRODNME == "BUTRANS"),]
                             
bup_csv        <- read.csv("NDCs_2020_Buprenex_Belbuca_Butrans.csv", colClasses = c(NDCNUM = "character"))
bup_trans      <- read_excel("ndcnum_bup_transdermal.xlsx")
bup_trans_ndcs <- sprintf("%011.0f", as.numeric(bup_trans$ndcnum))
bup_hcl        <- read_excel("ndcnum_BupHCl_powder_solution.xlsx")
bup_hcl_ndcs   <- sprintf("%011.0f", as.numeric(bup_hcl$ndcnum))
all_bup_pain_ndcs <- unique(c(as.character(cdc_bup_pain$NDC),
                              bup_csv$NDCNUM,
                              bup_trans_ndcs,
                              bup_hcl_ndcs))

# ndc is from NDCs_2020.csv — the MOUD NDC reference list loaded at top 
bup_overlap_ndcs <- all_bup_pain_ndcs[all_bup_pain_ndcs %in% ndc$ndcnum]
cat("Buprenorphine pain NDCs also in MOUD list (NDCs_2020.csv):",
    length(bup_overlap_ndcs), "of", length(all_bup_pain_ndcs), "\n")

# Standardize NDCs 
# sprintf uses float format — handles all 11-digit NDCs correctly
# formatC format="d" uses integer format — which maxes out above 2,147,483,647
# causing wrong values and false matches
# The below code leaves HCPCS codes unchanged but standardizes NDC codes via "grepl"
af_ext_df         <- as.data.frame(af_extract)
af_ext_df$ndc_std <- ifelse(grepl("^[0-9]+$", af_ext_df$ndc),
                            suppressWarnings(sprintf("%011.0f", as.numeric(af_ext_df$ndc))),
                            af_ext_df$ndc)

af_exp_moud <- af_ext_df[af_ext_df$BTrack == 90 &
                           (af_ext_df$Rx_B90_DOB == 1 | af_ext_df$Rx_DOB_P41 == 1) &
                           af_ext_df$StudyID %in% moud$studyid, ]


bup_overlap_hits <- af_exp_moud[af_exp_moud$ndc_std %in% bup_overlap_ndcs, ]
bup_pain_hits    <- af_exp_moud[af_exp_moud$ndc_std %in% all_bup_pain_ndcs, ]
cat("MOUD pregnancies with any buprenorphine pain NDC fill:",
    length(unique(bup_pain_hits$StudyID)), "\n")
cat("  Of those, from overlapping NDCs (ambiguous):",
    length(unique(bup_overlap_hits$StudyID)), "\n")
# 3 MOUD pregnancies had any buprenorphine pain NDC fill:
#   - 2 filled NDCs also in the MOUD list (NDCs_2020.csv) — use is ambiguous
#     (could be pain or OUD treatment), so excluded from pain flag
#   - 1 filled an NDC only on the pain list, not in the MOUD list —
#     captured in bup_pain_ndcs_clean below
# Note: bup_hcl has 52/52 NDCs in the MOUD list; the 2 ambiguous pregnancies
# come from bup_hcl fills specifically

# Non-overlapping buprenorphine pain NDCs — unambiguous pain fills
# (not in MOUD list); captures the 1 remaining pregnancy from the 3 above
bup_csv_clean       <- bup_csv$NDCNUM[!bup_csv$NDCNUM %in% ndc$ndcnum] #removes the NDCs that overlap with the MOUD list, keeps the remaining. 
bup_trans_clean     <- bup_trans_ndcs[!bup_trans_ndcs %in% ndc$ndcnum] #removes the NDCs that overlap with MOUD list, keeps the remaining. 
cdc_bup_pain_clean  <- as.character(cdc_bup_pain$NDC)[!as.character(cdc_bup_pain$NDC) %in% ndc$ndcnum] #CDC files, removes any that overlap with MOUD list, keeps the remaining.
bup_pain_ndcs_clean <- unique(c(bup_csv_clean, bup_trans_clean, cdc_bup_pain_clean)) #combines all three cleaned lists into one final set of unambiguous buprenorphien pain NDCs. 

# Final pain NDC list — includes non-overlapping buprenorphine pain NDCs
all_pain_ndcs <- unique(c(cdc_pain_ndcs, redbook_meth_ndcs, bup_pain_ndcs_clean))
cat("Total unique pain opioid NDCs:", length(all_pain_ndcs), "\n")


# Unique pregnancy IDs with pain opioid fills by period
pain_ids_overall <- unique(af_exp_moud$StudyID[af_exp_moud$ndc_std %in% all_pain_ndcs])
pain_ids_pre     <- unique(af_exp_moud$StudyID[af_exp_moud$ndc_std %in% all_pain_ndcs &
                                                 af_exp_moud$Rx_B90_DOB == 1])
pain_ids_post    <- unique(af_exp_moud$StudyID[af_exp_moud$ndc_std %in% all_pain_ndcs &
                                                 af_exp_moud$Rx_DOB_P41 == 1])

# Flags in moud
moud$pain_opioid_rx   <- as.integer(moud$studyid %in% pain_ids_overall)
label(moud$pain_opioid_rx) <- "Any pain opioid Rx in exposure window (not MOUD)"

moud$pain_opioid_pre  <- as.integer(moud$studyid %in% pain_ids_pre)
label(moud$pain_opioid_pre) <- "Pain opioid Rx pre-delivery (90-day window)"

moud$pain_opioid_post <- as.integer(moud$studyid %in% pain_ids_post)
label(moud$pain_opioid_post) <- "Pain opioid Rx post-delivery (0 to +41 days)"

moud$pain_opioid_timing <- ifelse(
  moud$pain_opioid_pre == 1 & moud$pain_opioid_post == 1, "Both",
  ifelse(moud$pain_opioid_pre  == 1, "Pre-delivery only",
         ifelse(moud$pain_opioid_post == 1, "Post-delivery only", "None")))
moud$pain_opioid_timing <- factor(moud$pain_opioid_timing,
                                  levels = c("None", "Pre-delivery only",
                                             "Post-delivery only", "Both"))
label(moud$pain_opioid_timing) <- "Pain opioid Rx timing relative to delivery"

# Summary
pct <- function(x) paste0(round(mean(x, na.rm = TRUE) * 100, 1), "%")
cat("\n Pain opioid Rx — MOUD cohort (n =", nrow(moud), ") \n")
cat("Any pain opioid:    ", sum(moud$pain_opioid_rx),
    paste0("(", pct(moud$pain_opioid_rx), ")\n"))
cat("Pre-delivery only:  ", sum(moud$pain_opioid_timing == "Pre-delivery only"), "\n")
cat("Post-delivery only: ", sum(moud$pain_opioid_timing == "Post-delivery only"), "\n")
cat("Both:               ", sum(moud$pain_opioid_timing == "Both"), "\n")
cat("None:               ", sum(moud$pain_opioid_timing == "None"), "\n")








# TEST CODE re. 1,983 pregnancies that were added: 

# Breakdown of the 1,983 pregnancies added via MOUD pharmacy fill only                                                                                                       
# (Expose1==0 but had an MOUD fill in extract — not captured by diagnosis alone)                                                                                           
added_ids <- basic$studyid[basic$expose1 == 0 & basic$filledrx == 1]                                                                                                         
cat("Pregnancies added via MOUD pharmacy fill only:", length(added_ids), "\n")

# Their MOUD fills from extract                                                                                                                                            
extract_added <- extract[extract$StudyID %in% added_ids, ]

# For HCPCS codes, gennme is NA after the ndc merge — label them by code
extract_added$med_type <- ifelse(!is.na(extract_added$gennme), extract_added$gennme,
                                 ifelse(extract_added$ndc %in% c("J0571","J0572","J0573",
                                                                 "J0574","J0575","G2068",
                                                                 "G2069","G2070","G2072",
                                                                 "G2079"), "Buprenorphine (HCPCS)",
                                        ifelse(extract_added$ndc %in% c("J2315"), "Naltrexone (HCPCS)",
                                               ifelse(extract_added$ndc %in% c("Q9991","Q9992"), "Buprenorphine (HCPCS)",
                                                      ifelse(extract_added$ndc %in% c("S0109"), "Buprenorphine (HCPCS)",
                                                             "Other")))))

# Summarise per pregnancy — which medication types did they fill?
med_summ <- extract_added %>%
  group_by(StudyID) %>%
  summarise(
    bup_mono   = as.integer(any(grepl("Buprenorphine", med_type, ignore.case = TRUE) &
                                  !grepl("Naloxone",   med_type, ignore.case = TRUE))),
    bup_nal    = as.integer(any(grepl("Buprenorphine", med_type, ignore.case = TRUE) &
                                  grepl("Naloxone",   med_type, ignore.case = TRUE))),
    naltrexone = as.integer(any(grepl("Naltrexone",   med_type, ignore.case = TRUE))),
    .groups = "drop"
  ) %>%
  mutate(rx_type = case_when(
    bup_nal  == 1 & naltrexone == 0                        ~ "Buprenorphine/Naloxone only",
    bup_mono == 1 & bup_nal == 0 & naltrexone == 0         ~ "Buprenorphine mono only",
    naltrexone == 1 & bup_mono == 0 & bup_nal == 0         ~ "Naltrexone only",
    TRUE                                                   ~ "Multiple medication types"
  ))

cat("\nMedication breakdown — mutually exclusive:\n")
print(table(med_summ$rx_type))

cat("\nNaltrexone StudyIDs:\n")
print(med_summ$StudyID[grepl("Naltrexone", med_summ$rx_type)])

# Were any of the 1,983 added using NDCs that also appear in the buprenorphine pain list?
# i.e. added to MOUD cohort via an NDC that is ambiguous (pain vs OUD)
extract_added$ndc_std <- ifelse(grepl("^[0-9]+$", extract_added$ndc),
                                suppressWarnings(sprintf("%011.0f", as.numeric(extract_added$ndc))),
                                extract_added$ndc)

bup_overlap_in_added <- extract_added[extract_added$ndc_std %in% bup_overlap_ndcs, ]
cat("\nOf the", length(added_ids), "added pregnancies:\n")
cat("Added using NDCs overlapping with buprenorphine pain list:",
    length(unique(bup_overlap_in_added$StudyID)), "\n")
if (length(unique(bup_overlap_in_added$StudyID)) > 0) {
  cat("Their StudyIDs:\n")
  print(unique(bup_overlap_in_added$StudyID))
}


# Product names for buprenorphine in the pharmacy-only added pregnancies
added_pharma_ids <- moud$studyid[moud$expose1 == 0 & moud$filledrx == 1]

# extract is already BTrack==90 and has prodnme merged from ndc
added_fills_exp <- extract[
  tolower(extract$StudyID) %in% added_pharma_ids &
    (extract$Rx_B90_DOB == 1 | extract$Rx_DOB_P41 == 1), ]

bup_fills <- added_fills_exp[
  !is.na(added_fills_exp$gennme) &
    grepl("buprenorphine", added_fills_exp$gennme, ignore.case = TRUE), ]

cat("Product names for buprenorphine MOUD fills in added pregnancies:\n")
print(sort(table(bup_fills$prodnme), decreasing = TRUE))







## MORE TEST CODE TO CONFIRM CHAD'S ERROR

############################################################
# Added pregnancies (Expose1==0 + pharmacy fill)
############################################################

added_ids <- moud$studyid[moud$expose1 == 0 & moud$filledrx == 1]                                                                                                               

# Pull exposure-window fills for added pregnancies only                                                                                                                           
added_fills <- extract[                                                                                                                                                           
  tolower(extract$StudyID) %in% added_ids &                                                                                                                                       
    (extract$Rx_B90_DOB == 1 | extract$Rx_DOB_P41 == 1), ]                                                                                                                        

bup_added_fills <- added_fills[
  !is.na(added_fills$gennme) &
    grepl("buprenorphine", added_fills$gennme, ignore.case = TRUE), ]
nal_added_fills <- added_fills[
  !is.na(added_fills$gennme) &
    grepl("naltrexone", added_fills$gennme, ignore.case = TRUE), ]

bup_added_ids <- unique(tolower(bup_added_fills$StudyID))
nal_added_ids <- unique(tolower(nal_added_fills$StudyID))

############################################################
# (1) BUP fills: confirm no OUD diagnosis
############################################################

cat("=== (1) BUP fills: OUD diagnosis check ===\n")
cat("Added pregnancies with BUP fills:", length(bup_added_ids), "\n")
cat("Of those, with Expose1==1 (OUD diagnosis per SAS pipeline):",
    sum(moud$expose1[moud$studyid %in% bup_added_ids] == 1), "\n")
# Expected: 0 — expose1==0 is the definition of the 'added' group,
# meaning the SAS pipeline found no OUD dx for these individuals.
# The reason they were missed was the p5 bug (BUP_TWindow not stacked),
# not absence of the SAS OUD lookup — so this confirms they have no OUD dx.

############################################################
# (2) Naltrexone: NAL-only vs concurrent with BUP
############################################################
cat("\n=== (2) Naltrexone: concurrent BUP? ===\n")
cat("Added pregnancies with NAL fills:", length(nal_added_ids), "\n")

nal_with_bup <- nal_added_ids[nal_added_ids %in% bup_added_ids]
nal_only     <- nal_added_ids[!nal_added_ids %in% bup_added_ids]

cat("  NAL + concurrent BUP in exposure window:", length(nal_with_bup), "\n")
cat("  NAL only (no BUP fills in window):",       length(nal_only),     "\n")

# Show product names for NAL-only group
if (length(nal_only) > 0) {
  nal_only_fills <- nal_added_fills[tolower(nal_added_fills$StudyID) %in% nal_only, ]
  cat("  NAL-only product names:\n")
  print(sort(table(nal_only_fills$prodnme), decreasing = TRUE))
}

############################################################
# (3) Naltrexone: OUD dx, AUD dx, or both
############################################################
cat("\n=== (3) Naltrexone: OUD and/or AUD diagnosis ===\n")
nal_moud <- moud[moud$studyid %in% nal_added_ids, ]

cat("NAL added pregnancies — OUD diagnosis (Expose1==1):",
    sum(nal_moud$expose1 == 1, na.rm = TRUE),
    "— expected 0 (added group is Expose1==0 by definition)\n")

cat("NAL added pregnancies — AUD diagnosis (alc_usedis==1):",
    sum(nal_moud$alc_usedis == 1, na.rm = TRUE), "\n")

cat("\nCross-tab OUD x AUD for NAL added pregnancies:\n")
print(table(
  OUD = ifelse(nal_moud$expose1 == 1, "Yes", "No"),
  AUD = ifelse(!is.na(nal_moud$alc_usedis) & nal_moud$alc_usedis == 1, "Yes", "No"),
  useNA = "ifany"
))
# All 60 naltrexone pregnancies have OUD = No
# 55 have neither OUD nor AUD diagnoses 
# 5 have an AUD diagnosis but no OUD diagnosis

## PLAN: (1) Drop the 5 with AUD dx; Keep the 55 without a dx but do a SA by dropping them 




############################################################################                                                                                                                    
# Naltrexone: concurrent BUP — full moud cohort                                                                                                                                   
# How many NAL users also have BUP fills in the exposure window?                                                                                                                  
# The SAS NAL_BUP_MET logic would have captured NAL + concurrent BUP.
############################################################################

# All exposure-window fills for the full moud cohort
moud_fills_all <- extract[
  tolower(extract$StudyID) %in% moud$studyid &
    (extract$Rx_B90_DOB == 1 | extract$Rx_DOB_P41 == 1), ]

# NAL and BUP fill IDs across full cohort
nal_all_ids <- unique(tolower(moud_fills_all$StudyID[
  !is.na(moud_fills_all$gennme) &
    grepl("naltrexone", moud_fills_all$gennme, ignore.case = TRUE)]))

bup_all_ids <- unique(tolower(moud_fills_all$StudyID[
  !is.na(moud_fills_all$gennme) &
    grepl("buprenorphine", moud_fills_all$gennme, ignore.case = TRUE)]))

nal_all_with_bup <- nal_all_ids[nal_all_ids %in% bup_all_ids]
nal_all_only     <- nal_all_ids[!nal_all_ids %in% bup_all_ids]

cat("=== Naltrexone: concurrent BUP — full moud cohort (n =", nrow(moud), ") ===\n")
cat("Total NAL users in exposure window:", length(nal_all_ids), "\n")
cat("  NAL + concurrent BUP:", length(nal_all_with_bup), "\n")
cat("  NAL only (no BUP):  ", length(nal_all_only), "\n\n")

# Break out the NAL-only group: added vs. original (Expose1==1)
nal_only_added    <- nal_all_only[nal_all_only %in%
                                    moud$studyid[moud$expose1 == 0 & moud$filledrx == 1]]
nal_only_original <- nal_all_only[!nal_all_only %in% nal_only_added]

cat("  Of the NAL-only pregnancies:\n")
cat("    Added via pharmacy fill (Expose1==0):", length(nal_only_added), "\n")
cat("    Original cohort (Expose1==1):        ", length(nal_only_original), "\n")

cat("=== NAL-only: OUD diagnosis breakdown ===\n")
nal_only_orig_moud <- moud[moud$studyid %in% nal_only_original, ]

cat("NAL-only, original cohort (Expose1==1):", length(nal_only_original), "\n")
cat("  With OUD dx (Expose1==1):", sum(nal_only_orig_moud$expose1 == 1), "\n")
cat("  Without OUD dx (Expose1==0):", sum(nal_only_orig_moud$expose1 == 0), "\n")
# Expected: all 172 have Expose1==1 from OUD diagnosis since no BUP means
# NAL_BUP_MET could not have qualified them — OUD dx is the only path

cat("\nNAL-only, added via pharmacy fill (Expose1==0):", length(nal_only_added), "\n")
cat("  With AUD dx:", sum(moud$alc_usedis[moud$studyid %in% nal_only_added] == 1, na.rm = TRUE), "\n")
cat("  No OUD or AUD dx:", sum(moud$alc_usedis[moud$studyid %in% nal_only_added] == 0 |
                                 is.na(moud$alc_usedis[moud$studyid %in% nal_only_added])), "\n")


## 335 total NAL users
# 60 added + 172 original = 232 NAL only 
# + 103 = concurrent with BUP (picked up in original SAS cohort) (103 + 232 = 335 total)
# The 172 original was in the original cohort because they had an OUD DX. 


############################################################################
# Drop all 60 added naltrexone pregnancies: 
# 5 have AUD diagnosis only (no OUD dx) - clear misclassification 
# 55 have neither OUD or AUD diagnosis -- insufficient evidence for OUD
############################################################################

cat("\n=== Naltrexone added pregnancies: OUD/AUD breakdown ===\n")
cat("StudyIDs with AUD dx only (n = 5):\n")
print(nal_aud_only)

cat("\nStudyIDs with no OUD or AUD dx (n = 55):\n")
print(nal_sensitivity_drop_ids)

moud <- moud[!moud$studyid %in% nal_added_ids, ]

n_step7        <- nrow(moud)
n_excl_nal_aud <- n_step6 - n_step7
cat("\nAfter excluding all added naltrexone pregnancies (no clear OUD evidence):",
    formatC(n_step7, format = "d", big.mark = ","),
    "| Excluded:", n_excl_nal_aud, "\n")




# save moud data
write.csv(moud, file="moud_cohort_202601.csv", row.names = FALSE)
save(moud, basic,
     n_step1,
     n_step2, n_excl_twins,
     n_step3, n_excl_od_exposure,
     n_step6, n_excl_nomoud,
     n_step7, n_excl_nal_aud,  # naltrexone AUD-only exclusion
     all_pain_ndcs,            # pain NDC list
     nal_aud_only,             # 5 dropped: naltrexone + AUD dx only, no OUD dx
     nal_rest_drop_ids, # 55 to drop, didn't have either an OUD or AUD dx, not definitely OUD
     file="moud_cohort_202601-ashley.Rdata") #this file has moud, basic, and flow diagram counts
save.image(file="moud_cohort_all_tables_202601.Rdata") # this has everything, moud, basic, af_extract, af_cov, af_master, af_outcomes, etc.







