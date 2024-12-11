#######################################################################################
###
### File name : 01_HMM_endpoint_wrangling.R
### Program developer: aardepi1, Piet Aarden based on work of all below:
###                    S. Gardiner (Oxford BDI, stephen.gardiner@ndm.ox.ac.uk)
###                    H. Ganjgahi (Oxford BDI, habib.ganjgahi@bdi.ox.ac.uk)
###                    S. Pendleton (Oxford BDI, samantha.pendleton@ndm.ox.ac.uk)
### Date: 06-07-2023
### Project/Trial code: NO.MS - HMM
### Description: Compiling the input data for the HMM modeling
### Application: R 4.1.0
### Libraries used: tidyverse_1.3.2, haven_2.5.1, lubridate_1.8.0
### Source table locations: /vob/CNOMSANON/anon/anon_2/anon_source/
###                         /funstorage/NSGDD_MS_MRI/HMM/00_source_data/
### Input: 
### Output: 
###
#######################################################################################

# Load of libraries -------------------------------------------------------

library(tidyverse)
library(haven)
library(lubridate)

# Functions ---------------------------------------------------------------

# Source function script
#source("/data/users/qs9f68/HMM/Piet/HMM/00a_source_scripts/00_HMM_wrangling_functions.R")
source("00a_source_scripts/00_HMM_wrangling_functions.R")

# Load of data ------------------------------------------------------------

# define user

#source_saspath <- paste0("/view/", user, "_view/vob/CNOMSANON/anon/anon_2/anon_source/")
#source_saspath <- '/data/ms/unprocessed/clinical/NOMS_Version2_20211022-OXF-ANALYTICS/'
#base <- read_sas(paste0(source_saspath, "base.sas7bdat")) %>% as.data.frame()
#msfc <- read_sas(paste0(source_saspath, "msfc.sas7bdat")) %>% as.data.frame()
#edss <- read_sas(paste0(source_saspath, "edss.sas7bdat")) %>% as.data.frame()
#mri <- read_sas(paste0(source_saspath, "mri.sas7bdat")) %>% as.data.frame()
#sdmt <- read_sas(paste0(source_saspath, "sdmt.sas7bdat")) %>% as.data.frame()
#follow <- read_sas(paste0(source_saspath, "follow.sas7bdat")) %>% as.data.frame()
#cdwevents <- read_sas(paste0(source_saspath, "cdwevents.sas7bdat")) %>% as.data.frame()


base=read.csv('../data/base.csv')
msfc=read.csv('../data/fcs_pasat.csv') %>% rename(DATE=fcs_DATE, DAY=fcs_DAY, PASAT=PASAT3)
edss=read.csv('../data/edss.csv')
mri=read.csv('../data/mri.csv') %>% rename(DATE=mri_DATE, DAY=mri_DAY)
relapse = read.csv('../data/00_results.csv')
# Load in previously wrangled data
#relapse <- read_csv("/funstorage/NSGDD_MS_MRI/HMM/interim_tables/00_results.csv", show_col_types = F) %>% as.data.frame()

#relapse <- read_csv("/data/users/qs9f68/HMM/Piet/HMM/interim_tables/00_results.csv", col_types = cols('c','c','D','D','i','i','f','D','i','f','f','f','f')) %>% as.data.frame()

#relapse <- read_csv("../interim_tables/00_results.csv", col_types = cols('c','c','D','D','i','i','f','D','i','f','f','f','f')) %>% as.data.frame()

# studies to include
#studies <- c("CFTY720D2201", "CFTY720D2301","CFTY720D2302", "CFTY720D2309", "CFTY720D2312",
#             "CFTY720D2306", "CBAF312A2304", "COMB157G2301", "COMB157G2302")

# Concatenating to one baseline and attach follow-up visits ---------------

# edss_mod: table with EDSS data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst EDSS)
# - Only non-missing EDSS observations are kept
edss_mod <- edss %>%
  filter(DAY <= 1, !is.na(EDSS)) %>%
  select(c(USUBJID, STUDYID, DATE, DAY, EDSS)) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., edss[which((edss$DAY > 1) & (!is.na(edss$EDSS))), 
                    c("USUBJID", "STUDYID", "DATE", "DAY", "EDSS")]) %>%
  arrange(USUBJID, DATE, desc(EDSS)) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  mutate(DAY = ifelse(DAY <= 1, 1, DAY))

# df_month: Table with the month assigned (addition to df_wtrt)
# - Using Dieter Haering's visit window function we get the interval for each month
# - Assign the monthly visit using data.table foverlap
# - Assign baseline month (DAY == 1) to -1
mco_monthly <- vw(seq(0, 175, 1))$dat
edss_mod <- data.table::data.table(edss_mod)
edss_mod[, DAYTMP := DAY]
mcotmp <- data.table::data.table(mco_monthly[c("Month", "upper", "lower", "target")])
data.table::setkey(mcotmp, lower, upper)

edss_month <- data.table::foverlaps(edss_mod, mcotmp, by.x = c('DAY', 'DAYTMP'), 
                                  by.y = c('lower', 'upper')) %>% 
  as.data.frame() %>%
  rename(MONTH = Month) %>%
  mutate(MONTH = ifelse(DAY == 1, -1, MONTH)) %>%
  select(c(USUBJID, STUDYID, DATE, DAY, target, MONTH, EDSS))

# - Only one observation per month is kept, eliminate multiple observations
#   per month by keeping the once closest to the target date
edss_month <- edss_month %>%
  mutate(DIFF = abs(target-DAY)) %>%
  arrange(USUBJID, MONTH, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "MONTH")])) %>%
  select(-c(target, DIFF))

# t25fwm_mod: table with timed 25-foot walk data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst T25FWM)
# - Only non-missing T25FWM observations are kept
# - rename DATE column for later
t25fw_mod <- msfc %>%
  select(c(USUBJID, STUDYID, DATE, DAY, T25FWM))%>%
  filter(!is.na(DAY), DAY <= 1, !is.na(T25FWM)) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., msfc[which((msfc$DAY > 1) & (!is.na(msfc$T25FWM)) & (!is.na(msfc$DAY))), 
                    c("USUBJID", "STUDYID", "DATE", "DAY", "T25FWM")]) %>%
  arrange(USUBJID, DATE, desc(T25FWM)) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(T25DAY = DAY) %>%
  select(c(USUBJID, T25DAY, T25FWM)) %>%
  mutate(T25DAY = ifelse(T25DAY <= 1, 1, T25DAY))

# hpt9_mod: table with 9-hole peg test data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst HPT9M)
# - Only non-missing HPT9M observations are kept
# - rename DATE column for later
hpt9_mod <- msfc %>%
  select(c(USUBJID, STUDYID, DATE, DAY, HPT9M)) %>%
  filter(!is.na(DAY), DAY <= 1, !is.na(HPT9M)) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., msfc[which((msfc$DAY > 1) & (!is.na(msfc$HPT9M)) &  (!is.na(msfc$DAY))), 
                    c("USUBJID", "STUDYID", "DATE", "DAY", "HPT9M")]) %>%
  arrange(USUBJID, DATE, desc(HPT9M)) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(HPTDAY = DAY) %>%
  select(c(USUBJID, HPTDAY, HPT9M)) %>%
  mutate(HPTDAY = ifelse(HPTDAY <= 1, 1, HPTDAY))

# pasat_mod: table with PASAT data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date (keeping worst PASAT)
# - Only non-missing PASAT observations are kept
# - rename DATE column for later
pasat_mod <- msfc %>%
  select(c(USUBJID, STUDYID, DATE, DAY, PASAT)) %>%
  filter(!is.na(DAY), DAY <= 1, !is.na(PASAT)) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., msfc[which((msfc$DAY > 1) & (!is.na(msfc$PASAT)) & (!is.na(msfc$DAY))), 
                    c("USUBJID", "STUDYID", "DATE", "DAY", "PASAT")]) %>%
  arrange(USUBJID, DATE, PASAT) %>%
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(PASATDAY = DAY) %>%
  select(c(USUBJID, PASATDAY, PASAT)) %>%
  mutate(PASATDAY = ifelse(PASATDAY <= 1, 1, PASATDAY))

# mri_mod: table mri data
# - Take only one baseline observation (before FDD, DAY == 1)
# - Add remaining observations
# - Remove multiple observations on one date
# - rename DATE column for later
mri_mod <- mri %>%
  dplyr::select(c(USUBJID, STUDYID, DATE, DAY, NUMGDT1, VOLT2,NBV)) %>%
  filter(DAY <= 1, !is.na(DAY)) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., mri[which((mri$DAY > 1) & (mri$USUBJID != "") ), 
                   c("USUBJID", "STUDYID", "DATE", "DAY", "NUMGDT1", "VOLT2","NBV")]) %>%
  arrange(USUBJID, DATE) %>% 
  filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
  rename(MRIDT = DATE, MRIDAY = DAY) %>%
  select(c(USUBJID, MRIDT, MRIDAY, NUMGDT1, VOLT2, NBV)) %>%
  mutate(MRIDAY = ifelse(MRIDAY <= 1, 1, MRIDAY))

# t25fw_edss_map: table containg mapping of t25fw data to edss observations
# - Aligns all observations in t25fw_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference or 30 days
# - Makes sure that mapped edss is only used once
t25fw_edss_map <- t25fw_mod %>%
  select(c(USUBJID, T25DAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(T25DAY - DAY)) %>%
  filter(DIFF <= 30) %>%
  arrange(USUBJID, T25DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "T25DAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, T25DAY) %>%
  select(-DIFF)

# hpt9_edss_map: table containg mapping of hpt9 data to edss observations
# - Aligns all observations in hpt9_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference or 30 days
# - Makes sure that mapped edss is only used once
hpt9_edss_map <- hpt9_mod %>%
  select(c(USUBJID, HPTDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(HPTDAY - DAY)) %>%
  filter(DIFF <= 30) %>%
  arrange(USUBJID, HPTDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "HPTDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, HPTDAY) %>%
  select(-DIFF)

# pasat_edss_map: table containg mapping of pasat data to edss observations
# - Aligns all observations in pasat_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference or 30 days
# - Makes sure that mapped edss is only used once
pasat_edss_map <- pasat_mod %>%
  select(c(USUBJID, PASATDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(PASATDAY - DAY)) %>%
  filter(DIFF <= 30) %>%
  arrange(USUBJID, PASATDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "PASATDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, PASATDAY) %>%
  select(-DIFF)

# mri_edss_map: table containg mapping of MRI data to edss observations
# - Aligns all observations in mri_mod to a unique edss observation in edss_month
# - Maps observations to closest edss observations with a maximum absolute difference or 90 days
# - Makes sure that mapped edss is only used once
mri_edss_map <- mri_mod %>% 
  select(c(USUBJID, MRIDAY)) %>% 
  merge(., edss_month[c("USUBJID", "DAY")], by = "USUBJID", all.x = T) %>%
  mutate(DIFF = abs(MRIDAY - DAY)) %>%
  filter(DIFF <= 90) %>%
  arrange(USUBJID, MRIDAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "MRIDAY")])) %>%
  arrange(USUBJID, DAY, DIFF) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  arrange(USUBJID, MRIDAY) %>%
  select(-DIFF)

# cdwevents_mod: table containing the 6-month confirmed PIRA events
# - Subset cdwevents to 6 month confirmed PIRA events
cdwevents_mod <- cdwevents %>%
  filter(DISEVENT == "T6MCDW", PIRAFLG == "Yes", STUDY %in% studies)

# df: Table containing all endpoints and a flag for a PIRA onset
# - Use previous *_mod tables and *_edss_map tables to merge everything
# - Merge the PIRA flag on top of the observations using cdwevents_mod
# - Keep MRIDT as its used later for merging addition MRI data

df <- edss_month %>%
  merge(., t25fw_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., hpt9_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., pasat_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., mri_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
  merge(., t25fw_mod, by.x = c("USUBJID", "T25DAY"), all.x = T) %>%
  merge(., hpt9_mod, by = c("USUBJID", "HPTDAY"), all.x = T) %>%
  merge(., pasat_mod, by = c("USUBJID", "PASATDAY"), all.x = T) %>%
  merge(., mri_mod, by = c("USUBJID", "MRIDAY"), all.x = T) %>%
  select(c(USUBJID, STUDYID, MRIDT, DATE, DAY, MONTH, EDSS, T25FWM, HPT9M, PASAT, NUMGDT1, VOLT2,NBV)) %>%
  arrange(USUBJID, DATE)


#df <- edss_month %>%
#  merge(., t25fw_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
#  merge(., hpt9_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
#  merge(., pasat_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
#  merge(., mri_edss_map, by = c("USUBJID", "DAY"), all.x = T) %>%
#  merge(., t25fw_mod, by.x = c("USUBJID", "T25DAY"), all.x = T) %>%
#  merge(., hpt9_mod, by = c("USUBJID", "HPTDAY"), all.x = T) %>%
#  merge(., pasat_mod, by = c("USUBJID", "PASATDAY"), all.x = T) %>%
#  merge(., mri_mod, by = c("USUBJID", "MRIDAY"), all.x = T) %>%
#  select(c(USUBJID, STUDY, MRIDT, DATE, DAY, MONTH, EDSS, T25FWM, HPT9M, PASAT, NUMGDT1, VOLT2)) %>%
#  merge(., cdwevents_mod[c("USUBJID", "ONSTDAY", "PIRAFLG")], by.x = c("USUBJID", "DAY"), by.y = c("USUBJID", "ONSTDAY"), all.x = T) %>%
#  mutate(PIRAFLG = ifelse(is.na(PIRAFLG), 0, 1)) %>%
#  arrange(USUBJID, DATE)

# df_wrel: Table with the relapse info for each observation (addition to df)
# - use AddRelapseInfo to check whether observation was measured during a relapse
df_wrel <- AddRelapseInfo(df, relapse)

# df_wtrt: Table containing the treatment at time of observation (addition to df_wrel)
# - use AddTRTInfo to check what treatment patient was on during observation
df_wtrt <- AddTRTInfo(df_wrel, follow)

# df_final: final table which includes demographic data
# - merge demographic data on table
# - transform ARM to ACTVTRT (active treatment 0/1)
# - offset AGE and DURFS longitudinally (AGELG, DURFSLG respectively)
# - Impute a longitudinal MS phenotype (MSTYPELG) by the condition of someone having
#   a PIRA event and a EDSS >= 3 (captured in df_spmsonset constructed below)
df_spmsonset <- df_wtrt %>%
  filter(PIRAFLG == 1, EDSS >= 3) %>%
  arrange(USUBJID, MONTH) %>%
  filter(!duplicated(.$USUBJID)) %>%
  rename(PIRAMONTH = MONTH) %>%
  select(c(USUBJID, PIRAMONTH))

df_final <- df_wtrt %>%
  merge(., base[c("USUBJID", "SEX", "AGE", "MSTYPE", "DURFS", "RELPST1Y", "RELPST2Y")], by = "USUBJID", all.x = T) %>%
  merge(., df_spmsonset, by = c("USUBJID"), all.x = T) %>%
  mutate(ACTVTRT = ifelse(ARM %in% c("No treatment", "Placebo"), 0, 1),
         YEARS = DAY/365.25,
         DURFS = plyr::mapvalues(DURFS, c("[0,2)", "[2,5)", "[5,10)", "[10,30)", "[30,50)"), c("1", "3.5", "7.5", "20", "40")),
         DURFS = as.numeric(DURFS),
         DURFSLG = DURFS + YEARS,
         AGELG = floor(AGE + YEARS),
         MSTYPELG = case_when(
           MSTYPE == "SPMS" ~ "SPMS",
           MSTYPE == "PPMS" ~ "PPMS",
           (MSTYPE == "RRMS") & (is.na(PIRAMONTH)) ~ "RRMS",
           (MSTYPE == "RRMS") & (!is.na(PIRAMONTH)) & (MONTH < PIRAMONTH) ~ "RRMS",
           (MSTYPE == "RRMS") & (!is.na(PIRAMONTH)) & (MONTH >= PIRAMONTH) ~ "SPMS")
         ) %>%
  arrange(USUBJID, MONTH) %>%
  select(c(USUBJID, STUDY, MRIDT, MONTH, DAY, YEARS,
           AGE, AGELG, SEX, MSTYPE, MSTYPELG, DURFS, DURFSLG, 
           RELPST1Y, RELPST2Y, 
           PIRAFLG, RELAPSE, ACTVTRT,
           EDSS, T25FWM, HPT9M, PASAT,
           NUMGDT1, VOLT2))

# Save table
#write.csv(df_final, "/data/users/qs9f68/HMM/Piet/HMM/interim_tables/01_results.csv", row.names = F)
write.csv(df_final, "../interim_tables/01_results.csv", row.names = F)
