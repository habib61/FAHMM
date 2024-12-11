#######################################################################################
###
### File name : 00_HMM_relapse_wrangling.R
### Program developer: aardepi1, Piet Aarden based on work of all below:
###                    S. Gardiner (Oxford BDI, stephen.gardiner@ndm.ox.ac.uk)
###                    H. Ganjgahi (Oxford BDI, habib.ganjgahi@bdi.ox.ac.uk)
###                    S. Pendleton (Oxford BDI, samantha.pendleton@ndm.ox.ac.uk)
### Date: 18-04-2023
### Project/Trial code: NO.MS - HMM
### Description: Compiling the relapse input data for the HMM modeling
### Application: R 4.1.0
### Libraries used: tidyverse_1.3.2, haven_2.5.1, reshape2_1.4.4
### Source table locations: /vob/CNOMSANON/anon/anon_2/anon_source/
### Input: base.sas7bdat, relapse.sas7bdat
### Output: 00_results.csv
###
#######################################################################################

# Load of libraries -------------------------------------------------------

library(tidyverse)
library(haven)
library(reshape2)

# Functions ---------------------------------------------------------------

# Source function script
#source("/data/users/qs9f68/HMM/Piet/HMM/00a_source_scripts/00_HMM_wrangling_functions.R")
source("../00a_source_scripts/00_HMM_wrangling_functions.R")

# Load of data ------------------------------------------------------------

# define user
user <- system("whoami", intern = TRUE)
#source_saspath <- paste0("/view/", user, "_view/vob/CNOMSANON/anon/anon_2/anon_source/")
source_saspath <- '/data/ms/unprocessed/clinical/NOMS_Version2_20211022-OXF-ANALYTICS/'
base <- read_sas(paste0(source_saspath, "base.sas7bdat")) %>% as.data.frame()
relapse <- read_sas(paste0(source_saspath, "relapse.sas7bdat")) %>% as.data.frame()

# studies to include
studies <- c("CFTY720D2201", "CFTY720D2301","CFTY720D2302", "CFTY720D2309", "CFTY720D2312",
             "CFTY720D2306", "CBAF312A2304", "COMB157G2301", "COMB157G2302")


# Wrangle and clean up relapse table --------------------------------------

# 1. Some relapses have neither start nor end date and are not confirmed (i.e. no confirmation date)
#    These can be removed.
rel_mod <- relapse %>%
  filter(STUDY %in% studies, !(is.na(STRTDATE) & is.na(ENDDATE) & is.na(CONFDATE)))

# 2. Some patients have multiple rows for the same relapse start date but with information spread out over the rows.
#    Concatenate these rows as much as possible. Following arranging and filtering takes care of it.
rel_mod <- rel_mod %>%
  arrange(USUBJID, STRTDATE, desc(ENDDATE), CONFDATE, desc(STEROID)) %>% 
  filter(!duplicated(.[c("USUBJID", "STRTDATE")]))

# 3. Some patients have multiple rows for the same relapse end date but with information spread out over the rows.
#    Concatenate these rows as much as possible. Following arranging and filtering takes care of it.
rel_mod <- rel_mod %>%
  filter(!is.na(ENDDATE)) %>%
  arrange(USUBJID, ENDDATE, CONFDATE, STRTDATE) %>%
  filter(!duplicated(.[c("USUBJID", "ENDDATE")])) %>%
  bind_rows(rel_mod[which(is.na(rel_mod$ENDDATE)),]) %>%
  arrange(USUBJID, STRTDATE, ENDDATE)

# 4. Some relapse have neither start nor end date but do have a confirmation date
#    Check whether for these, another row with the exact same confirmation date exists with non-missing
#    start and end date (this is due to partial dates which have not gone through imputation, all CBAF312A2304). 
#    Will impute them with 10 days before and after relapse
rel_mod <- rel_mod %>%
  filter(is.na(STRTDATE) & is.na(ENDDATE) & !is.na(CONFDATE)) %>%
  mutate(STRTDATE = CONFDATE - 10,
         ENDDATE = CONFDATE + 10) %>%
  bind_rows(rel_mod[which(!(is.na(rel_mod$STRTDATE) & is.na(rel_mod$ENDDATE) & !is.na(rel_mod$CONFDATE))),]) %>%
  arrange(USUBJID, STRTDATE, ENDDATE)

# 5. Duration of relapses should be truncated at 90 days (checked!). Missing end dates are imputed according to the
#    90 days rule.
rel_mod <- rel_mod %>%
  mutate(ENDDATE = ifelse(is.na(ENDDATE), STRTDATE + 89, ENDDATE),
         ENDDATE = as.Date(ENDDATE, origin = "1970-01-01"))

# 6. Some relapses are overlapping. For those relapse make them adjacent.
for(i in seq(1, nrow(rel_mod) - 1)){
 rel_mod[i, "ENDDATE"] = case_when(
   (rel_mod[i + 1, "USUBJID"] != rel_mod[i, "USUBJID"]) ~ rel_mod[i, "ENDDATE"],
   (rel_mod[i + 1, "USUBJID"] == rel_mod[i, "USUBJID"]) & (rel_mod[i, "ENDDATE"] < rel_mod[i + 1, "STRTDATE"]) ~ rel_mod[i, "ENDDATE"],
   (rel_mod[i + 1, "USUBJID"] == rel_mod[i, "USUBJID"]) & (rel_mod[i, "ENDDATE"] >= rel_mod[i + 1, "STRTDATE"]) ~ rel_mod[i+1, "STRTDATE"] - 1
   )
}

# 7. With dates changed/imputed, we recalculate the DAY columns in this table
rel_mod <- calculateDays(rel_mod, base[c("USUBJID", "FDD")])

# Save relapse table as intermediate output to be used later
#write.csv(rel_mod,"/data/users/qs9f68/HMM/Piet/HMM/interim_tables/00_results.csv", row.names = F)
write.csv(rel_mod,"../interim_tables/00_results.csv", row.names = F)
