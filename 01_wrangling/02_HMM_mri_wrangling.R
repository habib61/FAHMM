#######################################################################################
###
### File name : 02_HMM_mri_wrangling.R
### Program developer: aardepi1, Piet Aarden based on work of all below:
###                    S. Gardiner (Oxford BDI, stephen.gardiner@ndm.ox.ac.uk)
###                    H. Ganjgahi (Oxford BDI, habib.ganjgahi@bdi.ox.ac.uk)
###                    S. Pendleton (Oxford BDI, samantha.pendleton@ndm.ox.ac.uk)
### Date: 06-03-2023
### Project/Trial code: NO.MS - HMM
### Description: Compiling the input data for the HMM modeling
### Application: R 4.1.0
### Libraries used: tidyverse_1.3.2
### Source table locations: /funstorage/NSGDD_MS_MRI/HMM/interim_tables
###                         /funstorage/NSGDD_MS_MRI/HMM/00_source_data/
### Input: 01_results.csv, SIENA_FINAL_NormNoGM.txt, SIENA_FINAL_NormNoGM_SEQ.txt,
###        SIENAX_FINAL_NormNoGM.txt, cro_Native_T2vol.txt
### Output: 
###
#######################################################################################

# Load of libraries -------------------------------------------------------

library(tidyverse)

# Functions ---------------------------------------------------------------

# Source function script
#source("/data/users/qs9f68/HMM/Piet/HMM/00a_source_scripts/00_HMM_wrangling_functions.R")
source("../00a_source_scripts/00_HMM_wrangling_functions.R")

# Load of data ------------------------------------------------------------

# Load in previously wrangled data
#df <- read_csv("/funstorage/NSGDD_MS_MRI/HMM/interim_tables/01_results.csv", show_col_types = F) %>% as.data.frame()

#df <- read_csv("/data/users/qs9f68/HMM/Piet/HMM/interim_tables/01_results.csv", col_types = cols('c','c','D','i','i','d','d','d','f','f','f','d','d','d','d','f','f','f','d','d','d','d','d','d')) %>% as.data.frame()

df <- read_csv("../interim_tables/01_results.csv", col_types = cols('c','c','D','i','i','d','d','d','f','f','f','d','d','d','d','f','f','f','d','d','d','d','d','d')) %>% as.data.frame()

sna_path = "/data/ms/processed/mri/IDP/IDP_v2/SIENA_FINAL_NormNoGM.txt"
snaseq_path = "/data/ms/processed/mri/IDP/IDP_v2/SIENA_FINAL_NormNoGM_SEQ.txt"
snax_path = "/data/ms/processed/mri/IDP/IDP_v2/SIENAX_FINAL_NormNoGM.txt"
cro_volt2_path = "/data/ms/processed/mri/IDP/IDP_v2/cro_Native_T2vol.txt"


# studies to include
studies <- c("CFTY720D2201", "CFTY720D2301","CFTY720D2302", "CFTY720D2309", "CFTY720D2312",
             "CFTY720D2306", "CBAF312A2304", "COMB157G2301", "COMB157G2302")

# Replace VOLT2 of NeuroRx to native space --------------------------------

t2vol <- read.table(cro_volt2_path, header = FALSE,sep = '',na.strings = "",fill=TRUE)
colnames(t2vol) = c('STUDY','USUBJID','VISIT','vox','VOLT2_native')

t2vol <- t2vol %>%
  mutate(USUBJID = str_replace_all(str_remove_all(USUBJID, 'sub-'), "x", "_")) %>%
  rowwise() %>%
  mutate(MRIDT = rev(str_split(VISIT, "x", simplify = T))[1]) %>%
  as.data.frame() %>%
  select(c(USUBJID, STUDY, MRIDT, VOLT2_native)) %>%
  mutate(MRIDT = as.Date(MRIDT, "%Y%m%d")) %>%
  arrange(USUBJID, MRIDT, desc(VOLT2_native)) %>%
  filter(!duplicated(.[c("USUBJID", "MRIDT")]))

# Merge onto endpoints dataset and replace values for patients in question
df <- df %>%
  merge(., t2vol, by = c("USUBJID", "STUDY", "MRIDT"), all.x = T) %>%
  mutate(VOLT2 = ifelse(STUDY %in% unique(t2vol$STUDY), VOLT2_native, VOLT2)) %>%
  select(-VOLT2_native)
  
# wrangle mri measures from bdi pipeline ----------------------------------

# Regular siena
siena <- read.table(sna_path, header = FALSE, sep = '', na.strings = "", 
                    fill = TRUE, quote = '') %>% 
  select(-'V7')
colnames(siena) <- c('STUDY','USUBJID','base','VISIT','DAY','pbvc_b')

siena <- siena %>% 
  mutate(pbvc_b = as.numeric(pbvc_b)) %>% 
  filter(DAY > 0)

# siena sequentially
siena_norm <- read.table(snaseq_path, header = FALSE, sep = '', na.strings = "", fill=TRUE, quote='') %>% 
  select(-'V7')
colnames(siena_norm) <- c('STUDY','USUBJID','base','VISIT','DAY','pbvc_seq')

siena_norm <- siena_norm %>%
  mutate(pbvc_seq = as.numeric(pbvc_seq))

# Combining both and calculating difference in regular siena from visit to following visit
siena_all <- siena %>% 
  merge(., siena_norm, by = c("USUBJID", "STUDY", "VISIT", "DAY"), all.x = T) %>%
  group_by(USUBJID) %>%
  arrange(DAY, .by_group = TRUE) %>%
  mutate(pbvc_diff = -lag(pbvc_b) + pbvc_b)

# Get "bad" timepoints
cutoff <- quantile(siena_all$pbvc_diff[siena_all$pbvc_diff > 0], na.rm = TRUE)[[4]]

siena_all$pbvc_b[siena_all$pbvc_diff > cutoff | siena_all$pbvc_b > 1] = NA
siena_all$pbvc_seq[siena_all$pbvc_diff > cutoff | siena_all$pbvc_seq > 1] = NA

siena_all$year <- siena_all$DAY/365
kk <- lme4::lmer(data= siena_all, pbvc_seq ~ year + (year | USUBJID), 
                 REML = FALSE, control = lme4::lmerControl(optimizer ="Nelder_Mead"))
row.has.na <- is.na(siena_all$pbvc_seq)

#get the standardise residuals and replace the outliers with the predicted value
siena_all$res               = NA
siena_all$res[!row.has.na]  = residuals(kk, type="pearson", scaled = TRUE)

siena_all$p                 = siena_all$pbvc_seq
siena_all$p[row.has.na]     = predict(kk, siena_all[row.has.na, ], allow.new.levels = TRUE)

siena_all <- siena_all %>%
  group_by(USUBJID) %>% 
  arrange(DAY, .by_group = TRUE) %>%
  mutate(pbvc_prd = cumsum(p))

siena_all$pbvc_b2 <- siena_all$pbvc_b
siena_all$pbvc_b2[is.na(siena_all$pbvc_b2)] <- siena_all$pbvc_prd[is.na(siena_all$pbvc_b2)]
siena_all <- siena_all %>% 
  group_by(USUBJID) %>%
  arrange(DAY, .by_group = TRUE) %>%
  select(c(USUBJID, STUDY, VISIT, pbvc_b2))

# SIENAX
sienax <- read.table(snax_path, header = FALSE, sep = '', na.strings = "", fill=TRUE, quote='')
colnames(sienax) <- c('STUDY','USUBJID','VISIT','DAY','NGM','GM','NWM','WM','NBV','BV','scale')

# Combine SIENAX and SIENA (both versions)
mri_bdi <- sienax %>% 
  mutate(NBV = as.numeric(NBV),
         NBV2 = NBV) %>%
  merge(., siena_all, by = c("USUBJID", "STUDY", "VISIT"), all.x = T)

# Get baseline value for each USUBJID
BLNBV <- mri_bdi %>%
  filter(DAY == 0, !is.na(NBV2)) %>%
  select(c(USUBJID, NBV2)) %>%
  arrange(USUBJID, desc(NBV2)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  rename(BLNBV = NBV2)

mri_bdi <- mri_bdi %>%
  arrange(USUBJID, DAY, desc(NBV2)) %>%
  filter(!duplicated(.[c("USUBJID", "DAY")])) %>%
  merge(., BLNBV, by = "USUBJID", all.x = T) %>%
  filter(!is.na(BLNBV)) %>%
  arrange(USUBJID, DAY) %>%
  mutate(NBV2 = ifelse(!is.na(pbvc_b2), (BLNBV * (pbvc_b2/100 + 1)), NBV2)) %>%
  select(-BLNBV)

# Remove extra's from conventional USUBJID and create the date column
mri_bdi <- mri_bdi %>%
  mutate(USUBJID = str_replace_all(str_remove_all(USUBJID, 'sub-'), "x", "_")) %>%
  rowwise() %>%
  mutate(MRIDT = rev(str_split(VISIT, "x", simplify = T))[1]) %>%
  as.data.frame() %>%
  select(c(USUBJID, STUDY, MRIDT, NBV2)) %>%
  mutate(MRIDT = as.Date(MRIDT, "%Y%m%d"))

# Merge BDI mri endpoints onto general endpoint table ---------------------

# Merge the two dataframes and remove patients with a abnormal small NBV
df_final <- df %>%
  merge(., mri_bdi, by = c("USUBJID", "STUDY", "MRIDT"), all.x = T) %>%
  filter(!(USUBJID %in% unique(mri_bdi[which(mri_bdi$NBV2 < 1150000), "USUBJID"]))) %>%
  arrange(USUBJID, MONTH) %>%
  select(-MRIDT)

# Save table
#write.csv(df_final, "/data/users/qs9f68/HMM/Piet/HMM/interim_tables/02_results.csv", row.names = F)
write.csv(df_final, "../interim_tables/02_results.csv", row.names = F)
