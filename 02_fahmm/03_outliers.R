####################################################################################
#############################baseline relapse 
####################################################################################
####################################################################################
library(tidyverse)
library(lme4)

#follow           = read.csv('/data/users/qs9f68/HMM/Piet/HMM/interim_tables/02_results.csv') 
follow           = read.csv('../interim_tables/02_results.csv') 
####below is for defining baseline relapse from next 3 months
relB             = follow%>%group_by(USUBJID)%>%filter(DAY<=100)%>%summarise(rel=sum(RELAPSE))%>%filter(rel>0)
#split the follow up table into two parts: subjects with baseline Relapse and not
SubRel           = follow%>%filter(USUBJID %in% relB$USUBJID)
SubRel$RELAPSE[SubRel$DAY==1]=1
SubNoRel         = follow%>%filter(!USUBJID %in% relB$USUBJID)
df               = rbind(SubRel,SubNoRel)

# Check for missing values at baseline, for complete missing values at follow-up visits, and for patients
# who only contributed one visit (baseline)
baseline <- df %>%
  rowwise() %>%
  mutate(MISSFULL = all(is.na(EDSS), is.na(T25FWM), is.na(HPT9M), is.na(PASAT), 
                        is.na(NUMGDT1), is.na(VOLT2), is.na(NBV2))) %>%
  filter(MISSFULL == FALSE) %>%
  as.data.frame() %>%
  group_by(USUBJID) %>%
  mutate("n" = n()) %>%
  as.data.frame() %>%
  filter(MONTH == -1) %>%
  rowwise() %>%
  mutate(MISSING = any(is.na(EDSS), is.na(T25FWM), is.na(HPT9M), is.na(NUMGDT1), is.na(VOLT2), is.na(NBV2))) %>%
  as.data.frame() %>%
  filter(n > 1, MISSING == F) %>%
  select(-c(MISSFULL, n, MISSING))

# include only the patients in table "baseline" and remove large gaps
df_mod <- df %>%
  filter(USUBJID %in% baseline$USUBJID)%>%
  mutate(MONTHN = MONTH,
         MONTHN = ifelse(MONTHN == -1, 0, MONTHN)) %>%
  group_by(USUBJID) %>% 
  arrange(MONTHN, .by_group = TRUE) %>% 
  mutate(deltaM = -lag(MONTHN) + MONTHN) %>% 
  ungroup() %>%
  mutate(deltaM = ifelse(is.na(deltaM), 0, deltaM),
         USUBJID = as.factor(USUBJID),
         VOLT2 = VOLT2^(1/3)) %>%
  filter(deltaM < 41)

baseline%>%group_by(STUDY)%>%summarise(N=length(USUBJID))
#write.csv(baseline,'/data/users/qs9f68/HMM/Piet/HMM/interim_tables/03_baseline_results.csv',row.names = FALSE)
write.csv(baseline,'../interim_tables/03_baseline_results.csv',row.names = FALSE)

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################## remove outliers; rebaseline large gaps and predict missings ###
################################################################################################
######## find and remove the outliers (HPT9M and T25FWM, brain) and replace them with marginal predictions using GAM 
######## the outliers are defined based on https://pubmed.ncbi.nlm.nih.gov/22052713/ paper    
######## which removing timepoints with standardised residulas more than 10 

#T25 outliers
ff           = lmer(T25FWM~YEARS+AGE+SEX+DURFS+RELPST1Y+MSTYPE+RELAPSE+ACTVTRT+(YEARS|USUBJID),data=df_mod,REML = FALSE,control = lmerControl(optimizer ="Nelder_Mead"))
#in below columns 6,8,11,13,19 are AGE, SEX, DURFS, RELPST1Y and T25FWM respectively  
row.has.na   = apply(df_mod[,c(6,8,11,13,19)], 1, function(x){any(is.na(x))})
#get the standardised residuals and replace the outliers with NA
df_mod$res = NA
df_mod$res[!row.has.na]  = residuals(ff,type="pearson", scaled=TRUE)
Ind=which(abs(df_mod$res)>10)
df_mod$T25FWM[Ind]=NA

#HPT9M outliers
ff           = lmer(HPT9M~YEARS+AGE+SEX+DURFS+RELPST1Y+MSTYPE+RELAPSE+ACTVTRT+(YEARS|USUBJID),data=df_mod,REML = FALSE)
#in below columns 6,8,11,13,20 are AGE, SEX, DURFS, RELPST1Y and HPT9M respectively  
row.has.na   = apply(df_mod[,c(6,8,11,13,20)], 1, function(x){any(is.na(x))})
#get the standardised residuals and replace the outliers with NA
df_mod$res = NA
df_mod$res[!row.has.na]  = residuals(ff,type="pearson", scaled=TRUE)
Ind=which(abs(df_mod$res)>10)
df_mod$HPT9M[Ind]=NA

#T2 Volume: Outliers
ff           = lmer(VOLT2~YEARS+AGE+SEX+DURFS+RELPST1Y+MSTYPE+RELAPSE+ACTVTRT+(YEARS|USUBJID),data=df_mod,REML = FALSE,control = lmerControl(optimizer ="Nelder_Mead"))
#in below columns 6,8,11,13,23 are AGE, SEX, DURFS, RELPST1Y and VOLT2 respectively  
row.has.na   = apply(df_mod[,c(6,8,11,13,23)], 1, function(x){any(is.na(x))})
#get the standardise residuals and replace the outliers with the predicted value
df_mod$res=NA
df_mod$res[!row.has.na]  = residuals(ff,type="pearson", scaled=TRUE)
#follow$pred[!row.has.na] = predict(ff)
Ind=which(abs(df_mod$res)>10)
df_mod$VOLT2[Ind]=NA


#write.csv(df_mod,'/data/users/qs9f68/HMM/Piet/HMM/interim_tables/03_results_noOutlier.csv',row.names = FALSE)
write.csv(df_mod,'../interim_tables/03_results_noOutlier.csv',row.names = FALSE)
