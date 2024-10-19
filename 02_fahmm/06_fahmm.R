###############################################
###############################################
###############################################
###############################################
###############################################
###############################################
###############################################
############### fit the HMM model############

library(tidyverse)
source('src/02_HMM_HMM_functions.R')

i=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
####load and prepare the data
follow       = read.csv('interim_tables/follow_FAhmm_Dis.csv')
#read the composite scores (disability,brain,relapse,Gd)
yy           = as.matrix(follow[,c("V34","V35","V36","V37")])
#yy           = scale(yy,center = FALSE,scale = TRUE)
yy=scale(yy)
follow  <- follow  %>%
  mutate(MONTHN = MONTH,
         MONTHN = ifelse(MONTHN == -1, 0, MONTHN)) %>%
  group_by(USUBJID) %>% 
  arrange(MONTHN, .by_group = TRUE) %>% 
  mutate(deltaM = -lag(MONTHN) + MONTHN) %>% 
  ungroup() %>%
  mutate(deltaM = ifelse(is.na(deltaM), 0, deltaM))

Time   = follow$deltaM
ss     = follow%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq=ss$s

##########################################
##########run multiple states for AIC, BIC
Kk       = c(2,3,4,5,6,7,8,9,10,11,12,13)
K        = Kk[i]

print(paste('fitting hmm model with states = ', K))

thresh   = .1
colnames(yy)

#init     = InitKM_mv(yy,K)
init     = InitKM_mv2(yy,K,follow$USUBJID,Time)
#init=Init_clara2(yy,K,follow$USUBJID,Time,'manhattan')
hmm_CTDC = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)

save(file=paste('interim_tables/Dis_InitKM_mv2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''),list=c('hmm_CTDC'))
