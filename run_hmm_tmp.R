source('FAHMM_Roche/00a_source_scripts/02_HMM_HMM_functions.R')

follow       = read.csv('interim_tables/follow_FAhmm_Relapse_ALLC.csv')

#read the composite scores (disability,brain,relapse,Gd)
yy           = as.matrix(follow[,c("V32","V33","V34","V35")])
yy           = scale(yy,center = FALSE,scale = TRUE)

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


thresh   = .01
colnames(yy)

init     = InitKM_mv(yy,K)
hmm_CTDC = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_km.RData',list = c('hmm_CTDC','init'))


init     = InitKM_mv2(yy,K,follow$USUBJID,Time)
hmm_CTDC_KM_mv2 = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_kmv2.RData',list = c('hmm_CTDC_KM_mv2','init'))

init     = Init_clara2(yy,K,follow$USUBJID,Time,'euclidean')
hmm_CTDC_clara2_euc = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara2_euc.RData',list = c('hmm_CTDC_clara2_euc','init'))



init     = Init_clara2(yy,K,follow$USUBJID,Time,'manhattan')
hmm_CTDC_clara2_manh = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara2_manh.RData',list = c('hmm_CTDC_clara2_manh','init'))


init     = Init_clara(yy,K,'euclidean')
hmm_CTDC_clara_euc = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara_euc.RData',list = c('hmm_CTDC_clara_euc','init'))



init     = Init_clara(yy,K,'manhattan')
hmm_CTDC_clara_manh = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara_manh.RData',list = c('hmm_CTDC_clara_manh','init'))



load("interim_tables/NOMSv2_FAHMM_Relapse_Indp_8.RData")
init = hmm_CTDC
ind=c(1,4,3,2)
init$mu = init$mu[ind,]
init$sigma = init$sigma[ind,ind,]

hmm_CTDC_noms = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_noms.RData',list = c('hmm_CTDC_noms','init'))








##########################################
##########################################
##########################################
##########################################


follow       = read.csv('interim_tables/follow_FAhmm_Relapse_ALL.csv')

#read the composite scores (disability,brain,relapse,Gd)
yy           = as.matrix(follow[,c("V32","V33","V34","V35")])
yy           = scale(yy,center = FALSE,scale = TRUE)

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

init     = InitKM_mv(yy,K)
hmm_CTDC_all = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_km_all.RData',list = c('hmm_CTDC_all','init'))


init     = InitKM_mv2(yy,K,follow$USUBJID,Time)
hmm_CTDCKM_mv2_all = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_kmv2_all.RData',list = c('hmm_CTDCKM_mv2_all','init'))


init     = Init_clara2(yy,K,follow$USUBJID,Time,'euclidean')
hmm_CTDC_clara2_euc_all = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara2_euc_all.RData',list = c('hmm_CTDC_clara2_euc_all','init'))



init     = Init_clara2(yy,K,follow$USUBJID,Time,'manhattan')
hmm_CTDC_clara2_manh_all = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara2_manh_all.RData',list = c('hmm_CTDC_clara2_manh_all','init'))


init     = Init_clara(yy,K,'euclidean')
hmm_CTDC_clara_euc_all = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara_euc_all.RData',list = c('hmm_CTDC_clara_euc_all','init'))



init     = Init_clara(yy,K,'manhattan')
hmm_CTDC_clara_manh_all = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_clara_manh_all.RData',list = c('hmm_CTDC_clara_manh_all','init'))



load("interim_tables/NOMSv2_FAHMM_Relapse_Indp_8.RData")
init = hmm_CTDC
ind=c(1,4,3,2)
init$mu = init$mu[ind,]
init$sigma = init$sigma[ind,ind,]

hmm_CTDC_noms_all = hmm_DTMultSubj_MultV_Bayes(yy,init,seq,K,Time,thresh)
save(file='interim_tables/hmm_s8_init_noms_all.RData',list = c('hmm_CTDC_noms_all','init'))
