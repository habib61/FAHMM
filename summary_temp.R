hmm=hmm_CTDC

hmm_CTDC=hmm

hmm_CTDC = hmm_CTDC_KM_mv2

hmm_CTDC = hmm_CTDC_all

hmm_CTDC = hmm_CTDC_clara2_euc 

hmm_CTDC = hmm_CTDC_noms

hmm_CTDC = hmm_CTDC_clara2_manh


follow       = read.csv('interim_tables/follow_FAhmm_Relapse_ALLC.csv')
follow=follow%>%filter(! USUBJID %in% 'WA21093-208650-1932631')
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



hmm_CTDC = hmm_CTDC_clara2_manh
l = order(hmm_CTDC$mu[1,],decreasing = TRUE)
hmm_CTDC$mu = hmm_CTDC$mu[,l]
hmm_CTDC$A  = hmm_CTDC$A[l,l] 
hmm_CTDC$pi = hmm_CTDC$pi[l]
hmm_CTDC$sigma = hmm_CTDC$sigma[,,l]


Z = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)

kk= follow%>%dplyr::select(c( "EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV","NUMGDT1","RELAPSE"))
spl=split(kk,Z)
Emp_mu_S9 = sapply(spl,function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]=Emp_mu_S9[6,]/1000

AD=hmm_CTDC$A
AD[AD<.005]=0
AD
Emp_mu_S9


hmm_CTDC = hmm_CTDC_clara2_euc
l = order(hmm_CTDC$mu[1,],decreasing = TRUE)
hmm_CTDC$mu = hmm_CTDC$mu[,l]
hmm_CTDC$A  = hmm_CTDC$A[l,l] 
hmm_CTDC$pi = hmm_CTDC$pi[l]
hmm_CTDC$sigma = hmm_CTDC$sigma[,,l]


Z = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)

kk= follow%>%dplyr::select(c( "EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV","NUMGDT1","RELAPSE"))
spl=split(kk,Z)
Emp_mu_S9 = sapply(spl,function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]=Emp_mu_S9[6,]/1000

AD=hmm_CTDC$A
AD[AD<.005]=0
AD
Emp_mu_S9






hmm_CTDC = hmm_CTDC_clara2_manh_all
l = order(hmm_CTDC$mu[1,],decreasing = TRUE)
hmm_CTDC$mu = hmm_CTDC$mu[,l]
hmm_CTDC$A  = hmm_CTDC$A[l,l] 
hmm_CTDC$pi = hmm_CTDC$pi[l]
hmm_CTDC$sigma = hmm_CTDC$sigma[,,l]


Z = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)

kk= follow%>%dplyr::select(c( "EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV","NUMGDT1","RELAPSE"))
spl=split(kk,Z)
Emp_mu_S9 = sapply(spl,function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]=Emp_mu_S9[6,]/1000

AD=hmm_CTDC$A
AD[AD<.005]=0
AD
Emp_mu_S9




hmm_CTDC = hmm_CTDC_clara_euc_all
l = order(hmm_CTDC$mu[1,],decreasing = TRUE)
hmm_CTDC$mu = hmm_CTDC$mu[,l]
hmm_CTDC$A  = hmm_CTDC$A[l,l] 
hmm_CTDC$pi = hmm_CTDC$pi[l]
hmm_CTDC$sigma = hmm_CTDC$sigma[,,l]


Z = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)

kk= follow%>%dplyr::select(c( "EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV","NUMGDT1","RELAPSE"))
spl=split(kk,Z)
Emp_mu_S9 = sapply(spl,function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]=Emp_mu_S9[6,]/1000

AD=hmm_CTDC$A
AD[AD<.005]=0
AD
Emp_mu_S9
