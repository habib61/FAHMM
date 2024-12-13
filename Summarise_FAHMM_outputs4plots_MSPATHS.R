###################################################
###################################################
###################################################
#####CV on discovery sample
Kk       = c(2,3,4,5,6,7,8,9,10,11,12,13)
LogLike22=matrix(0,5,length(Kk))
LogLiketrain=matrix(0,5,length(Kk))
for (v in 1:5) {
  for (i in 1:length(Kk)) {
    K        = Kk[i]
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_InitAllsample/NOMSv2_FAHMM_Relapse_CV_',v,'_loglike_',K,'.RData',sep=''))
    #tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
    LogLike22[v,i] = hmm_CTDC_test_LogLike$loglik
    
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_InitAllsample/NOMSv2_FAHMM_Relapse_CV_',v,'_train_',K,'.RData',sep=''))
    LogLiketrain[v,i]= hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
  }
}

AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLiketrain
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(5156)-2*LogLiketrain


plot(Kk,colMeans(AIC_1),xlab = 'Number of States',ylab = 'CV Discovery sample training sample AIC')
plot(Kk,colMeans(BIC_1),xlab = 'Number of States',ylab = 'CV Discovery sample training sample BIC')
plot(Kk,colMeans(LogLike22),xlab = 'Number of States',ylab = 'CV Discovery sample out of sample Log Likelihood')

###################################################
###################################################
####summarise the discovery and replication studies
#sumarise the discovery study
#read the discovery sample file
library(gtsummary)
library(tidyverse)
source('hmmsrc_CTDTApprox.R')

followDis         = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_NoImpT1gd3_Dis.csv')
followDis$Month   = floor(followDis$DAY/30)
followDis         = followDis %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()

Time              = followDis$deltaM
Time[is.na(Time)] = 0

ss                = followDis%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq               = ss$s

yy                = as.matrix(followDis[,c(39:42)])
yy                = scale(yy,center = FALSE,scale = TRUE)
P                 = ncol(yy)
Kk       = c(2,3,4,5,6,7,8,9,10,11,12,13)
LogLike2   = matrix(0,3,length(Kk))
for (i in 1:length(Kk)) {
  K        = Kk[i]
  #load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Dis_ClaraManh2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Dis_ClaraManh2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  LogLike2[3,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
  #load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Rep_Dis2_testLogLike/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Dis2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  LogLike2[1,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Rep_Dis2_ConvRel/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  LogLike2[2,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
}
LogLike = LogLike2
AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLike
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(nrow(followDis))-2*LogLike
#BIC_2 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(6442)-2*LogLike
#save them 4x6
plot(Kk,AIC_1,xlab = 'Number of States',ylab = 'Discovery sample AIC')
plot(Kk,BIC_1,xlab = 'Number of States',ylab = 'Discovery sample BIC')
#plot(Kk,BIC_2)

#summarise the discovery sample
followDis         = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_NoImpT1gd3_Dis.csv')
followDis$Month   = floor(followDis$DAY/30)
followDis         = followDis %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time              = followDis$deltaM
Time[is.na(Time)] = 0

ss                = followDis%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq               = ss$s

yy                = as.matrix(followDis[,c(39:42)])
yy                = scale(yy,center = FALSE,scale = TRUE)

K        = 8
#load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_Conj_noOutlier_ImpT1gd2_InitKM_diffPrior_Dis_InitDis2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
#load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Dis_ClaraManh2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Dis2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
l              = order(hmm_CTDC$mu[4,],decreasing = TRUE)
l=c(6,8,7,9,1,4,2,3,5)
#l=c(6,1,3,5,2,7,9,4,8)
#l              = c(3,8,1,6,7,5,9,2,4)
#l              = c(3,8,1,9,5,7,6,2,4)
#l              = c(5,4,1,8,9,6,3,10,2,7)
hmm_CTDC$mu    = hmm_CTDC$mu[,l]
hmm_CTDC$A     = hmm_CTDC$A[l,l]
hmm_CTDC$pi    = hmm_CTDC$pi[l]
hmm_CTDC$sigma = hmm_CTDC$sigma[,,l]


Z              = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)
spl            = split(followDis[,c(6:13)], Z)
#spl            = split(followDis[,c(6,30:33,37:38,13)], Z)
Emp_mu_S9      = sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]  = Emp_mu_S9[6,]/1000000
Emp_mu_S9D     = round(Emp_mu_S9,2)
Emp_mu_S9D
AD=hmm_CTDC$A
AD[AD<.01]     = 0
AD
followDis$Z=Z

#followDis%>%group_by(MSTYPE)%>%filter(Z==9& NUMGDT1>0)%>%summarise(mean(T25FWM,na.rm=TRUE),mean(HPT9M,na.rm=TRUE),mean(PASAT,na.rm=TRUE),mean(VOLT2,na.rm=TRUE),mean(NBV2,na.rm=TRUE),mean(NUMGDT1,na.rm=TRUE),mean(EDSS,na.rm=TRUE),length(MSTYPE))

#followDis%>%group_by(MSTYPE)%>%filter(Z==9& NUMGDT1==0)%>%summarise(mean(T25FWM,na.rm=TRUE),mean(HPT9M,na.rm=TRUE),mean(PASAT,na.rm=TRUE),mean(VOLT2,na.rm=TRUE),mean(NBV2,na.rm=TRUE),mean(NUMGDT1,na.rm=TRUE),mean(EDSS,na.rm=TRUE),length(MSTYPE))
mu            = hmm_CTDC$mu
row.names(mu) = c('relapseT1gd','relapse','brainCog','disability')
write.csv(mu,'NOMSv2/Final_Plots_Cro2/Discovery/Dis_mu_s9.csv')
write.csv(Emp_mu_S9,'NOMSv2/Final_Plots_Cro2/Discovery/Emp_Dis_mu_s9.csv')
write.csv(hmm_CTDC$A,'NOMSv2/Final_Plots_Cro2/Discovery/A_Dis_s9.csv')
write.csv(hmm_CTDC$mu,'NOMSv2/Final_Plots_Cro2/Discovery/States_mean_Dis_s9.csv')
write.csv(hmm_CTDC$pi,'NOMSv2/Final_Plots_Cro2/Discovery/Pi_Dis_s9.csv')


ff2<-followDis%>%dplyr::select(c("AGE","SEX","MSTYPE","DURFS","RELPST1Y","Z"))%>%tbl_summary(by=Z, missing = "no",
                                                                                         label = list(AGE~'Age',SEX~'Gender',MSTYPE='MS Type',DURFS='Duration from First Symptom',RELPST1Y="Number of Relapses during past year"),
                                                                                         type = list(all_continuous() ~ "continuous2", DURFS~"continuous2",RELPST1Y~"continuous2"),
                                                                                         statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})")) %>%modify_header(label ~ "**Variable**")
                                                                                                                                                                      
                                                                                                                                                                                                       
tt=ff2%>%gtsummary::as_tibble()
write.csv(tt,'NOMSv2/Final_Plots_Cro2/Discovery/Summary_demo_s9.csv')

ff3<-followDis%>%dplyr::select(c("MSTYPE","EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV2","NUMGDT1","RELAPSE","Z"))%>%tbl_summary(by=Z, missing = "no",
                                                                                                                        label = list(HPT9M~'Hand Coordination',T25FWM~'Walking Time',VOLT2~'T2 lesion volume',NBV2='Normalised brain volume',NUMGDT1='Number of T1 Gd lesions'),
                                                                                                                                type = list(all_continuous() ~ "continuous2"),
                                                                                                                                statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})","({p10}, {p90})")) %>%modify_header(label ~ "**Variable**")
                                                                                                                                                                        
tt=ff3%>%gtsummary::as_tibble()
write.csv(tt,'NOMSv2/Final_Plots_Cro2/Discovery/Summary_var_s9.csv')


#meta states
followDis$Ms                               = followDis$Z
followDis$Ms[followDis$Ms %in% c(1,2,3,4)]   = 1
followDis$Ms[followDis$Ms==5]              = 2
followDis$Ms[followDis$Ms==6]              = 3
followDis$Ms[followDis$Ms %in% c(9,7,8)]   = 4

A_Ms        = matrix(0,4,4)

AD          = hmm_CTDC$A
A_Ms[1,1]   = .962
A_Ms[1,2]   = .038
A_Ms[1,3]   = 0
A_Ms[1,4]   = 0


A_Ms[2,1]   = sum(AD[5,1:4])
A_Ms[2,2]   = AD[5,5]
A_Ms[2,3]   = AD[5,6]       
A_Ms[2,4]   = sum(AD[5,7:9])

A_Ms[3,1]   = sum(AD[6,1:4])
A_Ms[3,2]   = AD[6,5]
A_Ms[3,3]   = AD[6,6]       
A_Ms[3,4]   = sum(AD[6,7:9])

A_Ms[4,1]   = 0
A_Ms[4,2]   = .035
A_Ms[4,3]   = .024
A_Ms[4,4]   = .941

write.csv(A_Ms,'NOMSv2/Final_Plots_Cro2/Discovery/A_Meta_Dis_S9.csv')

TP1<- melt(A_Ms)
#thresholding
TP1 <- TP1[TP1$value!=0,]
TP1 <- TP1[TP1$value>=0.02,]
colnames(TP1) <- c("from", "to", "thickness")
TP1 <- subset(TP1,from!=to)
#Below is the final save pdf 6x6
#aa=matrix( c(-4.41,0.85,-1.59,-1.95,-1.00,4.0,6.34,7.00,9.62,0.98,0.65,2.00,-0.14,-1.50,-1.00,2.00,0.00,-2.00),9,2)

#qgraph(TP1, esize=5 ,edge.labels=T,label.cex=1,edge.label.cex=.5,edge.label.margin = 0.005,
#       edge.label.color='#3a338f',edge.color="#3a338f",layoutScale=c(1,1),GLratio=4, theme="Borkulo",layout = "spring",
#       groups=list(`Late MS`=c(6),`Early MS`=c(1),`Acute Relapse`=c(2),`Transition States`=c(3,4,5)),palette='ggplot2',legend.cex=.5,cut=.03,
#       directed=TRUE, legend=TRUE, legend.cex=.6)

#qgraph(TP1, esize=5 ,edge.labels=T,label.cex=1,edge.label.cex=.5,edge.label.margin = 0.005,
#       edge.label.color='#3a338f',edge.color="#3a338f",layoutScale=c(1,1),GLratio=4, theme="Borkulo",layout = "spring",
#       groups=list(`Late MS`=c(4),`Symptomatic Transition`=c(2),`Symptomatic Transition`=c(3),`Early MS`=c(1)),
#       palette='ggplot2',legend.cex=.5,cut=.03,
#       directed=TRUE, legend=TRUE)

aa=matrix(c(-4,-2,-2,-0.2,12.5,9,14,12.5),4,2)
qgraph(TP1, esize=5 ,edge.labels=T,label.cex=1,edge.label.cex=.4,edge.label.margin = 0.005,
       edge.label.color='#3a338f',edge.color="#3a338f",layoutScale=c(1,1),GLratio=6, theme="Borkulo",layout = aa,
       groups=list(`Early MS`=c(1),`Acute Relapse`=c(2),`Transition`=c(3),`Late MS`=c(4)),
       palette='ggplot2',legend.cex=.4,cut=.03,
       directed=TRUE, legend=TRUE)


aa2=matrix(c(-5,-5,-5,-2,-2,1,1,1,12,10,8,7,14,12,10,8),8,2)
TP1<- melt(AD)
#thresholding
TP1 <- TP1[TP1$value!=0,]
TP1 <- TP1[TP1$value>=0.05,]
colnames(TP1) <- c("from", "to", "thickness")
TP1 <- subset(TP1,from!=to)
qgraph(TP1, esize=5 ,edge.labels=T,label.cex=1,edge.label.cex=.5,edge.label.margin = 0.005,
       edge.label.color='#3a338f',edge.color="#3a338f",layoutScale=c(1,1),GLratio=4, theme="Borkulo",layout = aa2,
       groups=list(`Late MS`=c(6,7,8),`Asymptomatic Transition`=c(4),`Symptomatic Transition`=c(5),`Early MS`=c(1,2,3)),
       palette='ggplot2',legend.cex=.4,cut=.03,
       directed=TRUE, legend=TRUE)



qgraph(TP1, esize=5 ,edge.labels=T,label.cex=1,edge.label.cex=.5,edge.label.margin = 0.005,
       edge.label.color='#3a338f',edge.color="#3a338f",layoutScale=c(1,1),GLratio=4, theme="Borkulo",layout = "spring",
       groups=list(`Late MS`=c(7,8,9),`Acute Relapse`=c(4),`Transition`=c(6),`Early MS`=c(1,2,3,4)),
       palette='ggplot2',legend.cex=.4,cut=.03,
       directed=TRUE, legend=TRUE)
aa2=matrix(c(-5,-5,-5,-5,-2,-2,1,1,1,12.5,10.5,8,6,5,14,12,10,8),9,2)

aa2=matrix(c(-5,-5,-5,-5,-2,-2,1,1,1,12.5,10.5,8,6,5,14,11.5,9.5,7.5),9,2)
qgraph(TP1, esize=5 ,edge.labels=T,label.cex=1,edge.label.cex=.5,edge.label.margin = 0.005,
       edge.label.color='#3a338f',edge.color="#3a338f",layoutScale=c(1,1),GLratio=4, theme="Borkulo",layout = aa2,
       groups=list(`Late MS`=c(7,8,9),`Acute Relapse`=c(5),`Transition`=c(6),`Early MS`=c(1,2,3,4)),
       palette='ggplot2',legend.cex=.4,cut=.03,
       directed=TRUE, legend=TRUE)


spl <- split(followDis[,c(6:13)], followDis$Ms)
#followDis$VOLT2p=followDis$VOLT2p^3
#followDis$NUMGDT1p2=followDis$NUMGDT1p2^2
#spl <- split(followDis[,c(6,30:33,37:38,13)], Z)

Emp_Meta_S9=sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_Meta_S9[6,]=Emp_Meta_S9[6,]/1000000
round(Emp_Meta_S9,2)

write.csv(Emp_Meta_S9,'NOMSv2/Final_Plots_Cro2/Discovery/Dis_metaMean_s9.csv')


ff2<-followDis%>%dplyr::select(c("AGE","SEX","MSTYPE","DURFS","RELPST1Y","Ms"))%>%tbl_summary(by=Ms, missing = "no",
                                                                                             label = list(AGE~'Age',SEX~'Gender',MSTYPE='MS Type',DURFS='Duration from First Symptom',RELPST1Y="Number of Relapses during past year"),
                                                                                             type = list(all_continuous() ~ "continuous2", DURFS~"continuous2",RELPST1Y~"continuous2"),
                                                                                             statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})")) %>%modify_header(label ~ "**Variable**")


tt=ff2%>%gtsummary::as_tibble()
write.csv(tt,'NOMSv2/Final_Plots_Cro2/Discovery/Summary_Meta_demo_s9.csv')

ff3<-followDis%>%dplyr::select(c("MSTYPE","EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV2","NUMGDT1","RELAPSE","Ms"))%>%tbl_summary(by=Ms, missing = "no",
                                                                                                                                 label = list(HPT9M~'Hand Coordination',T25FWM~'Walking Time',VOLT2~'T2 lesion volume',NBV2='Normalised brain volume',NUMGDT1='Number of T1 Gd lesions'),
                                                                                                                                 type = list(all_continuous() ~ "continuous2"),
                                                                                                                                 statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})","({p10}, {p90})")) %>%modify_header(label ~ "**Variable**")

tt=ff3%>%gtsummary::as_tibble()
write.csv(tt,'NOMSv2/Final_Plots_Cro2/Discovery/Summary_Meta_var_s9.csv')


##############################################################
##############################################################
# summarise the replication study
followRep         = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_NoImpT1gd3_rep.csv')
followRep$Month   = floor(followRep$DAY/30)
followRep         = followRep %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time              = followRep$deltaM
Time[is.na(Time)] = 0

ss                = followRep%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq               = ss$s

yy                = as.matrix(followRep[,c(39:42)])
yy                = scale(yy,center = FALSE,scale = TRUE)

i        = 7
K        = Kk[i]
#replication optimisation was initialised from discovery sample
#load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Rep_ClaraManh2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))

load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3_Rep_Dis2/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
l  = order(hmm_CTDC$mu[4,],decreasing = TRUE)
l=c(8,6,7,9,1,4,3,2, 5)
#l  = c(6,8,3,7,5,4,2,1)
hmm_CTDC$mu=hmm_CTDC$mu[,l]
hmm_CTDC$A=hmm_CTDC$A[l,l]
hmm_CTDC$pi=hmm_CTDC$pi[l]
hmm_CTDC$sigma=hmm_CTDC$sigma[,,l]



Z        = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)

followRep$Z=Z
spl <- split(followRep[,c(6:13)], Z)
Emp_mu_S9_rep=sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9_rep[6,]=Emp_mu_S9_rep[6,]/1000000
Emp_mu_S9_rep[7,2]=0
round(Emp_mu_S9_rep,2)
AR=hmm_CTDC$A
AR[AR<.01]=0
AR

followRep$Z=Z

mu    = hmm_CTDC$mu
row.names(mu)=c('relapseT1gd','relapse','brainCog','disability')
write.csv(mu,'NOMSv2/Final_Plots_Cro2/Replication/Rep_mu_s9.csv')
write.csv(Emp_mu_S9_rep,'NOMSv2/Final_Plots_Cro2/Replication/Emp_Rep_mu_s9_initDis.csv')
write.csv(hmm_CTDC$A,'NOMSv2/Final_Plots_Cro2/Replication/A_Rep_s9_initDis.csv')
write.csv(hmm_CTDC$mu,'NOMSv2/Final_Plots_Cro2/Replication/States_mean_Rep_initDis_s9.csv')
write.csv(hmm_CTDC$pi,'NOMSv2/Final_Plots_Cro2/Replication/Pi_Rep_initDis_s9.csv')

ff2<-followRep%>%dplyr::select(c("AGE","SEX","MSTYPE","DURFS","RELPST1Y","Z"))%>%tbl_summary(by=Z, missing = "no",
                                                                                             label = list(AGE~'Age',SEX~'Gender',MSTYPE='MS Type',DURFS='Duration from First Symptom',RELPST1Y="Number of Relapses during past year"),
                                                                                             type = list(all_continuous() ~ "continuous2", DURFS~"continuous2",RELPST1Y~"continuous2"),
                                                                                             statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})")) %>%modify_header(label ~ "**Variable**")


tt=ff2%>%gtsummary::as_tibble()
write.csv(tt,'NOMSv2/Final_Plots_Cro2/Replication/Summary_demo_s9.csv')

ff3<-followRep%>%dplyr::select(c("MSTYPE","EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV2","NUMGDT1","RELAPSE","Z"))%>%tbl_summary(by=Z, missing = "no",
                                                                                                                                 label = list(HPT9M~'Hand Coordination',T25FWM~'Walking Time',VOLT2~'T2 lesion volume',NBV2='Normalised brain volume',NUMGDT1='Number of T1 Gd lesions'),
                                                                                                                                 type = list(all_continuous() ~ "continuous2"),
                                                                                                                                 statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})","({p10}, {p90})")) %>%modify_header(label ~ "**Variable**")

tt=ff3%>%gtsummary::as_tibble()
write.csv(tt,'NOMSv2/Final_Plots_Cro2/Replication/Summary_var_s9.csv')


########################################################################
####summarise the whole sample
library(gtsummary)
source('hmmsrc_CTDTApprox.R')

follow         = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_NoImpT1gd3.csv')
follow$Month   = floor(follow$DAY/30)
follow         = follow %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time           = follow$deltaM
Time[is.na(Time)] = 0

ss                = follow%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq               = ss$s

yy                = as.matrix(follow[,c(39:42)])
yy                = scale(yy,center = FALSE,scale = TRUE)

K        = 9
load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_NoImpT1gd3/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
l              = order(hmm_CTDC$mu[4,],decreasing = TRUE)
#l              = c(8,7,3,1,4,5,6,2)
hmm_CTDC$mu    = hmm_CTDC$mu[,l]
hmm_CTDC$A     = hmm_CTDC$A[l,l]
hmm_CTDC$pi    = hmm_CTDC$pi[l]
hmm_CTDC$sigma = hmm_CTDC$sigma[,,l]


Z              = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)
spl            = split(follow[,c(6:13)], Z)
Emp_mu_S9      = sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]  = Emp_mu_S9[6,]/1000000
Emp_mu_S9     = round(Emp_mu_S9,2)
Emp_mu_S9
A=hmm_CTDC$A
A[A<.01]     = 0
A
follow$Z=Z



########################################################################

#calculate the loglikelihood in the rep study
LogLike2   = matrix(0,1,length(Kk))
LogPost2   = matrix(0,1,length(Kk))
followRep  = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_ImpT1gd2_rep.csv')
yyR=followRep[,39:42]
yyR=scale(yyR,center = FALSE,scale=TRUE)
followRep$Month=round(followRep$DAY/30)
followRep = followRep %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time   = followRep$deltaM

Time[is.na(Time)]=0

ss     = followRep%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq    = ss$s

for (i in 1:length(Kk)) {
  K        = Kk[i]
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_Conj_noOutlier_ImpT1gd2_InitKM_diffPrior_Dis/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
  LogLike2[1,i] = tmp$loglik
  LogPost2[1,i] = tmp$logPost
}

AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLike2
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(nrow(yyR))-2*LogLike2
BIC_2 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(length(seq))-2*LogLike2
plot(Kk,AIC_1,xlab = 'Number of States',ylab = 'Out of sample AIC')
plot(Kk,BIC_1,xlab = 'Number of States',ylab = 'Out of sample BIC')
plot(Kk,BIC_2)

######################################################################################
######################################################################################
#allocate the timepoints into the meta-states, meta-states transition probability and summary
#read the fit from all sample
follow            = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_ImpT1gd2_CountNoTruncT1gd.csv') 
follow$Month      = floor(follow$DAY/30)
follow            = follow %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time              = follow$deltaM 
Time[is.na(Time)] = 0

ss                = follow%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq               = ss$s

yy                = as.matrix(follow[,39:42])
yy                = scale(yy,center = FALSE,scale = TRUE)


K                 = 9
load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_ImpT1gd2_InitClaraManh_CountNoTruncT1gd/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))

l  = order(hmm_CTDC$mu[4,],decreasing = TRUE)
l  = c(9,8,7,6,4,1,5,2,3)
hmm_CTDC$mu=hmm_CTDC$mu[,l]
hmm_CTDC$A=hmm_CTDC$A[l,l]
hmm_CTDC$pi=hmm_CTDC$pi[l]
hmm_CTDC$sigma=hmm_CTDC$sigma[,,l]


Z        = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)
spl <- split(follow[,c(6:13)], Z)
#followDis$VOLT2p=followDis$VOLT2p^3
#followDis$NUMGDT1p2=followDis$NUMGDT1p2^2
#spl <- split(followDis[,c(6,30:33,37:38,13)], Z)

Emp_mu_S9=sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]=Emp_mu_S9[6,]/1000000
round(Emp_mu_S9,2)
A=hmm_CTDC$A
A[A<.01]=0
A

follow$Z=Z

follow$Ms                             = Z
follow$Ms[follow$Ms %in% c(1,2,3)]   = 1
follow$Ms[follow$Ms==4]              = 2
follow$Ms[follow$Ms==5]              = 3
follow$Ms[follow$Ms==6]              = 4
follow$Ms[follow$Ms==7]              = 5
follow$Ms[follow$Ms %in% c(8,9)]     = 6

spl <- split(follow[,c(6:13)], follow$Ms)
#followDis$VOLT2p=followDis$VOLT2p^3
#followDis$NUMGDT1p2=followDis$NUMGDT1p2^2
#spl <- split(followDis[,c(6,30:33,37:38,13)], Z)

Emp_Meta_S9=sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_Meta_S9[6,]=Emp_Meta_S9[6,]/1000000
round(Emp_Meta_S9,2)


follow$Ms2                             = Z
follow$Ms2[follow$Ms2 %in% c(1,2,3)]  = 1
follow$Ms2[follow$Ms2==4]              = 2
follow$Ms2[follow$Ms2%in% c(5,6,7)]    = 3
follow$Ms2[follow$Ms2 %in% c(8,9)]     = 4
write.csv(follow,'derived_tables/follow_MetaStates.csv',row.names = FALSE)
A_Ms2=matrix(0,4,4)
A_Ms2[1,1]=sum(AD[1,1:3])
A_Ms2[1,2]=AD[1,4]
A_Ms2[1,3]=1-sum(A_Ms2[1,])
 
A_MS[4,4]=sum(AD[9,8:9])

write.csv(hazard.msm(ms.msm),file='hmm_trt_Ms.csv')
#covariate and treatment effect
load('NOMSv2/hmm_trt_All_Base_Ms.RData')

load('NOMSv2/hmm_trt_Ms.RData')
sojourn.msm(ms.msm)
totlos.msm(ms.msm,tot=10)
envisits.msm
efpt.msm
write.csv(hazard.msm(ms.msm),file='NOMSv2/Final_Plots_Cro/Discovery/hmm_trt_Ms.csv')
write.csv(sojourn.msm(ms.msm),file='NOMSv2/Final_Plots_Cro/Discovery/hmm_trt_sojourn_Ms.csv')

load('NOMSv2/hmm_trt_All_Base_Ms2.RData')

follow$year=follow$DAY/365
A = statetable.msm(follow$Ms,follow$USUBJID)
A = A/rowSums(A)
Q<-crudeinits.msm(Ms ~ year, subject=USUBJID, data = follow,qmatrix = A)
MS.msm <- msm( Ms ~ year, subject=USUBJID, data = follow,qmatrix =Q,control = list(fnscale = 63589,trace=1, REPORT=1,reltol = 1e-12,maxit = 10000),covariates = ~ ACTVTRT)



A = statetable.msm(follow$Ms2,follow$USUBJID)
A = A/rowSums(A)
Q<-crudeinits.msm(Ms2 ~ year, subject=USUBJID, data = follow,qmatrix = A)
MS.msm <- msm( Ms2 ~ year, subject=USUBJID, data = follow,qmatrix =Q,control = list(fnscale = 63589,trace=1, REPORT=1,reltol = 1e-12,maxit = 10000),covariates = ~ ACTVTRT)


followG%>%group_by(Ms,ACTVTRT)%>%summarise(length(ARM))

followG%>%group_by(Ms,ACTVTRT)%>%summarise(length(ACTVTRT))

A=statetable.msm(follow$Ms,followG$USUBJID)
A=A/rowSums(A)






follow            = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_ImpT1gd2_CountT1gd.csv') 
follow$Month      = floor(follow$DAY/30)
follow            = follow %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time              = follow$deltaM 
Time[is.na(Time)] = 0

ss                = follow%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq               = ss$s

yy                = as.matrix(follow[,39:42])
yy                = scale(yy,center = FALSE,scale = TRUE)

Kk       = c(2,3,4,5,6,7,8,9,10,11,12,13)
#Kk       = c(2,3,4,5,6,7,8,9,10,11,12)
#Kk       = c(2,3,4,5,6,7,8,9,10)
LogLike2   = matrix(0,1,length(Kk))
for (i in 1:length(Kk)) {
  K        = Kk[i]
  #load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_Conj_noOutlier_ImpT1gd2_InitKM_diffPrior_Rep/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  #LogLike2[1,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_Conj_noOutlier_ImpT1gd2_InitKM_diffPrior_Rep_InitDis22_CountT1gd/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  LogLike2[1,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
}


LogLike = LogLike2
AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLike
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(6442)-2*LogLike
#BIC_2 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(6442)-2*LogLike
plot(Kk,AIC_1,xlab = 'Number of States',ylab = 'Discovery sample AIC')
plot(Kk,BIC_1,xlab = 'Number of States',ylab = 'Discovery sample BIC')


i                 = 8
K                 = Kk[i]
load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_Conj_noOutlier_ImpT1gd2_InitKM_diffPrior_Rep_InitDis22_CountT1gd/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))

l  = order(hmm_CTDC$mu[4,],decreasing = TRUE)
#l  = c(4,5,8,3,2,1,9,7,6)
hmm_CTDC$mu=hmm_CTDC$mu[,l]
hmm_CTDC$A=hmm_CTDC$A[l,l]
hmm_CTDC$pi=hmm_CTDC$pi[l]
hmm_CTDC$sigma=hmm_CTDC$sigma[,,l]


Z        = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)
spl <- split(follow[,c(6:13)], Z)
#followDis$VOLT2p=followDis$VOLT2p^3
#followDis$NUMGDT1p2=followDis$NUMGDT1p2^2
#spl <- split(followDis[,c(6,30:33,37:38,13)], Z)

Emp_mu_S9=sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]=Emp_mu_S9[6,]/1000000
round(Emp_mu_S9,2)
A=hmm_CTDC$A
A[A<.01]=0
A

follow$Z=Z

ff4<-follow%>%dplyr::select(c('MSTYPE',"EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV2","NUMGDT1","RELAPSE","Z"))%>%tbl_summary(by=Z, missing = "no",
                                                                                                                        label = list(HPT9M~'Hand Coordination',T25FWM~'Walking Time',VOLT2~'T2 lesion volume',NBV2='Normalised brain volume',NUMGDT1='Number of T1 Gd lesions'),
                                                                                                                        type = list(all_continuous() ~ "continuous2"),
                                                                                                                        statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})")) %>%modify_header(label ~ "**Variable**")






follow            = read.csv('derived_tables/follow_FAhmm_Relapse_ALL_t1GdT2croNew_noOutlier_ImpT1gd2_CountT1gd.csv') 
follow$Month      = floor(follow$DAY/30)
follow            = follow %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time              = follow$deltaM 
Time[is.na(Time)] = 0

ss                = follow%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq               = ss$s

yy                = as.matrix(follow[,39:42])
yy                = scale(yy,center = FALSE,scale = TRUE)

Kk       = c(2,3,4,5,6,7,8,9,10,11,12,13)
#Kk       = c(2,3,4,5,6,7,8,9,10,11,12)
#Kk       = c(2,3,4,5,6,7,8,9,10)
LogLike2   = matrix(0,1,length(Kk))
for (i in 1:length(Kk)) {
  K        = Kk[i]
  #load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_Conj_noOutlier_ImpT1gd2_InitKM_diffPrior_Rep/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  #LogLike2[1,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_ImpT1gd2_InitKM_CountT1gd/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))
  LogLike2[1,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
}


LogLike = LogLike2
AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLike
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(6442)-2*LogLike
#BIC_2 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(6442)-2*LogLike
plot(Kk,AIC_1,xlab = 'Number of States',ylab = 'Discovery sample AIC')
plot(Kk,BIC_1,xlab = 'Number of States',ylab = 'Discovery sample BIC')


i                 = 8
K                 = Kk[i]
load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_t1gdT2cro_noOutlier_ImpT1gd2_InitKM_CountT1gd/NOMSv2_FAHMM_Relapse_Indp_',K,'.RData',sep=''))

l  = order(hmm_CTDC$mu[4,],decreasing = TRUE)
#l  = c(4,5,8,3,2,1,9,7,6)
hmm_CTDC$mu=hmm_CTDC$mu[,l]
hmm_CTDC$A=hmm_CTDC$A[l,l]
hmm_CTDC$pi=hmm_CTDC$pi[l]
hmm_CTDC$sigma=hmm_CTDC$sigma[,,l]


Z        = ctdthmm_MultSubj_viterbi(hmm_CTDC,yy,seq,Time)
table(Z)
spl <- split(follow[,c(6:13)], Z)
#followDis$VOLT2p=followDis$VOLT2p^3
#followDis$NUMGDT1p2=followDis$NUMGDT1p2^2
#spl <- split(followDis[,c(6,30:33,37:38,13)], Z)

Emp_mu_S9=sapply(spl, function(x){colMeans(x,na.rm = TRUE)})
Emp_mu_S9[6,]=Emp_mu_S9[6,]/1000000
round(Emp_mu_S9,2)
A=hmm_CTDC$A
A[A<.01]=0
A

follow$Z=Z

ff3<-follow%>%dplyr::select(c('MSTYPE',"EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV2","NUMGDT1","RELAPSE","Z"))%>%tbl_summary(by=Z, missing = "no",
                                                                                                                              label = list(HPT9M~'Hand Coordination',T25FWM~'Walking Time',VOLT2~'T2 lesion volume',NBV2='Normalised brain volume',NUMGDT1='Number of T1 Gd lesions'),
                                                                                                                              type = list(all_continuous() ~ "continuous2"),
                                                                                                                              statistic = all_continuous() ~ c("{mean} ({sd})","{median} ({p25}, {p75})")) %>%modify_header(label ~ "**Variable**")






LogLike2=matrix(0,2,length(Kk))
for (i in 1:length(Kk)) {
  K        = Kk[i]
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_rep2/NOMSv2_FAHMM_Relapse_Rep_',K,'.RData',sep=''))
  #tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
  LogLike2[1,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
  load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_rep/NOMSv2_FAHMM_Relapse_Dis_',K,'.RData',sep=''))
  LogLike2[2,i] = hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
}

AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLike2
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(nrow(yy))-2*LogLike2
BIC_2 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(length(seq))-2*LogLike2
plot(Kk,AIC_1,xlab = 'Number of States',ylab = 'Out of sample AIC')
plot(Kk,BIC_1,xlab = 'Number of States',ylab = 'Out of sample BIC')
plot(Kk,LogLike2)



########################################################################
#######################################################################
#######some work on CV
#calculate the loglikelihood in the rep study
LogLike2   = matrix(0,5,length(Kk))
LogPost2   = matrix(0,5,length(Kk))
followRep$Month=round(followRep$DAY/30)
followRep = followRep %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
Time   = followRep$deltaM

Time[is.na(Time)]=0

ss     = followRep%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
seq    = ss$s
ID     = list()
for (v in 1:5) {
  ID[[v]]=sample(unique(followDis$USUBJID),1280,replace = FALSE)
  tmp    = followDis%>%filter(USUBJID %in% ID[[v]])
  tmp    = tmp %>% group_by(USUBJID)%>%arrange(DAY, .by_group = TRUE)%>%mutate(deltaM=-lag(Month)+Month)%>%ungroup()
  Time   = tmp$deltaM
  Time[is.na(Time)]=0
  
  ss     = tmp%>%group_by(USUBJID)%>%summarise(s=length(USUBJID))
  seq    = ss$s
  yyR    = tmp[,21:24]
  yyR    = scale(yyR,center = FALSE,scale=TRUE)
  for (i in 1:length(Kk)) {
    K        = Kk[i]
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_dis/NOMSv2_FAHMM_Relapse_Dis_',K,'.RData',sep=''))
    tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
    LogLike2[v,i] = tmp$loglik
    LogPost2[v,i] = tmp$logPost
  }
}
plot(Kk,colMeans(LogLike2))

write.csv(ID[[1]],file = 'derived_tables/Discovery_CV_1.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[2]],file = 'derived_tables/Discovery_CV_2.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[3]],file = 'derived_tables/Discovery_CV_3.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[4]],file = 'derived_tables/Discovery_CV_4.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[5]],file = 'derived_tables/Discovery_CV_5.csv',row.names = FALSE,col.names = FALSE)



#CV using caret package
FF=createFolds(unique(followDis$USUBJID),k=5)
ll=unique(followDis$USUBJID)
ID     = list()
ID[[1]] = ll[FF$Fold1]
ID[[2]] = ll[FF$Fold2]
ID[[3]] = ll[FF$Fold3]
ID[[4]] = ll[FF$Fold4]
ID[[5]] = ll[FF$Fold5]

write.csv(ID[[1]],file = 'derived_tables/Discovery_CV_1.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[2]],file = 'derived_tables/Discovery_CV_2.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[3]],file = 'derived_tables/Discovery_CV_3.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[4]],file = 'derived_tables/Discovery_CV_4.csv',row.names = FALSE,col.names = FALSE)
write.csv(ID[[5]],file = 'derived_tables/Discovery_CV_5.csv',row.names = FALSE,col.names = FALSE)


AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLike2
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(nrow(yy))-2*LogLike2
BIC_2 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(length(seq))-2*LogLike2
plot(Kk,AIC_1,xlab = 'Number of States',ylab = 'Out of sample AIC')
plot(Kk,BIC_1,xlab = 'Number of States',ylab = 'Out of sample BIC')
plot(Kk,BIC_2)

##################################################################
##################################################################
######CV summary
LogLike2=matrix(0,5,length(Kk))
LogLiketest=matrix(0,5,length(Kk))
LogLiketrain=matrix(0,5,length(Kk))
for (v in 1:5) {
  for (i in 1:length(Kk)) {
    K        = Kk[i]
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_/NOMSv2_FAHMM_Relapse_CV_',v,'_loglike_',K,'.RData',sep=''))
    #tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
    LogLike2[v,i] = hmm_CTDC_test_LogLike$loglik
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV/NOMSv2_FAHMM_Relapse_CV_',v,'_train_',K,'.RData',sep=''))
    LogLiketrain[v,i]= hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
    
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV/NOMSv2_FAHMM_Relapse_CV_',v,'_test_',K,'.RData',sep=''))
    LogLiketest[v,i] = hmm_CTDC_test$loglik[length(hmm_CTDC_test$loglik)]
  }
}


Kk       = c(2,3,4,5,6,7,8,9,10,11,12,13)
LogLike2=matrix(0,5,length(Kk))
LogLike22=matrix(0,5,length(Kk))
LogLike=matrix(0,5,length(Kk))
#LogLiketest=matrix(0,5,length(Kk))
LogLiketrain=matrix(0,5,length(Kk))
for (v in 1:5) {
  for (i in 1:length(Kk)) {
    K        = Kk[i]
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_caret_randInit/NOMSv2_FAHMM_Relapse_CV_',v,'_loglike_',K,'.RData',sep=''))
    #tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
    LogLike2[v,i] = hmm_CTDC_test_LogLike$loglik
    
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_caret/NOMSv2_FAHMM_Relapse_CV_',v,'_loglike_',K,'.RData',sep=''))
    #tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
    LogLike22[v,i] = hmm_CTDC_test_LogLike$loglik
    
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_InitAllsample/NOMSv2_FAHMM_Relapse_CV_',v,'_loglike_',K,'.RData',sep=''))
    #tmp=hmm_DTMultSubj_MultV_BayesLoglike(yyR,hmm_CTDC,seq,K,Time,thresh)
    LogLike22[v,i] = hmm_CTDC_test_LogLike$loglik
    
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_caret_randInit/NOMSv2_FAHMM_Relapse_CV_',v,'_train_',K,'.RData',sep=''))
    LogLiketrain[v,i]= hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
    
    load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV_caret/NOMSv2_FAHMM_Relapse_CV_',v,'_train_',K,'.RData',sep=''))
    LogLike[v,i]= hmm_CTDC$loglik[length(hmm_CTDC$loglik)]
    
    
    #load(paste('NOMSv2/Dep_RELAPSE_MONTH_ABS_ALL_CV/NOMSv2_FAHMM_Relapse_CV_',v,'_test_',K,'.RData',sep=''))
    #LogLiketest[v,i] = hmm_CTDC_test$loglik[length(hmm_CTDC_test$loglik)]
  }
}




AIC_1 = (P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1-LogLiket
BIC_1 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(nrow(yy))-2*LogLike2
BIC_2 = ((P)^2*Kk+(P)*Kk+Kk*(Kk-1)+Kk-1)*log(length(seq))-2*LogLike2
plot(Kk,AIC_1,xlab = 'Number of States',ylab = 'Out of sample AIC')
plot(Kk,BIC_1,xlab = 'Number of States',ylab = 'Out of sample BIC')
plot(Kk,BIC_2)


