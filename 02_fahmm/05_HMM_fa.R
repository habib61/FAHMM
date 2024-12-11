################################################################################################
################################################################################################
################################################################################################
################################################################################################
##################### predict the missing timepoints using marginal model#######################
########submit 04_prd script (sbatch PrdAll) to predict missing timepoints
################################################################################################
################################################################################################
################################################################################################
################################################################################################
#########now create a table with predicted values

#setwd("/data/users/qs9f68/HMM/final")
library(tidyverse)
#ff                              = read.csv('interim_tables/03_results_noOutlier.csv')
ff                              = read.csv('../interim_tables/03_results_noOutlier.csv')

List                            = c("T25FWM","HPT9M","PASAT","VOLT2","NBV2","t1gd")
List2                           = c("T25FWM","HPT9M","PASAT","VOLT2","NBV2","NUMGDT1")
Listp                           = c("T25FWMp","HPT9Mp","PASATp","VOLT2p","NBV2p","NUMGDT1p")
for (i in 1:5) {
  #load(paste('interim_tables/out_',List[i],'_re_pred_MissOnlyF.RData',sep=''))
  load(paste('../interim_tables/out_',List[i],'_re_pred_MissOnlyF.RData',sep=''))
    
  follow$USUBJID              = as.character(follow$USUBJID)
  follow                      = follow%>%dplyr::select(c("USUBJID","STUDY","DAY",List2[i],"pm"))
  follow[is.na(follow[,4]),4] = follow$pm[is.na(follow[,4])]
  follow                      = follow[,1:4]
  colnames(follow)            = c("USUBJID","STUDY","DAY",Listp[i])
  ff                          = left_join(ff,follow,by=c("USUBJID","STUDY","DAY"))
}

#write.csv(ff,file='interim_tables/longitudinal_noOutlier_NoMiss3.csv',row.names = FALSE)
write.csv(ff,file='../interim_tables/longitudinal_noOutlier_NoMiss3.csv',row.names = FALSE)

####################################################################################
####################################################################################
####################################################################################
#########################FA+HMM#####################################################
#2 step model: FA on baseline data to find the latent variable (MS composite scores) 
#and predict follw up latent variables and fit HMM
####################################################################################
#########################################################################
#########################################################################
####################################################################################

####baseline data FA analysis using the method from our beloved Veronica and copula

#source("/Users/h.ganjgahi/Desktop/HMM/NOMSv2/final/src/FACTOR_CODE_update.R")
source("../00a_source_scripts/01_HMM_FA_functions.R")


options(bitmapType='cairo-png')

#base         = read.csv('interim_tables/03_baseline_results.csv')
base         = read.csv('../interim_tables/03_baseline_results.csv')

base$NUMGDT1 = sqrt(base$NUMGDT1)
base$VOLT2   = base$VOLT2^(1/3)
bb           = base%>%filter(!is.na(PASAT))%>%dplyr::select(c( "EDSS","T25FWM","HPT9M","PASAT","VOLT2","NBV2","NUMGDT1","RELAPSE"))

Y            = scale(bb)
K            = ncol(Y)
p            = ncol(Y)
G            = ncol(Y)# Number of responses 	
startB       <- eigen(cov(Y))$vectors
alpha        <- 1/G
lambda1      <- 0.001
epsilon      <- 0.01

myImagePlot(abs(startB),F)
title("Initialization")

# PXEM: Dynamic Posterior Exploration (Approximate M-step)
start<-list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))

lambda0<-5
result_5<-FACTOR_ROTATE(Y,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-10
result_10<-FACTOR_ROTATE(Y,lambda0,lambda1,result_5,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-20
result_20<-FACTOR_ROTATE(Y,lambda0,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-30
result_30<-FACTOR_ROTATE(Y,lambda0,lambda1,result_20,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-40
result_40<-FACTOR_ROTATE(Y,lambda0,lambda1,result_30,K,epsilon,alpha,TRUE,TRUE,100,TRUE)
l=result_40$B
rownames(l)=colnames(Y)
l2=l
l2[result_40$P_star<.5]=0

#save(file = 'interim_tables/FA_baseline_Relapse.Rdata',list=c('result_5','result_10','result_20','result_30','result_40','start','Y'))
save(file = '../interim_tables/FA_baseline_Relapse.Rdata',list=c('result_5','result_10','result_20','result_30','result_40','start','Y'))

l=result_40$B
#order of latent varibales (disability,brain,relapse,Gd)
lambda=l[,c(7,4,3,2)]



################################################################################################
################################################################################################
########################project the data into MS Modes space and get the composite scores#######

#load('interim_tables/FA_baseline_Relapse.Rdata')
load('../interim_tables/FA_baseline_Relapse.Rdata')

#ff            = read.csv('interim_tables/longitudinal_noOutlier_NoMiss3.csv')
ff            = read.csv('../interim_tables/longitudinal_noOutlier_NoMiss3.csv')

#predict the missing number of Gd lesions from the FA model
ff$NUMGDT1p2  = sqrt((ff$NUMGDT1))
l             = result_40$B
#order of latent varibales (disability,brain,relapse,Gd)
lambda=l[,c(7,4,3,2)]

y             = as.matrix(ff%>%dplyr::select(c( "EDSS","T25FWMp","HPT9Mp","PASATp","VOLT2p","NBV2p","NUMGDT1p2","RELAPSE")))
y             = scale(y)
##find NA indices in each column
nT      = nrow(y)
p       = ncol(y)
IND     = matrix(FALSE,nT,p)
for (i in 1:p) {
  IND[,i]=is.na(y[,i])
}

#impue missing data and get composit scores
Sigma = result_40$B%*%t(result_40$B)+diag(result_40$sigma)

tmp_1 = solve(t(lambda)%*%diag(1/result_40$sigma)%*%lambda+diag(ncol(lambda)))
tmp_2 = t(lambda)%*%diag(1/result_40$sigma)
tmp   = tmp_1%*%tmp_2

fC    = matrix(0,nT,ncol(lambda))
yC    = y
for (i in 1:nT) {
  tmp2 = y[i,]
  if(any(IND[i,])){
    #tt            = yMean[IND[i,]]+Sigma[IND[i,],!IND[i,]]%*%solve(Sigma[!IND[i,],!IND[i,]])%*%(y[i,!IND[i,]]-yMean[!IND[i,]])
    tt            = Sigma[IND[i,],!IND[i,]]%*%solve(Sigma[!IND[i,],!IND[i,]])%*%(y[i,!IND[i,]])
    tmp2[IND[i,]] = tt
  }
  fC[i,]=tmp%*%tmp2
  yC[i,]=tmp2
  #fC[i,]=tmp2%*%tmp
}
ff[,34:37]=fC


#write.csv(ff,'interim_tables/follow_FAhmm_Relapse_ALL.csv',row.names = FALSE) 
write.csv(ff,'../interim_tables/follow_FAhmm_Relapse_ALL.csv',row.names = FALSE) 


################################################################################################
################################################################################################
################################################################################################
################################################################################################
#split the data into discovery and replication samples
library(tidyverse)
#follow              = read.csv('interim_tables/follow_FAhmm_Relapse_ALL.csv')
follow              = read.csv('../interim_tables/follow_FAhmm_Relapse_ALL.csv')

followB             = follow%>%group_by(USUBJID)%>%summarise(mV34=mean(V34),mV35=mean(V35),mV36=mean(V36),mV37=mean(V37))
#run kmeans clustering to find homogenous set of patients
#some functions
df2=scale(followB[,2:5])
wss <- function(k) {
  kmeans(df2, k, nstart = 30 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 20
k.values <- 1:20

# extract wss for 2-15 clusters
set.seed(1234)
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

ff          = kmeans(df2,3)
followB$val = ff$cluster
follow      = left_join(follow,followB,by='USUBJID')

followB%>%group_by(val)%>%summarise(length(USUBJID))

set.seed(1234)
RepID      = followB%>%group_by(val)%>%sample_frac(.2,replace = FALSE)%>%dplyr::select(USUBJID)
Dis        = follow%>%filter(! USUBJID %in% RepID$USUBJID)
Rep        = follow%>%filter(USUBJID %in% RepID$USUBJID)


#write.csv(Rep[,-c(27,38:43)],'interim_tables/follow_FAhmm_rep.csv',row.names = FALSE)
#write.csv(Dis[,-c(27,38:43)],'interim_tables/follow_FAhmm_Dis.csv',row.names = FALSE)
write.csv(Rep[,-c(27,38:43)],'../interim_tables/follow_FAhmm_rep.csv',row.names = FALSE)
write.csv(Dis[,-c(27,38:43)],'../interim_tables/follow_FAhmm_Dis.csv',row.names = FALSE)
