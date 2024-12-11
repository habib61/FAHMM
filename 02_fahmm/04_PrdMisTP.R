#setwd('/data/users/qs9f68/HMM/Piet/HMM')
Sys.setenv('R_MAX_VSIZE'=128000000000)
Sys.getenv('R_MAX_VSIZE')

set.seed(1234)

library(mgcv)
library(tidyverse)

FAMILY                              = c('gaussian','gaussian','poisson','nb','gaussian','gaussian')
i                                   = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
I                                   = i+18 

List                                = c("T25FWM","HPT9M","PASAT","NUMGDT1","VOLT2","NBV2")

#follow                              = read.csv('interim_tables/03_results_noOutlier.csv')
follow                              = read.csv('../interim_tables/03_results_noOutlier.csv')
follow$USUBJID                      = as.factor(follow$USUBJID)
#follow$VOLT2                       = (follow$VOLT2)^(1/3)

#split the complete dataset into train and test
#find the complete dataset with no missing, split into train and test
row.has.na                          = apply(follow[,c("AGE","SEX","MSTYPE","DURFSLG","RELPST1Y","RELAPSE","ACTVTRT",List[i])], 1, function(x){any(is.na(x))})
Train                               = follow[!row.has.na,]
Train$id                            = 1:nrow(Train)
kk                                  = Train%>%group_by(USUBJID)%>%summarise(ses=length(USUBJID))
Train_1                             = Train%>%filter(USUBJID %in% kk$USUBJID[kk$ses<3])
tmp                                 = Train%>%filter(!USUBJID %in% Train_1$USUBJID)
Train_2                             = tmp%>%group_by(USUBJID)%>%sample_frac(size=.7)

train                               = rbind(Train_1,Train_2)
test                                = tmp%>%filter(!id%in%Train_2$id)
#row.has.na                          = apply(follow[,c(6,8:9,12,13,16,17)], 1, function(x){any(is.na(x))})
#free some memory
rm('Train')
#now replace the missing T25 with the predictions
colnames(train[,c(6,I)])
dim(train)
print(I)
ff                                  = bam(train[,I]~s(YEARS)+s(AGE)+SEX+s(DURFSLG)+RELPST1Y+MSTYPE+RELAPSE+ACTVTRT+ s(USUBJID,YEARS, bs = "re")+ s(USUBJID, bs = "re"),
                                          family = FAMILY[i], data = train,discrete = TRUE)
summary(ff)
#save(file=paste('interim_tables/','out_',List[i],'_re.RData',sep=''),list=c('ff'))
save(file=paste('../interim_tables/','out_',List[i],'_re.RData',sep=''),list=c('ff'))

#evaluate train and test erros
dd                                  = predict.bam(ff,type='response',discrete = TRUE)
mean(abs(train[,I]-dd))
mean((train[,I]-dd)^2)
dd2                                 = predict.bam(ff,newdata=test,type='response',discrete = TRUE)
mean(abs(test[,I]-dd2))
mean((test[,I]-dd2)^2)
###########replace
#find missing covariates
row.has.na                          = apply(follow[,c("AGE","SEX","MSTYPE","DURFSLG","RELPST1Y","RELAPSE","ACTVTRT")], 1, function(x){any(is.na(x))})
follow$pm                           = NA
follow$pm[!row.has.na]              = predict.gam(ff,follow[!row.has.na,],allow.new.levels=TRUE,type='response',discrete = TRUE)

#save(file=paste('interim_tables/','out_',List[i],'_re_pred.RData',sep=''),list=c('test','train','dd','dd2','follow'))
save(file=paste('../interim_tables/','out_',List[i],'_re_pred.RData',sep=''),list=c('test','train','dd','dd2','follow'))

#save(file=paste('interim_tables/','out_',List[i],'_re_pred_MissOnlyF.RData',sep=''),list=c('follow'))
save(file=paste('../interim_tables/','out_',List[i],'_re_pred_MissOnlyF.RData',sep=''),list=c('follow'))






