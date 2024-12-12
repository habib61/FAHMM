library(mgcv)
library(tidyverse)

FAMILY                              = c('gaussian','gaussian','poisson','nb','gaussian','gaussian')
List                                = c("T25FWM","HPT9M","PASAT","NUMGDT1","VOLT2","NBV")
follow                              = read.csv('interim_tables/03_results_noOutlier.csv')
follow$USUBJID                      = as.factor(follow$USUBJID)


for (i in c(1,2,3,5,6)) {
  I                                   = i+16 
  row.has.na                          = apply(follow[,c("AGE","SEX","MSTYPE","DURFSLG","RELPST1Y","RELAPSE",List[i])], 1, function(x){any(is.na(x))})
  train                               = follow[!row.has.na,]
  print(colnames(train[,c(6,I)]))
  ff                                  = bam(train[,I]~s(YEARS)+s(AGE)+SEX+s(DURFSLG)+RELPST1Y+MSTYPE+RELAPSE+ s(USUBJID,YEARS, bs = "re")+ s(USUBJID, bs = "re"),
                                            family = FAMILY[i], data = train,discrete = TRUE)
  print(summary(ff))
  save(file=paste('interim_tables/','out_',List[i],'_re.RData',sep=''),list=c('ff'))
  #find missing covariates
  row.has.na                          = apply(follow[,c("AGE","SEX","MSTYPE","DURFSLG","RELPST1Y","RELAPSE")], 1, function(x){any(is.na(x))})
  follow$pm                           = NA
  follow$pm[!row.has.na]              = predict.gam(ff,follow[!row.has.na,],allow.new.levels=TRUE,type='response',discrete = TRUE)
  save(file=paste('interim_tables/','out_',List[i],'_re_pred_MissOnlyF.RData',sep=''),list=c('follow'))
}