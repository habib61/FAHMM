#######################################################################################
###
### File name : 00_HMM_Wranling_Functions.R
### Program developer: aardepi1, Piet Aarden based on work of all below:
###                    S. Gardiner (Oxford BDI, stephen.gardiner@ndm.ox.ac.uk)
###                    H. Ganjgahi (Oxford BDI, habib.ganjgahi@bdi.ox.ac.uk)
###                    S. Pendleton (Oxford BDI, samantha.pendleton@ndm.ox.ac.uk)
### Date: 06-03-2023
### Project/Trial code: NO.MS - HMM
### Description: Function script which accompanies 00_HMM* scripts
### Application: R 4.1.0
### Libraries used: tidyverse_1.3.2, haven_2.5.1, reshape2_1.4.4
### Source table locations: /vob/CNOMSANON/anon/anon_2/anon_source/
###                         /funstorage/NSGDD_MS_MRI/HMM/00_source_data/
### Input: 
### Output: 
###
#######################################################################################


# Table join function -----------------------------------------------------

# Table join functions which aligns different component visits into one visit if within size window in days
table_join <- function(table1, table2, join_key=c("USUBJID", "STUDY", "DAY"), window_size=30) {
  
  table1 <- table1 %>% arrange(USUBJID, DAY) %>% mutate(TAB1DAY = DAY) %>% select(-any_of(c("DATE")))
  table2 <- table2 %>% arrange(USUBJID, DAY) %>% mutate(TAB2DAY = DAY) %>% select(-any_of(c("DATE")))
  
  dt1 <- data.table::data.table(table1, key = join_key)
  dt2 <- data.table::data.table(table2, key = join_key)
  
  leftjn <- dt1[dt2, roll='nearest'] %>%
    as.data.frame() %>%
    mutate(DIFF = abs(TAB1DAY - TAB2DAY)) %>%
    arrange(USUBJID, DIFF) %>%
    filter((DIFF <= window_size | is.na(DIFF))) %>%
    bind_rows(., table1) %>%
    arrange(USUBJID, DIFF) %>%
    filter((!duplicated(.[c("USUBJID", "TAB1DAY")]) | is.na(.[c("TAB1DAY")]))) %>% 
    bind_rows(., table2) %>%
    filter((!duplicated(.[c("USUBJID", "TAB2DAY")]) | is.na(.[c("TAB2DAY")]))) %>%
    select(-c(DIFF))
  
  rightjn <- dt2[dt1, roll='nearest'] %>%
    as.data.frame() %>%
    mutate(DIFF = abs(TAB2DAY - TAB1DAY)) %>%
    arrange(USUBJID, DIFF) %>%
    filter((DIFF <= window_size | is.na(DIFF))) %>%
    bind_rows(., table2) %>%
    arrange(USUBJID, DIFF) %>%
    filter((!duplicated(.[c("USUBJID", "TAB2DAY")]) | is.na(.[c("TAB2DAY")]))) %>% 
    bind_rows(., table1) %>%
    filter((!duplicated(.[c("USUBJID", "TAB1DAY")]) | is.na(.[c("TAB1DAY")]))) %>%
    select(-c(DIFF))
  
  mergdf <- bind_rows(leftjn, rightjn) %>%
    arrange(USUBJID, DAY, TAB1DAY, TAB2DAY) %>%
    filter((!duplicated(.[c("USUBJID", "DAY")]) | is.na(.["DAY"]))) %>%
    arrange(USUBJID, TAB1DAY, TAB2DAY) %>%
    filter((!duplicated(.[c("USUBJID", "TAB1DAY")]) | is.na(.["TAB1DAY"]) | is.na(.["DAY"]))) %>%
    arrange(USUBJID, TAB2DAY) %>%
    filter((!duplicated(.[c("USUBJID", "TAB2DAY")]) | is.na(.["TAB2DAY"]) | is.na(.["DAY"]))) %>%
    arrange(USUBJID, DAY) %>%
    filter(!duplicated(.)) %>% 
    select(-c(TAB1DAY, TAB2DAY))
  
  return(mergdf)
}


# Days since FDD function -------------------------------------------------

calculateDays <- function(data, FDD){
  
  #######################
  
  # This function calculates the days since FDD for each patient in the table
  
  # input:
  # data; dataframe containing the date column(s)
  # FDD; dataframe with two columns; USUBJID and FDD
  # output:
  # data.days; data with additions days column(s)
  
  #######################    
  
  # Get FDD
  FDD <- FDD %>% rename(STRTDATE = FDD)
  
  # check which date columns are present and carry out calculation for each patient
  patients <- unique(data$USUBJID)
  
  # for STRTDATE and ENDDATE present in data
  if(("STRTDATE" %in% colnames(data)) & ("ENDDATE" %in% colnames(data))){
    for(i in 1:length(patients)){
      ptsFDD <- FDD[which(FDD$USUBJID == patients[i]), "STRTDATE"]
      if(length(ptsFDD) == 0){next}
      data[which(data$USUBJID == patients[i]), "STRTDAY"] <- (data[which(data$USUBJID == patients[i]), "STRTDATE"] - ptsFDD) + 1
      data[which(data$USUBJID == patients[i]), "ENDDAY"] <- (data[which(data$USUBJID == patients[i]), "ENDDATE"] - ptsFDD) + 1
    }
    
    # Adjust for Novartis structure of no zeros in DAY
    data$STRTDAY <- ifelse(data$STRTDAY <= 0, data$STRTDAY - 1, data$STRTDAY)
    data$ENDDAY <- ifelse(data$ENDDAY <= 0, data$ENDDAY - 1, data$ENDDAY)
  }
  
  # for DATE present in data
  if("DATE" %in% colnames(data)){
    for(i in 1:length(patients)){
      ptsFDD <- FDD[which(FDD$USUBJID == patients[i]), "STRTDATE"]
      if(length(ptsFDD) == 0){next}
      data[which(data$USUBJID == patients[i]), "DAY"] <- (data[which(data$USUBJID == patients[i]), "DATE"] - ptsFDD) + 1
    }
    
    # Adjust for Novartis structure of no zeros in DAY
    data$DAY <- ifelse(data$DAY <= 0, data$DAY - 1, data$DAY)  
    
  }
  
  # for CONFDATE present in data (relapse table)
  if("CONFDATE" %in% colnames(data)){
    for(i in 1:length(patients)){
      ptsFDD <- FDD[which(FDD$USUBJID == patients[i]), "STRTDATE"]
      if(length(ptsFDD) == 0){next}
      data[which(data$USUBJID == patients[i]), "CONFDAY"] <- (data[which(data$USUBJID == patients[i]), "CONFDATE"] - ptsFDD) + 1
    }
    
    # Adjust for Novartis structure of no zeros in DAY
    data$CONFDAY <- ifelse(data$CONFDAY <= 0, data$CONFDAY - 1, data$CONFDAY)  
    
  }
  
  # for DATESCN1 and DATESCN2 present in data
  if(("DATESCN1" %in% colnames(data)) & ("DATESCN2" %in% colnames(data))){
    for(i in 1:length(patients)){
      
      ptsFDD <- FDD[which(FDD$USUBJID == patients[i]), "STRTDATE"][[1]]
      if(length(ptsFDD) == 0){next}
      data[which(data$USUBJID == patients[i]), "DAYSCN1"] <- (data[which(data$USUBJID == patients[i]), "DATESCN1"] - ptsFDD) + 1
      data[which(data$USUBJID == patients[i]), "DAYSCN2"] <- (data[which(data$USUBJID == patients[i]), "DATESCN2"] - ptsFDD) + 1
    }
    
    # Adjust for Novartis structure of no zeros in DAY
    data$DAYSCN1 <- ifelse(data$DAYSCN1 <= 0, data$DAYSCN1 - 1, data$DAYSCN1)
    data$DAYSCN2 <- ifelse(data$DAYSCN2 <= 0, data$DAYSCN2 - 1, data$DAYSCN2)
  }
  
  data.days <- data
  
  return(data.days)
}  


# AddRelapseInfo ----------------------------------------------------------
# Takes an endpoint df as input and adds a RELAPSE column on to it where it flags
# whether the endpoint was captured during a relapse

AddRelapseInfo <- function(df, rel){
  
  # Add relapse information onto table, if relapse was present during endpoint capture RELAPSE == 1, else RELAPSE == 0
  for(i in seq(1, nrow(df))){
    
    tmprel <- rel[which(rel$USUBJID == df[i, "USUBJID"]), c("STRTDAY", "ENDDAY")]
    arrflg <- c()
    if(nrow(tmprel) > 0){
      
      for(j in seq(1, nrow(tmprel))){
        iatj <- df[i, "DAY"] %in% seq(tmprel[j, "STRTDAY"], tmprel[j, "ENDDAY"])
        arrflg <- c(arrflg, iatj)
      }
    }
    
    df[i, "RELAPSE"] <- ifelse(any(arrflg), 1, 0)
  }
  
  return(df)
}

# AddTRTInfo ----------------------------------------------------------
# Takes an endpoint df as input and adds a ARM column on to it where it assigns
# the treatment the patient was assigned to while the endpoints was captured

AddTRTInfo <- function(df, follow){

  df$DATE <- as.Date(df$DATE)
  
  follow <- follow %>%
    mutate(ENDDATE = ifelse(is.na(ENDDATE), as.Date(Sys.Date()), ENDDATE),
           ENDDATE = as.Date(ENDDATE, "1970-01-01"))
  
  follow_core <- follow %>%
    filter(CORE == "Yes") %>%
    select(c(USUBJID, ARM)) %>%
    rename(COREARM = ARM)
  
  df <- df %>% 
    merge(., follow[c("USUBJID", "STRTDATE", "ENDDATE", "ARM")], by = "USUBJID", all.x = T) %>%
    mutate(KEEP = DATE %within% interval(ymd(STRTDATE), ymd(ENDDATE))) %>%
    arrange(USUBJID, DATE, desc(KEEP)) %>%
    mutate(ARM = ifelse(KEEP == FALSE, "No treatment", ARM)) %>%
    filter(!duplicated(.[c("USUBJID", "DATE")])) %>%
    merge(., follow_core, by = "USUBJID", all.x = T) %>%
    mutate(ARM = ifelse(DAY == 1, COREARM, ARM)) %>%
    select(-c(STRTDATE, ENDDATE, KEEP, COREARM))
  
  return(df)
}


# EndpointImputation ------------------------------------------------------

# Carries out imputation for different endpoints using GAM's.
# The formula always stays the same with only the endpoint and the family changing.
EndpointImputation <- function(data, fam, endpoint){
  
  library(mgcv)
  
  # Checking for missingness
  row.has.na <- apply(data[,c(endpoint, "YEARS", "MSTYPE", "RELAPSE", "ACTVTRT", "AGE",
                                "SEX", "DURFSLG", "RELPST1Y")], 1, function(x){any(is.na(x))})
  
  df_train  <- data[!row.has.na,]
  
  mdl <- bam(as.formula(paste0(endpoint,' ~ s(YEARS) + s(AGE) + SEX + s(DURFSLG) + RELPST1Y + MSTYPE + RELAPSE + ACTVTRT + s(USUBJID,YEARS, bs = "re") + s(USUBJID, bs = "re")')),
                 family = fam, data = df_train, discrete = TRUE)
  
  return(mdl)
  
}


# VW ----------------------------------------------------------------------
# This function provides the cutoffs for each nominal visit months
# i.e. provide an array with the months and it will provide you the cutoffs.
vw <- function(x=6){
  exact   <- x*365.25/12
  target  <- floor(exact)
  
  lower   <- c(1, (floor((target + c(target[-1], 0))/2)+1)[-length(target)])
  last.elem <- target[length(target)]+target[length(target)]-lower[length(lower)]
  upper   <- c(floor((target+c(0, target[-length(target)]))/2)[-1],  last.elem)
  dat      <- data.frame(Month=round(x,1), Year=round(x/12,1),  lower, target, upper)
  return( list(exact=exact, dat=dat))
}


# AssignMonth -------------------------------------------------------------
# Takes the cutoffs from the vw() functions and assigns the Month to the endpoints
# based on the DAY column
AssignMonth <- function(data, mco){
  
  data$MONTH <- NA
  for(i in 1:nrow(data)){
    
    day <- data[i, "DAY"]
    tmpmco <- mco
    tmpmco$day <- day
    tmpmco <- tmpmco %>%
      rowwise() %>%
      mutate(FLG = day %in% c(lower:upper)) %>%
      filter(FLG == TRUE)
    
    data[i, "MONTH"] <- tmpmco$Month
    
  }
  
  return(data)
}
