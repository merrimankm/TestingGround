##### Coxnet analysis #########
## Katie Merriman 3/7/23 ##

install.packages('survival')
install.packages('survminer')
install.packages('ROCR')
install.packages('pROC')
install.packages('glmnet')
install.packages('precrec')
install.packages('dplyr')
install.packages('Coxnet')



library(survival)
library(survminer)
library(ROCR)
library(pROC)
library(glmnet)
library(precrec)
library(ROCR)
library(precrec)
library(survival)
library(dplyr)
library(readr)
library(Coxnet)

#######################################################
############ LOAD AND PRE-PROCESS DATA ################
#######################################################
g<-"C:\\Users\\merrimankm\\Desktop\\SegLesionRadiomics_cleaned.csv"
data2<- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

g<-"T:\\MIP\\Katie_Merriman\\RadiomicsProject\\RadTest.csv"
data1<- read.csv(g, header=TRUE, stringsAsFactors=FALSE)


#checks for missing data
  any(is.na(data1))

#if any data is missing, find where
  which(is.na(data1), arr.ind=TRUE)
  which(is.character(cleanedData), arr.ind=TRUE)

  


#### SELECT ONLY NUMERIC DATA, with TIME, then STATUS, first #####
  cleanedData = data1[c(2:3, 8:7838,)]
  cleanedData = data1


#### PRE-PROCESS ORDINAL FEATURES ###### 

# convert from int to factor
  cleanedData$Race <- as.factor(cleanedData$Race)
  cleanedData$MRI_EPE <- as.factor(cleanedData$MRI_EPE)
  cleanedData$MRI_SVI <- as.factor(cleanedData$MRI_SVI)
  cleanedData$EPEscore <- as.factor(cleanedData$EPEscore)
  cleanedData$PIRADS.score.of.indexed.lesion <- as.factor(cleanedData$PIRADS.score.of.indexed.lesion)
  cleanedData$Total_GG <- as.factor(cleanedData$Total_GG)



#### PRE-PROCESS CONTINUOUS DATA #####  

# Change % positive cores from character to numeric #
  cleanedData$Total_perc <- parse_number(cleanedData$Total_perc)/100
  
# scale all continuous data
  cleanedData[c(1:8722)]<-as.numeric(cleanedData[c(1:8722)])
  cleanedData[c(1:8722)]<-scale(cleanedData[c(1:8722)])
  cleanedData[c(5:7,10:11,14:114)]<-scale(cleanedData[c(5:7,10:11,14:114)])
  
########################################################################
################  RUN COXNET ##########################################
##############################################################
  
  
  y = as.matrix(cleanedData[c(2:3)])
  
  
  noRad = cleanedData[c(8:21)]
  x_noRad = as.matrix(noRad)
  fit_noRad = Coxnet(x_noRad, y, Omega = NULL, penalty = "Lasso", nfolds = 20)
  
  Rad = cleanedData[c(8:7838,7840:8716)]
  x_Rad = as.matrix(Rad)
  fit_Rad = Coxnet(x_Rad, y, Omega = NULL, penalty = "Lasso", nfolds = 20)
  
  
  vars_noRad = which(fit_noRad$Beta0!=0)
  colnames(noRad[,c(vars_noRad)])
  coefs_noRad = fit_noRad$Beta[c(which(fit_noRad$Beta0!=0))]
  
  
  vars_Rad = which(fit_Rad$Beta0!=0)
  colnames(Rad[,c(vars_Rad)])
  coefs_Rad = fit_Rad$Beta[c(which(fit_Rad$Beta0!=0))]
  
  
  
  cleanedData$pred_noRads = coefs_noRad[1]*noRad[,vars_noRad[1]]
  print(cleanedData$pred_noRads[1:3])
  for (i in 2:length(vars_noRad)){
    print(cat("i = ", i))
    print(colnames(noRad[vars_noRad[i]]))
    print(noRad[1:3,vars_noRad[i]])
    print(coefs_noRad[i])
    print(coefs_noRad[i]*noRad[1:3,vars_noRad[i]])
    cleanedData$pred_noRads = cleanedData$pred_noRads + coefs_noRad[i]*noRad[,vars_noRad[i]]
    print(cleanedData$pred_noRads[1:3])
    }
  
g1 = roc(cleanedData$status,cleanedData$pred_noRads)
auc1 = paste("Standard AUC = ",toString(round(g1$auc, digits = 3)))
best.coords <- coords(g1, "best", best.method="youden")


cleanedData$pred_Rads = coefs_Rad[1]*Rad[,vars_Rad[1]]
print(cleanedData$pred_Rads[1:3])
for (i in 2:length(vars_Rad)){
  print(cat("i = ", i))
  print(colnames(Rad[vars_Rad[i]]))
  print(cleanedData[1:3,vars_Rad[i]])
  print(coefs_Rad[i])
  print(coefs_Rad[i]*Rad[1:3,vars_Rad[i]])
  cleanedData$pred_Rads = cleanedData$pred_Rads + coefs_Rad[i]*Rad[,vars_Rad[i]]
  print(cleanedData$pred_Rads[1:3])
}

g2 = roc(cleanedData$status,cleanedData$pred_Rads)   
auc2 = paste("Radiomics = ",toString(round(g2$auc, digits = 3)))
best.coords <- coords(g2, "best", best.method="youden")


plot(g1, col = 'purple')
text(0.9,0.15,'--',  col = 'purple')
text(0.48,.15,auc1,  col = 'purple')

plot(g2, add=TRUE, col='#2166AC')
text(0.9,0.10,'--', col = '#2166AC')
text(0.54,0.10,auc2,col='#2166AC')

roc.test(g1, g2)
  