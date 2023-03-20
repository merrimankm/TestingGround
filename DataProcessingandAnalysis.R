install.packages('survival')
install.packages('survminer')
install.packages('ROCR')
install.packages('pROC')
install.packages('glmnet')
install.packages('precrec')


library(survival)
library(survminer)
library(ROCR)
library(pROC)
library(glmnet)
library(precrec)
library(ROCR)
library(precrec)
library(survival)

#######################################################
############ LOAD AND PRE-PROCESS DATA ################
#######################################################
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\UPDATEDdata\\ManuscriptAnalysis\\FullStagingData.csv"
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\Radiomics_BCR\\FULLPatientLevelFeatures.csv"
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\UPDATEDdata\\ManuscriptAnalysis\\BiopsyData.csv"
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\UPDATEDdata\\ManuscriptAnalysis\\FullStagingData_noMargins.csv"
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\UPDATEDdata\\ManuscriptAnalysis\\RadiomicsData_noMargins.csv"
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\UPDATEDdata\\ManuscriptAnalysis\\BiopsyData_noMargins.csv"
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\UPDATEDdata\\Radiomics.csv"
  g<-"T:\\MIP\\Katie_Merriman\\Project1Data\\Radiomics_BCR\\CLEANEDPatientLevelFeatures.csv"
  
  
  data1<- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

  #checks for missing data
    any(is.na(data1$P3vol2wp))
    
    
#### PRE-PROCESS STATUS ######    
    
     data1$status_ordinal <- as.factor(data1$status)
     
#### PRE-PROCESS CLINICAL DATA #####    
     data1$PSA_cat <- cut(data1$PSA, 
                         breaks=c(-Inf,6,10,20,Inf), right = TRUE, 
                         labels = c('0-6', '6.01-10', '10.01-20', '>20'))
     data1$PSA_scaled <- scale(data1$PSA)    
     data1$PSA_density_scaled <- scale(data1$PSA_density) 
    
     data1$PIRADS.score.of.indexed.lesion <- as.numeric(data1$PIRADS.score.of.indexed.lesion)
     data1$PIRADS_cat <- cut(data1$PIRADS.score.of.indexed.lesion,
                            breaks=c(-Inf,3,4,Inf), right = TRUE, 
                            labels = c('</= 3', '4', '5')) 
     data1$PIRADS <- as.factor(data1$PIRADS.score.of.indexed.lesion)
     data1$PIRADS_4cat <- cut(data1$PIRADS.score.of.indexed.lesion,
                              breaks=c(-Inf,2,3,4,Inf), right = TRUE, 
                              labels = c('1/2','3', '4', '5')) 
    
    
    
#### PRE-PROCESS STAGING DATA #####
    
      data1$EPEscore <- as.factor(data1$EPEscore)
      data1$MRI_EPE <- as.factor(data1$MRI_EPE)
      data1$MRI_SVI <- as.factor(data1$MRI_SVI)
     
      data1$path_EPE <- as.factor(data1$path_EPE)
      data1$path_SVI <- as.factor(data1$path_SVI)
      data1$path_LNI <- as.factor(data1$path_LNI)
      
      data1$Tstage_agreement <- factor(data1$Tstage_agreement, labels = c('Pathology T2/ MRI T2', 'Pathology T2/ MRI T3', 'Pathology T3/ MRI T2', 'Pathology T3/ MRI T3'))
      data1$EPE_agreement <- factor(data1$EPE_agreement, labels = c('no path. EPE/ no MRI EPE', 'no path. EPE/ MRI EPE', 'path. EPE/ no MRI EPE', 'path. EPE/ MRI EPE'))
      data1$SVI_agreement <- factor(data1$SVI_agreement, labels = c('no path. SVI/ no MRI SVI', 'no path. SVI/ MRI SVI', 'path. SVI/ no MRI SVI', 'path. SVI/ MRI SVI'))
      
      
      
#### PRE-PROCESS GLEASON DATA ####
      
      data1$Gleason <- cut(data1$Total_GG, 
                           breaks=c(-Inf,1,2,3,Inf), right = TRUE, 
                           labels = c('benign/1', '2', '3', '4/5')) 
      
      data1$SbxISUP <- cut(data1$GleasonSx, 
                           breaks=c(-Inf,1,2,3,Inf), right = TRUE, 
                           labels = c('benign/1', '2', '3', '4/5')) 
      data1$FbxISUP <- cut(data1$GleasonFx,
                           breaks=c(-Inf,1,2,3,Inf), right = TRUE, 
                           labels = c('benign/1', '2', '3', '4/5')) 
      data1$CbxISUP <- cut(data1$biopsyGleason,
                           breaks=c(-Inf,1,2,3,Inf), right = TRUE, 
                           labels = c('benign/1', '2', '3', '4/5')) 
      data1$RP_ISUP <- cut(data1$GleasonRp,
                           breaks=c(-Inf,1,2,3,Inf), right = TRUE, 
                           labels = c('benign/1', '2', '3', '4/5')) 
      
      
      
#### PRE-PROCESS BCR PREDICTION MODEL DATA ####
      
      data1$UCSF_CAPRA_combined <- as.factor(data1$UCSF_CAPRA_combined)
      data1$UCSF_CAPRA_Sbx <- as.factor(data1$UCSF_CAPRA_Sbx)
      data1$CAPRA_S <- as.factor(data1$CAPRA_S)
      data1$modified_CAPRA_S_combined <- as.factor(data1$modified_CAPRA_S_combined)
      data1$modified_CAPRA_S_Sbx <- as.factor(data1$modified_CAPRA_S_Sbx)
      data1$UCSF_CAPRA_combined_risk <- as.factor(data1$UCSF_CAPRA_combined_risk)
      data1$UCSF_CAPRA_Sbx_risk <- as.factor(data1$UCSF_CAPRA_Sbx_risk)
      data1$CAPRA_S_risk <- as.factor(data1$CAPRA_S_risk)
      data1$modified_CAPRA_S_combined_risk <- as.factor(data1$modified_CAPRA_S_combined_risk)
      data1$modified_CAPRA_S_Sbx_risk <- as.factor(data1$modified_CAPRA_S_Sbx_risk)
      
      
      
#### PRE-PROCESS RADIOMICS DATA #####      
      
      data1$NumLesions_scaled <- scale(data1$NumLesions)
      data1$NumLesions_ordinal <- as.factor(data1$NumLesions)
      
      data1$prostVol_scaled <- scale(data1$prostVol)
      data1$totalVol_all_scaled <- scale(data1$totalVol_all)
      data1$totalVol_P5_scaled <- scale(data1$totalVol_P5)
      data1$P5toTotalVol_scaled <- scale(data1$P5toTotalVol)
      data1$P5toProstVol_scaled <- scale(data1$P5toProstVol)
      data1$totalVol_P4_scaled <- scale(data1$totalVol_P4)
      data1$P4toTotalVol_scaled <- scale(data1$P4toTotalVol)
      data1$P4toProstVol_scaled <- scale(data1$P4toProstVol)
      data1$totalVol_P3_scaled <- scale(data1$totalVol_P3)
      data1$P3toTotalVol_scaled <- scale(data1$P3toTotalVol)
      data1$P3toProstVol_scaled <- scale(data1$P3toProstVol)
      data1$indexSA_scaled <- scale(data1$indexSA)
      data1$indexSphericity_scaled <- scale(data1$indexSphericity)
      data1$indexSVR_scaled <- scale(data1$indexSVR)
      data1$indexmax3Ddiameter_scaled <- scale(data1$indexmax3Ddiameter)
      data1$indexElongation_scaled <- scale(data1$indexElongation)
      data1$indexFlatness_scaled <- scale(data1$indexFlatness)
      data1$prostT2_ent_scaled <- scale(data1$prostT2_ent)
      data1$indexT2_ent_scaled <- scale(data1$indexT2_ent)
      data1$maxT2_ent_scaled <- scale(data1$maxT2_ent)
      data1$minT2_ent_scaled <- scale(data1$minT2_ent)
      data1$T2skewness_scaled <- scale(data1$T2skewness)
      data1$prostADC_ent_scaled <- scale(data1$prostADC_ent)
      data1$indexADC_ent_scaled <- scale(data1$indexADC_ent)
      data1$maxADC_ent_scaled <- scale(data1$maxADC_ent)
      data1$minADC_ent_scaled <- scale(data1$minADC_ent)
      data1$prostHighB_ent_scaled <- scale(data1$prostHighB_ent)
      data1$highBskew_scaled <- scale(data1$highBskew)
      data1$indexHighB_ent_scaled <- scale(data1$indexHighB_ent)
      data1$maxHighB_ent_scaled <- scale(data1$maxHighB_ent)
      data1$minHighB_ent_scaled <- scale(data1$minHighB_ent)
      data1$prostT2_uniform_scaled <- scale(data1$prostT2_uniform)
      data1$indexT2_uniform_scaled <- scale(data1$indexT2_uniform)
      data1$maxT2_uniform_scaled <- scale(data1$maxT2_uniform)
      data1$minT2_uniform_scaled <- scale(data1$minT2_uniform)
      data1$prostADC_uniform_scaled <- scale(data1$prostADC_uniform)
      data1$indexADC_uniform_scaled <- scale(data1$indexADC_uniform)
      data1$maxADC_uniform_scaled <- scale(data1$maxADC_uniform)
      data1$minADC_uniform_scaled <- scale(data1$minADC_uniform)
      data1$prostHighB_uniform_scaled <- scale(data1$prostHighB_uniform)
      data1$indexHighB_uniform_scaled <- scale(data1$indexHighB_uniform)
      data1$maxHighB_uniform_scaled <- scale(data1$maxHighB_uniform)
      data1$minHighB_uniform_scaled <- scale(data1$minHighB_uniform)
      data1$prostADC_mean_scaled <- scale(data1$prostADC_mean)
      data1$prostADC_median_scaled <- scale(data1$prostADC_median)
      data1$prostADC_10th_scaled <- scale(data1$prostADC_10th)
      data1$prostADC_90th_scaled <- scale(data1$prostADC_90th)
      data1$prostADC_IQR_scaled <- scale(data1$prostADC_IQR)
      data1$minADC_mean_scaled <- scale(data1$minADC_mean)
      data1$variance_minADCmean_scaled <- scale(data1$variance_minADCmean)
      data1$minADC_median_scaled <- scale(data1$minADC_median)
      data1$variance_minADCmedian_scaled <- scale(data1$variance_minADCmedian)
      data1$minADC_10th_scaled <- scale(data1$minADC_10th)
      data1$minADC_IQR_scaled <- scale(data1$minADC_IQR)
      data1$indexADC_mean_scaled <- scale(data1$indexADC_mean)
      data1$indexADC_median_scaled <- scale(data1$indexADC_median)
      data1$indexADC_10th_scaled <- scale(data1$indexADC_10th)
      data1$indexADC_IQR_scaled <- scale(data1$indexADC_IQR)
      data1$prostHighB_mean_scaled <- scale(data1$prostHighB_mean)
      data1$prostHighB_median_scaled <- scale(data1$prostHighB_median)
      data1$prostHighB_10th_scaled <- scale(data1$prostHighB_10th)
      data1$prostHighB_90th_scaled <- scale(data1$prostHighB_90th)
      data1$prostHighB_IQR_scaled <- scale(data1$prostHighB_IQR)
      data1$maxHighB_mean_scaled <- scale(data1$maxHighB_mean)
      data1$variance_maxHighBmean_scaled <- scale(data1$variance_maxHighBmean)
      data1$maxHighB_median_scaled <- scale(data1$maxHighB_median)
      data1$variance_maxHighBmedian_scaled <- scale(data1$variance_maxHighBmedian)
      data1$maxHighB_90th_scaled <- scale(data1$maxHighB_90th)
      data1$maxHighB_IQR_scaled <- scale(data1$maxHighB_IQR)
      data1$indexHighB_mean_scaled <- scale(data1$indexHighB_mean)
      data1$indexHighB_median_scaled <- scale(data1$indexHighB_median)
      data1$indexHighB_90th_scaled <- scale(data1$indexHighB_90th)
      data1$indexHighB_IQR_scaled <- scale(data1$indexHighB_IQR)
      data1$t2IDM_scaled <- scale(data1$t2IDM)
      data1$adcIDM_scaled <- scale(data1$adcIDM)
      data1$highbIDM_scaled <- scale(data1$highbIDM)
      data1$t2IDMN_scaled <- scale(data1$t2IDMN)
      data1$adcIDMN_scaled <- scale(data1$adcIDMN)
      data1$highbIDMN_scaled <- scale(data1$highbIDMN)
      data1$t2ID_scaled <- scale(data1$t2ID)
      data1$adcID_scaled <- scale(data1$adcID)
      data1$highbID_scaled <- scale(data1$highbID)
      data1$t2IDN_scaled <- scale(data1$t2IDN)
      data1$adcIDN_scaled <- scale(data1$adcIDN)
      data1$highbIDN_scaled <- scale(data1$highbIDN)
      data1$t2GLN_scaled <- scale(data1$t2GLN)
      data1$adcGLN_scaled <- scale(data1$adcGLN)
      data1$highbGLN_scaled <- scale(data1$highbGLN)
      data1$t2GLNN_scaled <- scale(data1$t2GLNN)
      data1$adcGLNN_scaled <- scale(data1$adcGLNN)
      data1$highbGLNN_scaled <- scale(data1$highbGLNN)
      data1$t2ZoneEnt_scaled <- scale(data1$t2ZoneEnt)
      data1$adcZoneEnt_scaled <- scale(data1$adcZoneEnt)
      data1$highbZoneEnt_scaled <- scale(data1$highbZoneEnt)
      data1$t2ZonePerc_scaled <- scale(data1$t2ZonePerc)
      data1$adcZonePerc_scaled <- scale(data1$adcZonePerc)
      data1$highbZonePerc_scaled <- scale(data1$highbZonePerc)
      data1$t2Coarseness_scaled <- scale(data1$t2Coarseness)
      data1$adcCoarseness_scaled <- scale(data1$adcCoarseness)
      data1$highbCoarseness_scaled <- scale(data1$highbCoarseness)
      
      
#############################################################################
################  Box plots and violin plots #####################
######################################################
      library(ggplot2)
      library(viridis)
      library(RColorBrewer)
      library(wesanderson)
      library(ggsci)
      # library(dyplr)
      
      x1 = data1$PIRADS
      x2 = data1$status_ordinal
      
      
      
      y = data1$PSA
      name = "PSA [Continuous]"
      
      plot(x1,y,xlab = "PI-RADS", ylab=name)
      boxplot(y~x2,col = x1,xlab = "BCR Status", ylab=name)
      
      ggplot(data1, aes(x1,y, fill=x1)) + 
        geom_boxplot(alpha=0.3) +
        theme(legend.position="none") +
        scale_fill_jama()+
        #scale_fill_brewer(palette="RdYlBu", direction=-1)+
        #scale_fill_brewer(palette="Dark2", direction=-1)+
        xlab("PI-RADS") + 
        ylab(name)

      ggplot(data1, aes(x2,y, fill=x2)) + 
        geom_boxplot(alpha=0.3) +
        theme(legend.position="none") +
        scale_fill_brewer(palette="Dark2")+
        xlab("BCR Status") + 
        ylab(name)
      
      ggplot(data1, aes(x1, y, color=x1), )+ 
        geom_violin(trim=FALSE)+
        stat_summary(fun.data = mean_se, geom="pointrange", size=.5, color="red") +
        #scale_color_manual(values = c("#313695","#4575B4", "#E69F00", "#FC4E07",  "#A50026"))+
        scale_color_jama()+
        #scale_color_manual(values = wes_palette("Zissou1", n = 5))+
        #scale_color_brewer(palette = "Dark2", direction=-1) +
        #scale_color_viridis(discrete = TRUE,option = "D") +
        xlab("PI-RADS") + 
        ylab(name)
      
      ggplot(data1, aes(x2, y, color=x2))+ 
        geom_violin(trim=FALSE)+
        stat_summary(fun.data = mean_se, geom="pointrange", size=.5, color="red") +
        scale_color_brewer(palette = "Dark2") +
        xlab("BCR Status") + 
        ylab(name)
    
      aggregate(y, list(data1$status_ordinal), function(x) c(mean = mean(x), sd = sd(x)))
      wilcox.test(y ~ data1$status_ordinal, alternative = "two.sided")

      
      
      
      ggplot(data1, aes(x=Gleason, fill=status_ordinal)) + geom_bar(position="fill") +
        scale_fill_jama() +
        xlab("Maximum Gleason Grade Group") +
        ylab("Scaled Count")+
        guides(fill=guide_legend(title="BCR Status"))
      
      library(lattice)
      print(histogram(~PIRADS | EPEscore, data = data1))
            
      ## variations:#####
      
      #ggplot(data1, aes(x=PIRADS, y=adcGLN, fill=PIRADS)) + 
       # geom_violin()
      
      #ggplot(data1, aes(x=status_ordinal, y=adcGLN)) + 
       # geom_violin()
      
      #ggplot(data1, aes(x=PIRADS, y=adcGLN))+ 
       # geom_violin(trim=FALSE)+
       # stat_summary(fun = mean, geom="point", shape=23, size=2, color="red")+
       # xlab("PI-RADS") + 
       # ylab(name)
 
      #ggplot(data1, aes(x=status_ordinal, y=adcGLN), )+
       # geom_violin(trim=FALSE) +
       # stat_summary(fun = mean, geom="point",size=2, color="red")
      

      
      

##############################################################################
################## Normality, Wilcoxon and Chi square   ######################
############################################################################## 
      
    #### Check for normality in order to select test ####
    library(ggpubr)
    ggqqplot(data1$age)
    # points should lie roughly along line
    ggdensity(data1$Age_at_MRI, 
              main = "Density plot",
              xlab = "Variable")
    ggdensity(subset(data1, status == 1)$Rptime, main = "Density plot", xlab = "Variable")
    # should look generally bell shaped
    shapiro.test(data1$Rptime) # checks full set
    shapiro.test(subset(data1, status == 1)$NumLesions) # checks a subset
    # null hypothesis is normal distribution
    # p-value > 0.05 = not significantly different from normal distribution. Can assume the normality if p is high.
    
    ################ Mean, Std dev, Wilcoxon, and t-test #################
    
    aggregate(data1$time, list(data1$status), function(x) c(mean = mean(x), sd = sd(x)))
    aggregate(data1$time, list(data1$status), function(x) c(median = median(x), quantile = quantile(x)))
    
    wilcox.test(data1$Rptime ~ data1$status, alternative = "two.sided")
        # only use t test if data for EACH SUBSET is normally distributed
    
    t.test(data1$Age_at_MRI ~ data1$status, var.equal = TRUE)
    
    
    ############# confusion matrix of variables and Chi square #################
    
    ftable(data1$LNI, data1$status)
    chisq.test(data1$Race, data1$status)
    
    # high p value = differences are purely by chance
    
    
    ############  sensitivity and specificity ######################
    library(caret)
    #threshold=0.5
    #predicted_values<-ifelse(predict(Fiberbits_model_1,type="response")>threshold,1,0)
    predicted_values<-data1$MRI_SVI
    actual_values<-data1$path_SVI
    conf_matrix<-table(predicted_values,actual_values)
    conf_matrix
    sensitivity(conf_matrix)
    specificity(conf_matrix)
    
    
    
    
    
    
    
    
    
######################################################
######## UNIVARIABLE COX-PH ANALYSIS #################
######################################################

    surv.object <- Surv(data1$time, data1$status)
    
    Variable1 <- data1$LNI
  
       
    cph = coxph(surv.object ~ Variable1,data=data1)
    cox.zph(cph) # checks if CoxPH appropriate to use. Should be >0.05
    summary(cph)    
    
    
    
    
    
######################################################
######## MULTIVARIABLE COX-PH ANALYSIS ###############
######################################################
    library("survival", "KMsurv", "Hmisc")
    
    
    #after you find how each variable does in univariate assessment (both continuous and categorical), now we look at how they all do together
    # "forward" selection here is based on our univariate assessment, you can either pick (1) all variables, (2) variables with p<0.1 or p<0.2 in univariate, or (3) known variables
    multiX <- data1[, names(data1) %in% c("time","status",
                                          'PSA',
                                          'PSA_density',
                                          'GleasonSx',
                                          'GleasonRp',
                                          'indexHighB_IQR',
                                          'indexHighB_90th',
                                          'prostHighB_IQR',
                                          'indexADC_10th',
                                          'minADC_IQR',
                                          'maxHighB_90th',
                                          'maxHighB_IQR',
                                          'indexADC_median',
                                          'indexHighB_mean',
                                          'indexADC_mean',
                                          'indexHighB_median',
                                          'prostHighB_90th',
                                          'maxHighB_mean',
                                          'maxHighB_median',
                                          'indexSA',
                                          'indexmax3Ddiameter',
                                          'indexSVR',
                                          'indexSphericity',
                                          'indexElongation',
                                          'indexFlatness',
                                          'path_EPE',
                                          'path_SVI',
                                          'MRI_EPE',
                                          'MRI_SVI',
                                          'margins',
                                          'adcGLN',
                                          't2GLN',
                                          'adcZoneEnt',
                                          'minHighB_uniform',
                                          'highbZoneEnt',
                                          't2ZoneEnt',
                                          'prostHighB_uniform',
                                          'indexHighB_uniform',
                                          'prostHighB_ent',
                                          'highbZonePerc',
                                          't2ZonePerc',
                                          'maxHighB_uniform',
                                          'adcCoarseness',
                                          'T2skewness',
                                          'prostT2_uniform',
                                          'prostT2_ent',
                                          'minHighB_ent',
                                          'minADC_ent',
                                          'P5toProstVol',
                                          'totalVol_P5',
                                          'P5toTotalVol',
                                          'PIRADS.score.of.indexed.lesion',
                                          'totalVol_all',
                                          'P3toTotalVol',
                                          'P4toTotalVol',
                                          'totalVol_P4',
                                          'P4toProstVol',
                                          'NumLesions')]
    
    
    # we only fit multivariable models in patients that have all variables observed, i.e. no "NA" values
    m_cox <- multiX[complete.cases(multiX), ]
    
    
    #fit the initial model
    cox.multi.pre <- coxph(Surv(m_cox$time, m_cox$status) ~  
                           GleasonSx +
                           PSA +
                           PSA_density +
                           indexHighB_IQR +
                           indexHighB_90th +
                           prostHighB_IQR +
                           indexADC_10th +
                           minADC_IQR +
                           maxHighB_90th +
                           maxHighB_IQR +
                           indexADC_median +
                           indexHighB_mean +
                           indexADC_mean +
                           indexHighB_median +
                           prostHighB_90th +
                           maxHighB_mean +
                           maxHighB_median +
                           indexSA +
                           indexmax3Ddiameter +
                           indexSVR +
                           indexSphericity +
                           indexElongation +
                           indexFlatness +
                           MRI_EPE +
                           MRI_SVI +
                           adcGLN +
                           t2GLN +
                           adcZoneEnt +
                           minHighB_uniform +
                           highbZoneEnt +
                           t2ZoneEnt +
                           prostHighB_uniform +
                           indexHighB_uniform +
                           prostHighB_ent +
                           highbZonePerc +
                           t2ZonePerc +
                           maxHighB_uniform +
                           adcCoarseness +
                           T2skewness +
                           prostT2_uniform +
                           prostT2_ent +
                           minHighB_ent +
                           minADC_ent +
                           P5toProstVol +
                           totalVol_P5 +
                           P5toTotalVol +
                           PIRADS.score.of.indexed.lesion +
                           totalVol_all +
                           P3toTotalVol +
                           P4toTotalVol +
                           totalVol_P4 +
                           P4toProstVol +
                           NumLesions
                          , data = m_cox)
    
    summary(cox.multi.pre)
    
    cox.multi.post <- coxph(Surv(m_cox$time, m_cox$status) ~  
                           GleasonRp +
                           PSA +
                           PSA_density +
                           indexHighB_IQR +
                           indexHighB_90th +
                           prostHighB_IQR +
                           indexADC_10th +
                           minADC_IQR +
                           maxHighB_90th +
                           maxHighB_IQR +
                           indexADC_median +
                           indexHighB_mean +
                           indexADC_mean +
                           indexHighB_median +
                           prostHighB_90th +
                           maxHighB_mean +
                           maxHighB_median +
                           indexSA +
                           indexmax3Ddiameter +
                           indexSVR +
                           indexSphericity +
                           indexElongation +
                           indexFlatness +
                           path_EPE +
                           path_SVI +
                           margins +
                           adcGLN +
                           t2GLN +
                           adcZoneEnt +
                           minHighB_uniform +
                           highbZoneEnt +
                           t2ZoneEnt +
                           prostHighB_uniform +
                           indexHighB_uniform +
                           prostHighB_ent +
                           highbZonePerc +
                           t2ZonePerc +
                           maxHighB_uniform +
                           adcCoarseness +
                           T2skewness +
                           prostT2_uniform +
                           prostT2_ent +
                           minHighB_ent +
                           minADC_ent +
                           P5toProstVol +
                           totalVol_P5 +
                           P5toTotalVol +
                           PIRADS.score.of.indexed.lesion +
                           totalVol_all +
                           P3toTotalVol +
                           P4toTotalVol +
                           totalVol_P4 +
                           P4toProstVol +
                           NumLesions
                         , data = m_cox)
    
    summary(cox.multi.post)
    
    cox.multi.combo <- coxph(Surv(m_cox$time, m_cox$status) ~  
                              GleasonRp +
                              GleasonSx+
                              PSA +
                              PSA_density +
                              indexHighB_IQR +
                              indexHighB_90th +
                              prostHighB_IQR +
                              indexADC_10th +
                              minADC_IQR +
                              maxHighB_90th +
                              maxHighB_IQR +
                              indexADC_median +
                              indexHighB_mean +
                              indexADC_mean +
                              indexHighB_median +
                              prostHighB_90th +
                              maxHighB_mean +
                              maxHighB_median +
                              indexSA +
                              indexmax3Ddiameter +
                              indexSVR +
                              indexSphericity +
                              indexElongation +
                              indexFlatness +
                              path_EPE +
                              path_SVI +
                              MRI_EPE + 
                              MRI_SVI +
                              margins +
                              adcGLN +
                              t2GLN +
                              adcZoneEnt +
                              minHighB_uniform +
                              highbZoneEnt +
                              t2ZoneEnt +
                              prostHighB_uniform +
                              indexHighB_uniform +
                              prostHighB_ent +
                              highbZonePerc +
                              t2ZonePerc +
                              maxHighB_uniform +
                              adcCoarseness +
                              T2skewness +
                              prostT2_uniform +
                              prostT2_ent +
                              minHighB_ent +
                              minADC_ent +
                              P5toProstVol +
                              totalVol_P5 +
                              P5toTotalVol +
                              PIRADS.score.of.indexed.lesion +
                              totalVol_all +
                              P3toTotalVol +
                              P4toTotalVol +
                              totalVol_P4 +
                              P4toProstVol +
                              NumLesions
                            , data = m_cox)
    
    summary(cox.multi.combo)
    
    cox.multi.radVol <- coxph(Surv(m_cox$time, m_cox$status) ~  
                               indexHighB_IQR +
                               indexHighB_90th +
                               prostHighB_IQR +
                               indexADC_10th +
                               minADC_IQR +
                               maxHighB_90th +
                               maxHighB_IQR +
                               indexADC_median +
                               indexHighB_mean +
                               indexADC_mean +
                               indexHighB_median +
                               prostHighB_90th +
                               maxHighB_mean +
                               maxHighB_median +
                               indexSA +
                               indexmax3Ddiameter +
                               indexSVR +
                               indexSphericity +
                               indexElongation +
                               indexFlatness +
                               adcGLN +
                               t2GLN +
                               adcZoneEnt +
                               minHighB_uniform +
                               highbZoneEnt +
                               t2ZoneEnt +
                               prostHighB_uniform +
                               indexHighB_uniform +
                               prostHighB_ent +
                               highbZonePerc +
                               t2ZonePerc +
                               maxHighB_uniform +
                               adcCoarseness +
                               T2skewness +
                               prostT2_uniform +
                               prostT2_ent +
                               minHighB_ent +
                               minADC_ent +
                               P5toProstVol +
                               totalVol_P5 +
                               P5toTotalVol +
                               PIRADS.score.of.indexed.lesion +
                               totalVol_all +
                               P3toTotalVol +
                               P4toTotalVol +
                               totalVol_P4 +
                               P4toProstVol +
                               NumLesions
                             , data = m_cox)
    
    summary(cox.multi.radVol)
    
    cox.multi.Vol <- coxph(Surv(m_cox$time, m_cox$status) ~  
                                P5toProstVol +
                                totalVol_P5 +
                                P5toTotalVol +
                                PIRADS.score.of.indexed.lesion +
                                totalVol_all +
                                P3toTotalVol +
                                P4toTotalVol +
                                totalVol_P4 +
                                P4toProstVol +
                                NumLesions
                              , data = m_cox)
    
    summary(cox.multi.Vol)
    
    
    cox.multi.rad <- coxph(Surv(m_cox$time, m_cox$status) ~  
                                indexHighB_IQR +
                                indexHighB_90th +
                                prostHighB_IQR +
                                indexADC_10th +
                                minADC_IQR +
                                maxHighB_90th +
                                maxHighB_IQR +
                                indexADC_median +
                                indexHighB_mean +
                                indexADC_mean +
                                indexHighB_median +
                                prostHighB_90th +
                                maxHighB_mean +
                                maxHighB_median +
                                indexSA +
                                indexmax3Ddiameter +
                                indexSVR +
                                indexSphericity +
                                indexElongation +
                                indexFlatness +
                                adcGLN +
                                t2GLN +
                                adcZoneEnt +
                                minHighB_uniform +
                                highbZoneEnt +
                                t2ZoneEnt +
                                prostHighB_uniform +
                                indexHighB_uniform +
                                prostHighB_ent +
                                highbZonePerc +
                                t2ZonePerc +
                                maxHighB_uniform +
                                adcCoarseness +
                                T2skewness +
                                prostT2_uniform +
                                prostT2_ent +
                                minHighB_ent +
                                minADC_ent
                              , data = m_cox)
    
    summary(cox.multi.rad)
    
    
    cox.multi.3 <- coxph(Surv(m_cox$time, m_cox$status) ~ RP_ISUP + pathEPE , data = m_cox)
    
    
    anova(cox.multi.2, cox.multi.3, test="LRT")
    anova(cox.multi.1, cox.multi.3, test="LRT")   
    
    
    
    #perform AIC selection
    cox.multi.AIC <- step(cox.multi.Vol, direction = "backward")
    
    #perform BIC selection
    n <- dim(m_cox)[1]
    cox.multi.BIC <- step(cox.multi.rad, k=log(n), direction = "backward")
    
    # BIC tends to produce smaller models (fewer variables)
    summary(cox.multi.BIC)
    summary(cox.multi.AIC)    
    
    # please note this multivariable model kept the continuous variables as continuous, if you want to make them categorical based on optimal cut-point create a new variable
    
    ############ LASSO ################
    
    
    y2<-data.matrix(m_cox)[1:nrow(m_cox),1:2]
    foldp<-data.matrix(m_cox)[1:nrow(m_cox),3]
    x2<-data.matrix(m_cox)[1:nrow(m_cox),4:length(m_cox)]
    
    
    
    
    fit2 <- glmnet(x2,y2,standardize = TRUE,family="cox")
    plot(fit2)
    
    
    cvfp <- cv.glmnet(x2, y2, family = "cox", foldid = foldp, type.measure = "C")
    plot(cvfp)
    
    
    
    #find optimal lambda (where CV-error curve hits its minimum)
    lam<-cvfp$lambda.min
    lam
    coef(fit2, s=lam)
    
    #find lambda for regularized model with CV-error within 1 st.dev of min
    lam_se<-cvfp$lambda.1se
    lam_se
    coef(fit2, s=.0966)
    
    
    CF <- as.matrix(coef(cvfp, min(cvfp$lambda[which(cvfp$nzero<=6)])))
    CF[CF!=0,]
    features_i=CF[CF!=0,]
    
    
################################################################    
###################    KM ANALYSIS   ###########################
################################################################
    
    
    surv.object <- Surv(data1$time, data1$status)
    Variable1 <- data1$EPEscore
    
    
    fit <- survfit(surv.object ~ EPEscore, data = data1)
    summary(fit)
    
    #print(fit)
    ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               pval.method = TRUE, 
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c("#2166AC", "orange", "#FF3300","#B2182B"),
               #"#33cc00""purple", "orange", "#FF3300"
               title="MRI/pathology T-stage agreement", xlab="Time [months]", ylab="BCR-free Survival")
    
    res <- pairwise_survdiff(Surv(data1$time, data1$status) ~ data1$FbxISUP,
                             data = data1)
    res
    res <- pairwise_survdiff(Surv(time, status) ~ UCSF_CAPRA_combined_risk,
                             data = data1)
    res
    
    summary(fit)
    
    
    
    
    
#####################################################################  
###################  ROC curves and Youden ##########################
#####################################################################
    
    
    library(ROCR)
    library(precrec)
    library(survival)   
    

    data1$Variable1 <- 
    #cox proportional-hazard
    cph = coxph(Surv(data1$time, data1$status) ~ data1$CAPRA_S_risk,data=data1)
    cox.zph(cph)
    summary(cph)
    
    fit1 <- glm(status ~ PSA_density + MRI_EPE + MRI_SVI + NumLesions  + GleasonSx,  data=m_cox, family = binomial)
    m_cox$prob1=predict(fit1)
    g1 <- roc(status~prob1, data=m_cox, ci=TRUE)
    auc1 = paste("Pre-surgical = ",toString(round(g1$auc, digits = 3)))
    best.coords <- coords(g1, "best", best.method="youden")
    
    fit2 <- glm(status ~ PSA_density + path_EPE + path_SVI + NumLesions + margins + NumLesions + GleasonRp,  data=m_cox, family = binomial)
    m_cox$prob2=predict(fit2)
    g2 <- roc(status~prob2, data=m_cox, ci=TRUE)
    auc2 = paste("Post-surgical  = ",toString(round(g2$auc, digits = 3)))
    best.coords <- coords(g2, "best", best.method="youden")
    
    fit3 <- glm(status ~ PSA_density + MRI_EPE + MRI_SVI + NumLesions + NumLesions + GleasonSx  + path_EPE + path_SVI + margins + GleasonRp , data=m_cox, family = binomial)
    m_cox$prob3=predict(fit3)
    g3 <- roc(status~prob3, data=m_cox, ci=TRUE)
    auc3 = paste("Combined  = ",toString(round(g3$auc, digits = 3)))
    best.coords <- coords(g3, "best", best.method="youden")
    
    fit4<- glm(status ~ indexHighB_90th + prostHighB_IQR + indexADC_10th + maxHighB_IQR + indexADC_median + indexADC_mean + indexHighB_median+ maxHighB_median + indexFlatness + t2GLN + highbZoneEnt + prostHighB_uniform + prostHighB_ent + highbZonePerc, data=m_cox, family = binomial)
    m_cox$prob4=predict(fit4)
    g4 <- roc(status~prob4, data=m_cox, ci=TRUE)
    auc4 = paste("Radiomics = ",toString(round(g4$auc, digits = 3)))
    best.coords <- coords(g4, "best", best.method="youden")
    
    fit5 <- glm(status ~ P5toProstVol + PIRADS.score.of.indexed.lesion,  data=m_cox, family = binomial)
    #fit5 <- glm(status ~ PSA +  NumLesions + PIRADS_cat + Prostate.Volume + NumLesions + SbxISUP,  data=m_cox, family = binomial)
    m_cox$prob5=predict(fit5)
    g5 <- roc(status~prob5, data=m_cox, ci=TRUE)
    auc5 = paste("Volumetrics = ",toString(round(g5$auc, digits = 3)))
    best.coords <- coords(g5, "best", best.method="youden")
    
    
    fit5 <- glm(status ~ 
                  indexHighB_90th + prostHighB_IQR + indexADC_10th + maxHighB_IQR + indexADC_median + indexADC_mean + indexHighB_median+ maxHighB_median + indexFlatness + t2GLN + highbZoneEnt + prostHighB_uniform + prostHighB_ent + highbZonePerc
                  +PSA_density + MRI_EPE + MRI_SVI + NumLesions + GleasonSx + 
                  P5toProstVol + PIRADS.score.of.indexed.lesion,  data=m_cox, family = binomial)
    #fit5 <- glm(status ~ PSA +  NumLesions + PIRADS_cat + Prostate.Volume + NumLesions + SbxISUP,  data=m_cox, family = binomial)
    m_cox$prob5=predict(fit5)
    g5 <- roc(status~prob5, data=m_cox, ci=TRUE)
    auc5 = paste("Pre-surgical + volumetrics = ",toString(round(g5$auc, digits = 3)))
    best.coords <- coords(g5, "best", best.method="youden")
    
    # fit2 <- glm(status ~ PSA_density + path_EPE + path_SVI + NumLesions + margins + NumLesions + GleasonRp,  data=m_cox, family = binomial)
    fit6 <- glm(status ~ indexHighB_90th + prostHighB_IQR + indexADC_10th + maxHighB_IQR + indexADC_median + indexADC_mean + indexHighB_median+ maxHighB_median + indexFlatness + t2GLN + highbZoneEnt + prostHighB_uniform + prostHighB_ent + highbZonePerc
                +PSA_density + MRI_EPE + MRI_SVI + NumLesions + GleasonSx 
                  ,  data=m_cox, family = binomial)
                  
                  
                  #indexHighB_90th + prostHighB_IQR + indexADC_10th + maxHighB_IQR + indexADC_median + indexADC_mean + indexHighB_median+ maxHighB_median + indexFlatness + t2GLN + highbZoneEnt + prostHighB_uniform + prostHighB_ent + highbZonePerc
                  #PSA_density + path_EPE + path_SVI + NumLesions + margins + NumLesions + GleasonRp  + 
                  #P5toProstVol + PIRADS.score.of.indexed.lesion,  data=m_cox, family = binomial)
    m_cox$prob6=predict(fit6)
    g6 <- roc(status~prob6, data=m_cox, ci=TRUE)
    auc6 = paste("Post-surgical + volumetrics = ",toString(round(g6$auc, digits = 3)))
    best.coords <- coords(g6, "best", best.method="youden")
    
    
    # fit3 <- glm(status ~ PSA_density + MRI_EPE + MRI_SVI + NumLesions + NumLesions + GleasonSx  + path_EPE + path_SVI + margins + GleasonRp , data=m_cox, family = binomial)
    fit7 <- glm(status ~ 
                  #indexHighB_90th + prostHighB_IQR + indexADC_10th + maxHighB_IQR + indexADC_median + indexADC_mean + indexHighB_median+ maxHighB_median + indexFlatness + t2GLN + highbZoneEnt + prostHighB_uniform + prostHighB_ent + highbZonePerc
                  PSA_density + MRI_EPE + MRI_SVI + NumLesions + NumLesions + GleasonSx  + path_EPE + path_SVI + margins + GleasonRp  + 
                  P5toProstVol + PIRADS.score.of.indexed.lesion, data=m_cox, family = binomial)
    m_cox$prob7=predict(fit7)
    g7 <- roc(status~prob7, data=m_cox, ci=TRUE)
    auc7 = paste("Combined + volumetrics = ",toString(round(g7$auc, digits = 3)))
    best.coords <- coords(g7, "best", best.method="youden")
    
    fit8<- glm(status ~ indexHighB_90th + prostHighB_IQR + indexADC_10th + maxHighB_IQR + indexADC_median + indexADC_mean + indexHighB_median+ maxHighB_median + indexFlatness + t2GLN + highbZoneEnt + prostHighB_uniform + prostHighB_ent + highbZonePerc
               + PSA_density + MRI_EPE + MRI_SVI + NumLesions + GleasonSx + P5toTotalVol, data=m_cox, family = binomial)
    m_cox$prob8=predict(fit8)
    g8 <- roc(status~prob8, data=m_cox, ci=TRUE)
    auc8 = paste("8. rad + vol = ",toString(round(g8$auc, digits = 3)))
    best.coords <- coords(g8, "best", best.method="youden")
 
    
#######################################################################    
    fit1 <- glm(status ~ PSA_density + MRI_EPE + MRI_SVI,  data=m_cox, family = binomial)
    m_cox$prob1=predict(fit1)
    g1 <- roc(status~prob1, data=m_cox, ci=TRUE)
    auc1 = paste("Pre-surcial = ",toString(round(g1$auc, digits = 3)))
    best.coords <- coords(g1, "best", best.method="youden")
    
    fit2 <- glm(status ~ PSA_density + path_EPE + path_SVI + margins,  data=m_cox, family = binomial)
    m_cox$prob2=predict(fit2)
    g2 <- roc(status~prob2, data=m_cox, ci=TRUE)
    auc2 = paste("Post-surgical  = ",toString(round(g2$auc, digits = 3)))
    best.coords <- coords(g2, "best", best.method="youden")
    
    fit3 <- glm(status ~ PSA_density + MRI_EPE + MRI_SVI + path_EPE + path_SVI + margins , data=m_cox, family = binomial)
    m_cox$prob3=predict(fit3)
    g3 <- roc(status~prob3, data=m_cox, ci=TRUE)
    auc3 = paste("Combined = ",toString(round(g3$auc, digits = 3)))
    best.coords <- coords(g3, "best", best.method="youden")
    
    fit4<- glm(status ~ indexHighB_90th + minADC_IQR + indexHighB_median + indexSphericity + indexFlatness + highbZoneEnt + prostHighB_uniform + indexHighB_uniform + highbZonePerc, data=m_cox, family = binomial)
    m_cox$prob4=predict(fit4)
    g4 <- roc(status~prob4, data=m_cox, ci=TRUE)
    auc4 = paste("Radiomics = ",toString(round(g4$auc, digits = 3)))
    best.coords <- coords(g4, "best", best.method="youden")
    
    fit5 <- glm(status ~ P5toProstVol + PIRADS.score.of.indexed.lesion,  data=m_cox, family = binomial)
    #fit5 <- glm(status ~ PSA +  NumLesions + PIRADS_cat + Prostate.Volume + NumLesions + SbxISUP,  data=m_cox, family = binomial)
    m_cox$prob5=predict(fit5)
    g5 <- roc(status~prob5, data=m_cox, ci=TRUE)
    auc5 = paste("Volumetrics = ",toString(round(g5$auc, digits = 3)))
    best.coords <- coords(g5, "best", best.method="youden")   
    
    
    
    
    
    fit5 <- glm(status ~ 
                  #indexHighB_90th + minADC_IQR + indexHighB_median + indexSphericity + indexFlatness + highbZoneEnt + prostHighB_uniform + indexHighB_uniform + highbZonePerc
                + PSA_density + MRI_EPE + MRI_SVI + P5toProstVol + PIRADS.score.of.indexed.lesion,  data=m_cox, family = binomial)
    #fit5 <- glm(status ~ PSA +  NumLesions + PIRADS_cat + Prostate.Volume + NumLesions + SbxISUP,  data=m_cox, family = binomial)
    m_cox$prob5=predict(fit5)
    g5 <- roc(status~prob5, data=m_cox, ci=TRUE)
    auc5 = paste("Pre-surgical + radiomics = ",toString(round(g5$auc, digits = 3)))
    best.coords <- coords(g5, "best", best.method="youden")
    
    # fit2 <- glm(status ~ PSA_density + path_EPE + path_SVI + NumLesions + margins + NumLesions + GleasonRp,  data=m_cox, family = binomial)
    fit6 <- glm(status ~ 
                  #indexHighB_90th + minADC_IQR + indexHighB_median + indexSphericity + indexFlatness + highbZoneEnt + prostHighB_uniform + indexHighB_uniform + highbZonePerc
                + PSA_density + path_EPE + path_SVI  + margins+ P5toTotalVol + P5toProstVol + PIRADS.score.of.indexed.lesion,  data=m_cox, family = binomial)
    m_cox$prob6=predict(fit6)
    g6 <- roc(status~prob6, data=m_cox, ci=TRUE)
    auc6 = paste("Post-surgical + radiomics + volumetrics = ",toString(round(g6$auc, digits = 3)))
    best.coords <- coords(g6, "best", best.method="youden")
    
    
    # fit3 <- glm(status ~ PSA_density + MRI_EPE + MRI_SVI + NumLesions + NumLesions + GleasonSx  + path_EPE + path_SVI + margins + GleasonRp , data=m_cox, family = binomial)
    fit7 <- glm(status ~ 
                  #indexHighB_90th + minADC_IQR + indexHighB_median + indexSphericity + indexFlatness + highbZoneEnt + prostHighB_uniform + indexHighB_uniform + highbZonePerc
                + PSA_density + MRI_EPE + MRI_SVI + path_EPE + path_SVI + margins + P5toTotalVol + P5toProstVol + PIRADS.score.of.indexed.lesion, data=m_cox, family = binomial)
    m_cox$prob7=predict(fit7)
    g7 <- roc(status~prob7, data=m_cox, ci=TRUE)
    auc7 = paste("Combined + radiomics + volumetrics = ",toString(round(g7$auc, digits = 3)))
    best.coords <- coords(g7, "best", best.method="youden")
    
    fit8<- glm(status ~ indexHighB_90th + minADC_IQR + indexHighB_median + indexSphericity + indexFlatness + highbZoneEnt + prostHighB_uniform + indexHighB_uniform + highbZonePerc
               + PSA_density + MRI_EPE + MRI_SVI + P5toTotalVol, data=m_cox, family = binomial)
    m_cox$prob8=predict(fit8)
    g8 <- roc(status~prob8, data=m_cox, ci=TRUE)
    auc8 = paste("8. rad + vol = ",toString(round(g8$auc, digits = 3)))
    best.coords <- coords(g8, "best", best.method="youden")
    
 
    
    plot(g1, col = 'purple')
    text(0.65,0.25,'--',  col = 'purple')
    text(0.32,.25,auc1,  col = 'purple')
    
    plot(g2, add=TRUE, col='#2166AC')
    text(0.65,0.20,'--', col = '#2166AC')
    text(0.30,0.20,auc2,col='#2166AC')
    
    plot(g3, add=TRUE, col='orange')
    text(0.65,0.15,'--', col = 'orange')
    text(0.35,0.15,auc3,col='orange')
    
    plot(g5, add=TRUE, col='#FF3300')
    text(0.65,0.20,'--', col = '#FF3300')
    text(0.37,0.20,auc5,col='#FF3300')
    
    plot(g4, add=TRUE, col='#B2182B')
    text(0.65,0.15,'--', col = '#B2182B')
    text(0.38,0.15,auc4,col='#B2182B')
    
    
    
    
    
    plot(g5, col = 'purple')
    text(0.85,0.25,'--',  col = 'purple')
    text(0.30,.25,auc5,  col = 'purple')
    
    plot(g6, add=TRUE, col='#2166AC')
    text(0.85,0.20,'--', col = '#2166AC')
    text(0.29,0.20,auc6,col='#2166AC')
    
    plot(g7, add=TRUE, col='orange')
    text(0.85,0.15,'--', col = 'orange')
    text(0.33,0.15,auc7,col='orange')
    
       
    
    plot(g5, col = 'purple')
    text(0.85,0.45,'--',  col = 'purple')
    text(0.20,.45,auc5,  col = 'purple')
    
    plot(g6, add=TRUE, col='#2166AC')
    text(0.85,0.30,'--', col = '#2166AC')
    text(0.19,0.30,auc6,col='#2166AC')
    
    plot(g7, add=TRUE, col='orange')
    text(0.85,0.15,'--', col = 'orange')
    text(0.24,0.15,auc7,col='orange')
    
    
    
    plot(g1, col = 'purple')
    text(0.9,0.35,'--',  col = 'purple')
    text(0.48,.35,auc1,  col = 'purple')
    
    plot(g2, add=TRUE, col='#2166AC')
    text(0.9,0.30,'--', col = '#2166AC')
    text(0.54,0.30,auc2,col='#2166AC')
    
    plot(g3, add=TRUE, col='#FF3300')
    text(0.9,0.25,'--', col = '#FF3300')
    text(0.43,0.25,auc3,col='#FF3300')
    
    plot(g4, add=TRUE, col='#B2182B')
    text(0.9,0.20,'--', col = '#B2182B')
    text(0.43,0.20,auc4,col='#B2182B')
    
    plot(g5, add=TRUE, col='#2166AC')
    text(0.9,0.15,'--', col = '#2166AC')
    text(0.54,0.15,auc5,col='#2166AC')
    
    plot(g6, add=TRUE, col='#FF3300')
    text(0.9,0.10,'--', col = '#FF3300')
    text(0.43,0.10,auc6,col='#FF3300')
    
    plot(g7, add=TRUE, col='#B2182B')
    text(0.9,0.05,'--', col = '#B2182B')
    text(0.43,0.05,auc7,col='#B2182B')    
    
    roc.test(g1, g2)
    roc.test(g1, g3)
    roc.test(g1, g4)
    roc.test(g1, g5)
    roc.test(g1, g6)
    roc.test(g1, g7)
    roc.test(g2, g3)
    roc.test(g2, g4)
    roc.test(g2, g5)
    roc.test(g2, g6)
    roc.test(g2, g7)
    roc.test(g3, g4)
    roc.test(g3, g5)
    roc.test(g3, g6)
    roc.test(g3, g7)
    roc.test(g4, g5)
    roc.test(g4, g6)
    roc.test(g4, g7)
    roc.test(g5, g6)
    roc.test(g5, g7)
    roc.test(g6, g7)
    
    # palette = c("#2166AC", "#B2182B", "orange", "#FF3300")
    #"#33cc00""purple", "orange", "#FF3300"
    
    
    ############## EXAMPLE AUC AND YOUDEN INDEX CALCULATIONS #################
    
    # write a function to calculate optimal cut-point based on Youden Index
    opt.cut = function(in.cut, in.sens, in.spec){
      d = in.sens + in.spec -1
      ind = which(d == max(d))
      ind = ind[1]
      opt.out <- c(cutoff = in.cut[[ind]])
      return(opt.out)
    }
    

    # remember that censor = 0 is not observed (no recurrence) and censor=1 is obersved (recurrence)
    pred <- prediction(data1$Vol_P4, data1$status) #prediction outcome from ADC data
    perf2 <- performance(pred, "auc") #calculate AUC from Variable1
    est.auc = perf2@y.values[[1]]
    
    # get acc, tpr, fpr plot information so we can see how performance changes at each Variable1 value
    perf1 <- performance(pred, "tpr", "fpr")
    perf3 <- performance(pred, "acc")
    # sensitivity and specificity values
    sensitivity = perf1@y.values[[1]]
    specificity = 1-perf1@x.values[[1]]
    # cutoffs (i.e. Variable1 values)
    cutoffs = perf3@x.values[[1]]
    # calculate optimal cut point based on our function
    est.cut = opt.cut(cutoffs, sensitivity, specificity)
    
    est.auc
    est.cut        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 #########  OLD CODE! MAY OR MAY NOT BE USEFUL! ################################   
    
    
    
    

######################################################
####### CATEGORICAL KM and COX-PH ANALYSIS ###########
######################################################
 
    surv.object <- Surv(data1$time, data1$status)
    Variable1 <- data1$EPEscore
        

    #code from Stephanie which includes legend
    # myfits.v4 <- survfit(surv.object~Variable1, data = data1)
    # myfits.v4.LR = survdiff(surv.object~Variable1,data=data1)
    # myfits.v4.p = 1-pchisq(myfits.v4.LR$chisq, 1)
    # plot(myfits.v4, col=c(1,2,4), lwd=2,main="P5 Relative Volume - Margins Included", xlab="Time", ylab="Progression Free Survival", cex.lab=2,cex.axis=2, mark.time = TRUE, mark =3, cex=2, cex.main = 2)

    # # to add a legend we need to know the category names
    # # in our fictious example it has 2 groups
    # legend("bottomright", c(paste("v4 low, N=", myfits.v4$n[1], sep=""),paste("v4 high, N=", myfits.v4$n[2], sep="")),lwd=2, col=c(1,2,4), cex=1.75)
    # #lets put p-value on the plot
    # text(5,0.05,paste("P =",signif(as.numeric(myfits.v4.p),2)), cex=1.2)
    
    data1$Variable1 <- cut(data1$P5vol2wp, 
                           breaks=c(-Inf,0.026,Inf), right = TRUE, 
                           labels = c('low relative volume', 'high relative volume'))
    data1$SVI_agreement <- data1$Tstage_agreement
    
    fit <- survfit(surv.object ~ Tstage_agreement, data = data1)
    #print(fit)
    ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               pval.method = TRUE, 
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c("#2166AC", "orange", "#FF3300","#B2182B"),
               #"#33cc00""purple", "orange", "#FF3300"
               title="MRI/pathology T-stage agreement", xlab="Time [months]", ylab="BCR-free Survival")
    
    res <- pairwise_survdiff(Surv(data1$time, data1$status) ~ data1$FbxISUP,
                             data = data1)
    res
    res <- pairwise_survdiff(Surv(time, status) ~ UCSF_CAPRA_combined_risk,
                             data = data1)
    res
    
    summary(fit)

    surv.object <- Surv(data1$time, data1$status)

    data1$Variable1 <- scale(data1$P3vol2wp)

    
    
    
    
    
    data1$Variable1 <- data1$path_SVI
    #cox proportional-hazard
    cph = coxph(Surv(data1$time, data1$status) ~ data1$Variable1,data=data1)
    cox.zph(cph)
    summary(cph)
### this is all added temporarily #####
    opt.cut = function(in.cut, in.sens, in.spec){
      d = in.sens + in.spec -1
      ind = which(d == max(d))
      ind = ind[1]
      opt.out <- c(cutoff = in.cut[[ind]])
      return(opt.out)
    }
    surv.object <- Surv(data1$time, data1$status)
    

    
    text = "Minimum single lesion ADC entropy\n[scaled cutoff = 0.5283]"
    data1$Variable1 <- data1$Min.ADC.Entropy
    fit1 <- glm(status ~ Variable1, data=data1, family = binomial)
    data1$prob1=predict(fit1)
    g1 <- roc(status~prob1, data=data1, ci=TRUE)
    auc1 = paste("Post-operative AUC = ",toString(round(g1$auc, digits = 3)))
    best.coords <- coords(g1, "best", best.method="youden")
    pred <- prediction(data1$Variable1, data1$status) #prediction outcome from ADC data
    perf2 <- performance(pred, "auc") #calculate AUC from Variable1
    est.auc = perf2@y.values[[1]]
    perf1 <- performance(pred, "tpr", "fpr")
    perf3 <- performance(pred, "acc")
    sensitivity = perf1@y.values[[1]]
    specificity = 1-perf1@x.values[[1]]
    cutoffs = perf3@x.values[[1]]
    est.cut = opt.cut(cutoffs, sensitivity, specificity)
    g1$auc
    g1$ci
    est.cut
    best.coords$sensitivity
    best.coords$specificity
    threshold<-as.numeric(est.cut)
    data1$Variable2<-cut(data1$Variable1, 
                         breaks=c(-Inf,threshold,Inf), right = TRUE, 
                         labels = c('below cutoff', 'above cutoff'))
    fit <- survfit(surv.object ~ Variable2, data = data1)
    summary(fit)
    
    ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               pval.method = TRUE, 
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c("#2166AC", "#B2182B", "orange"),
               #"#33cc00""purple", "orange", "#FF3300"
               title=text, xlab="Time [months]", ylab="BCR-free Survival")
    
    
    plot(g1)

    ##########################################################    
    
    
    fit <- survfit(surv.object ~ Variable1, data = data1)
    summary(fit)
    
    #adjusted cox proportional-hazard
    cph = coxph(surv.object ~ PSA + path_EPE + path_SVI + RP_ISUP + margins,data=data1)
    cox.zph(cph)
    summary(cph)
    
    cph = coxph(surv.object ~ PSA + MRI_EPE + MRI_SVI + CbxISUP,data=data1)
    cox.zph(cph)
    summary(cph)
    
    cph = coxph(surv.object ~ PSA + MRI_EPE + MRI_SVI + SbxISUP,data=data1)
    cox.zph(cph)
    summary(cph)
    
    
    EPEscore + SbxISUP + FbxISUP+ CbxISUP+ RP_ISUP + PSA + PSA_cat + MRI_EPE + MRI_SVI 
    + path_EPE + path_SVI + margins + PSA_density + PIRADS + Min..Mean.ADC + Min..Median.ADC
    + Min.ADC.Entropy + WP.T2.Entropy + PIRADS.5.relative.volume + Total.relative.lesion.volume 
    + Indexed.ADC.Mean + Indexed.ADC.Median + Indexed.ADC.10th.Percentile + Indexed.Lesion.Volume
    + Indexed.Relative.Lesion.Volume + Surface.to.Volume.Ratio + Sphericity + Maximum.3D.Diameter
    
    
    
    
######################################################
######## CONTINUOUS VARIABLES COX-PH ANALYSIS ########
######################################################
    # as a short-cut, here is a neat way to print all continous variable data into a table
    
    # example for selecting variables by name
    # I only select continuous variables for this
    #  x_cox <- inData[, names(inData) %in% c("Time","Outcome","Variable1","Variable2","Variable3")]
    
    x_cox <- data1[, names(data1) %in% c("time","status","age")]
    
    surv.object <- Surv(x_cox$time, x_cox$status)
    
    #do univariate regression on each variable individually
    cox.univ <- vector(mode="list", length=(dim(x_cox)[2]-1))
    cox.coef <- array(-1, dim=c(dim(x_cox)[2],8))
    colnames(cox.coef) <- c("name", "p_val", "coef", "exp(coef)","SE_coef","low_95CI","up_95CI", "rsq");
    
    # since my first two variables are Time and Outcome, I start at 3
    for (i in 3:dim(x_cox)[2]) {
      # I scale the data so that the HR are comparable
      cox.univ[[i]] <- coxph(surv.object ~ scale(x_cox[,i]))
      cox.coef[i-2,1] <- colnames(x_cox)[i] #variable name
      cox.coef[i-2,2] <- summary(cox.univ[[i]])$coefficients[1,5] #p-val
      cox.coef[i-2,3] <- summary(cox.univ[[i]])$coefficients[1,1] #the coefficient
      cox.coef[i-2,4] <- summary(cox.univ[[i]])$coefficients[1,2] #the exponential of the coefficient
      cox.coef[i-2,5] <- summary(cox.univ[[i]])$coefficients[1,3] #SE of coefficient
      cox.coef[i-2,6] <- summary(cox.univ[[i]])$conf.int[1,3] #lower 95% CI
      cox.coef[i-2,7] <- summary(cox.univ[[i]])$conf.int[1,4] #upper 95% CI
      cox.coef[i-2,8] <- summary(cox.univ[[i]])$rsq[1]
    }
    
    write.table(cox.coef, file = "clipboard", sep="\t", row.names = FALSE, col.names = TRUE) #copy to clipboard
    
    
    
    
#########################################################
########## NUMERICAL FIND OPTIMAL CUTOFF ################
#########################################################
    #finding optimal cutoff for imaging data
    # you can either make a new dataframe with only significant variables or you can just run them all
    #example for selecting variables by name
    sigVars <- x_cox[, names(x_cox) %in% c("time","status","PSA","age")]
    surv.object <- Surv(sigVars$time, sigVars$status)
    # PLOTTING CATEGORICAL DATA, make sure this folder exists
    # the for loop will save plots for all variables
    directory2save = "T:\MIP\Katie_Merriman\Project1Data\UPDATEDdata"
    
    # since my first two variables are Time and Outcome, I start at 3
    for (i in 3:dim(sigVars)[2]){
      optResults <- matrix(NA,dim(surv.object)[1],2)
      var <- sigVars[,i]
      var.sort <- var[order(var)]
      surv.sort <- surv.object[order(var)]  
      #go through values and find optimal cut-off based on survdiff fit
      for (j in 6:(length(var.sort)-6)){
        V <- var.sort
        optResults[j,1] <- var.sort[j]
        V[which(var.sort > var.sort[j])] <- 1
        V[which(var.sort <= var.sort[j])] <- 0
        M = survdiff(surv.sort~V)
        optResults[j,2] = 1-pchisq(M$chisq, 1)
      }
      #find best cutoff-model
      ind <- which.min(optResults[,2])
      Vf <- var.sort
      Vf[which(var.sort > optResults[ind,1])] <- 1
      Vf[which(var.sort <= optResults[ind,1])] <- 0
      M2 <- survfit(surv.sort~Vf)
      #save image
      tiff(filename=paste(directory2save, colnames(sigVars)[i],"_opt.tiff"), units="in", width=6, height = 6, pointsize=12, res=300)
      plot(M2, col=c(1,2), lty=1:2,lwd=2,main=colnames(sigVars)[i], xlab="Years", ylab="Progression-Free Survival", cex.lab=1.5,cex.axis=1.3)
      legend("topright", c(paste("\u2264",signif(optResults[ind,1],2),", N=",M2$n[1] , sep=""),paste(">",signif(optResults[ind,1],2), ", N=",M2$n[2] , sep="")),lty=c(1,2),lwd=2,  col=c(1,2), cex=1.2)
      
      #print p-value from cox!
      coxM <- coxph(surv.object ~ sigVars[,i])
      p = summary(coxM)$coefficients[1,5]
      text(5,0.05,paste("P =",signif(as.numeric(p),2)), cex=1.2) 
      dev.off()
    }
    
    
    ##########################################################################
    ############## EXAMPLE AUC AND YOUDEN INDEX CALCULATIONS #################
    ##########################################################################
    library(ROCR)
    library(precrec)
    library(survival)
    
    #### write a function to calculate optimal cut-point based on Youden Index
    opt.cut = function(in.cut, in.sens, in.spec){
      d = in.sens + in.spec -1
      ind = which(d == max(d))
      ind = ind[1]
      opt.out <- c(cutoff = in.cut[[ind]])
      return(opt.out)
    }
    
    g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData.csv"
    inData <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)
    inData$Censor = inData$Status
    inData$Variable1 = inData$Vol_P4
    
    # remember that censor = 0 is not observed (no recurrence) and censor=1 is obersved (recurrence)
    pred <- prediction(inData$Vol_P4, inData$Censor) #prediction outcome from ADC data
    perf2 <- performance(pred, "auc") #calculate AUC from Variable1
    est.auc = perf2@y.values[[1]]
    
    # get acc, tpr, fpr plot information so we can see how performance changes at each Variable1 value
    perf1 <- performance(pred, "tpr", "fpr")
    perf3 <- performance(pred, "acc")
    # sensitivity and specificity values
    sensitivity = perf1@y.values[[1]]
    specificity = 1-perf1@x.values[[1]]
    # cutoffs (i.e. Variable1 values)
    cutoffs = perf3@x.values[[1]]
    # calculate optimal cut point based on our function
    est.cut = opt.cut(cutoffs, sensitivity, specificity)
    
    est.auc
    est.cut    

###############################################################################    
################# Model comparison for grouped data  ##########################
###############################################################################
    
   
    multiX <- data1[, names(data1) %in% c("time","status","fold", "EPEscore", "SbxISUP", "FbxISUP", "CbxISUP", "RP_ISUP", "PSA_cat", "MRI_EPE", "MRI_SVI",                       
                                          "path_EPE", "path_SVI", "margins", "PSA_density", "PIRADS", "Min..Mean.ADC", "Min..Median.ADC",
                                          "Min.ADC.Entropy", "WP.T2.Entropy", "PIRADS.5.relative.volume", "Total.relative.lesion.volume",
                                          "Indexed.ADC.Mean", "Indexed.ADC.Median", "Indexed.ADC.10th.Percentile", "Indexed.Lesion.Volume",
                                          "Indexed.Relative.Lesion.Volume", "Surface.to.Volume.Ratio", "Sphericity", "Maximum.3D.Diameter")]
    
    
    
    
    
    
    
###################  ROC curves and Youden ##########################    
    
    
    
    
    
    
    
    
    data1$Variable1 <- 
    #cox proportional-hazard
    cph = coxph(Surv(data1$time, data1$status) ~ data1$CAPRA_S_risk,data=data1)
    cox.zph(cph)
    summary(cph)
    
    fit1 <- glm(status ~ PSA + MRI_EPE + MRI_SVI + NumLesions + PIRADS_cat + Prostate.Volume + NumLesions + SbxISUP,  data=m_cox, family = binomial)
    m_cox$prob1=predict(fit1)
    g1 <- roc(status~prob1, data=m_cox, ci=TRUE)
    auc1 = paste("UCSF-CAPRA AUC = ",toString(round(g1$auc, digits = 3)))
    best.coords <- coords(g1, "best", best.method="youden")
    
    fit2 <- glm(status ~ PSA + path_EPE + path_SVI + NumLesions + margins + Prostate.Volume + NumLesions + RP_ISUP,  data=m_cox, family = binomial)
    m_cox$prob2=predict(fit2)
    g2 <- roc(status~prob2, data=m_cox, ci=TRUE)
    auc2 = paste("CAPRA-S AUC = ",toString(round(g2$auc, digits = 3)))
    best.coords <- coords(g2, "best", best.method="youden")
    
    
    fit3 <- glm(status ~ PSA + MRI_EPE + MRI_SVI + NumLesions + PIRADS_cat + Prostate.Volume + NumLesions + SbxISUP  + path_EPE + path_SVI + margins + RP_ISUP , data=m_cox, family = binomial)
    m_cox$prob3=predict(fit3)
    g3 <- roc(status~prob3, data=m_cox, ci=TRUE)
    auc3 = paste("Modified CAPRA-S AUC = ",toString(round(g3$auc, digits = 3)))
    best.coords <- coords(g3, "best", best.method="youden")
    
    plot(g1, col = 'purple')
    text(0.9,0.15,'--',  col = 'purple')
    text(0.48,.15,auc1,  col = 'purple')
    
    plot(g2, add=TRUE, col='#2166AC')
    text(0.9,0.10,'--', col = '#2166AC')
    text(0.54,0.10,auc2,col='#2166AC')
    
    plot(g3, add=TRUE, col='#FF3300')
    text(0.9,0.05,'--', col = '#FF3300')
    text(0.43,0.05,auc3,col='#FF3300')
    
    
    roc.test(g1, g2)
    roc.test(g2, g3)
    roc.test(g1, g3)
    
    palette = c("#2166AC", "#B2182B", "orange", "#FF3300"),
    #"#33cc00""purple", "orange", "#FF3300"

    
    
    
    fit1 <- glm(status ~ PSA.density + pathEPE + pathSVI + RP_ISUP+ EPEscore, data=m_cox, family = binomial)
    m_cox$prob1=predict(fit1)
    g1 <- roc(status~prob1, data=m_cox, ci=TRUE)
    auc1 = paste("Post-operative, no radiomics\nAUC = ",toString(round(g1$auc, digits = 3)))
    auc11 = paste("Post-operative, no radiomics")
    auc12 = paste("AUC = ",toString(round(g1$auc, digits = 3)))
    best.coords <- coords(g1, "best", best.method="youden")
    
    fit2 <- glm(status ~ PSA.density + pathEPE + pathSVI + RP_ISUP + EPEscore 
    #            + indexVol2total, data=m_cox, family = binomial)
    + P5vol2Total + indexVol2total+ Indexed.ADC.Median + WP.T2.Entropy + Max.T2.Entropy + Flatness + Sphericity, data=m_cox, family = binomial)
    m_cox$prob2=predict(fit2)
    g2 <- roc(status~prob2, data=m_cox, ci=TRUE)
    auc2 = paste("Post-operative + radiomics\nAUC = ",toString(round(g2$auc, digits = 3)))
    auc21 = paste("Post-operative + radiomics")
    auc22 = paste("AUC = ",toString(round(g2$auc, digits = 3)))
    best.coords <- coords(g2, "best", best.method="youden")
    
    fit3 <- glm(status ~ PSA.density + MRI_EPE + MRI_SVI + CbxISUP, data=m_cox, family = binomial)
    m_cox$prob3=predict(fit3)
    g3 <- roc(status~prob3, data=m_cox, ci=TRUE)
    auc3 = paste("Pre-operative - standard EPE, no radiomics\nAUC = ",toString(round(g3$auc, digits = 3)))
    auc31 = paste("Pre-operative - standard EPE, no radiomics")
    auc32 = paste("AUC = ",toString(round(g3$auc, digits = 3)))
    best.coords <- coords(g3, "best", best.method="youden")
    
    fit4 <- glm(status ~ PSA.density + MRI_EPE + MRI_SVI + CbxISUP 
    #            + indexVol2total, data=m_cox, family = binomial)
    + P5vol2Total + indexVol2total+ Indexed.ADC.Median + WP.T2.Entropy + Max.T2.Entropy + Flatness + Sphericity, data=m_cox, family = binomial)
    m_cox$prob4=predict(fit4)
    g4 <- roc(status~prob4, data=m_cox, ci=TRUE)
    auc4 = paste("Pre-operative - standard EPE,+ radiomics\nAUC = ",toString(round(g4$auc, digits = 3)))
    auc41 = paste("Pre-operative - standard EPE,+ radiomics")
    auc42 = paste("AUC = ",toString(round(g4$auc, digits = 3)))
    best.coords <- coords(g4, "best", best.method="youden")
    
    
    fit5 <- glm(status ~ PSA.density + EPEscore + MRI_SVI + CbxISUP, data=m_cox, family = binomial)
    m_cox$prob5=predict(fit5)
    g5 <- roc(status~prob5, data=m_cox, ci=TRUE)
    auc5 = paste("Pre-operative - NCI-EPE grade, no radiomics\nAUC = ",toString(round(g5$auc, digits = 3)))
    auc51 = paste("Pre-operative - EPE grade, no radiomics")
    auc52 = paste("AUC = ",toString(round(g5$auc, digits = 3)))
    best.coords <- coords(g5, "best", best.method="youden")
    
    fit6 <- glm(status ~ PSA.density + EPEscore + MRI_SVI + CbxISUP 
    #            + indexVol2total, data=m_cox, family = binomial)
    + P5vol2Total + indexVol2total+ Indexed.ADC.Median + WP.T2.Entropy + Max.T2.Entropy + Flatness + Sphericity, data=m_cox, family = binomial)
    m_cox$prob6=predict(fit6)
    g6 <- roc(status~prob6, data=m_cox, ci=TRUE)
    auc6 = paste("Pre-operative - NCI-EPE grade + radiomics\nAUC = ",toString(round(g6$auc, digits = 3)))
    auc61 = paste("Pre-operative - EPE grade + radiomics")
    auc62 = paste("AUC = ",toString(round(g6$auc, digits = 3)))
    best.coords <- coords(g6, "best", best.method="youden")
   
  
    plot(g1, col = 'black')
    text(0.68,0.55,'__', cex = 0.9, col = 'black')
    text(0.18,0.55,auc1, cex = 0.9, col = 'black')
    
    plot(g2A, add=TRUE, col='#9D9D9D')
    text(0.68,0.45,'__', cex = 0.9,col = '#9D9D9D')
    text(0.18,0.45,auc2A, cex = 0.9,col='#9D9D9D')
    
    plot(g2A, add=TRUE, col='#CB1515')
    text(0.68,0.35,'__', cex = 0.9,col = '#CB1515')
    text(0.18,0.35,auc2B, cex = 0.9,col='#CB1515')
    
    plot(g3, add=TRUE, col='#CB1515')
    text(0.68,0.35,'__', cex = 0.9,col = '#CB1515')
    text(0.18,0.35,auc3, cex = 0.9,col='#CB1515')
    
    plot(g4A, add=TRUE, col='#F38026')
    text(0.68,0.25,'__', cex = 0.9,col = '#F38026')
    text(0.18,0.25,auc4A, cex = 0.9,col='#F38026')
    
    plot(g4B, add=TRUE, col='#F38026')
    text(0.68,0.25,'__', cex = 0.9,col = '#F38026')
    text(0.18,0.25,auc4B, cex = 0.9,col='#F38026')
    
    
    plot(g5, add=TRUE, col='#1515CB')
    text(0.68,0.15,'__', cex = 0.9,col = '#1515CB')
    text(0.18,0.15,auc5, cex = 0.9,col='#1515CB')
    
    plot(g6, add=TRUE, col='#15CBCB')
    text(0.68,0.05,'__', cex = 0.9,col = '#15CBCB')
    text(0.18,0.05,auc6, cex = 0.9,col='#15CBCB')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
      
    
    
 
    plot(g1, col = 'black')
    text(0.61,0.55,'__', cex = 0.9, col = 'black')
    text(0.18,0.55,auc1, cex = 0.9, col = 'black')
     
    plot(g2, add=TRUE, col='#9D9D9D')
    text(0.61,0.45,'__', cex = 0.9,col = '#9D9D9D')
    text(0.18,0.45,auc2, cex = 0.9,col='#9D9D9D')
    
    plot(g3, add=TRUE, col='#CB1515')
    text(0.61,0.35,'__', cex = 0.9,col = '#CB1515')
    text(0.18,0.35,auc3, cex = 0.9,col='#CB1515')
    
    plot(g4, add=TRUE, col='#F38026')
    text(0.61,0.25,'__', cex = 0.9,col = '#F38026')
    text(0.18,0.25,auc4, cex = 0.9,col='#F38026')
     
    plot(g5, add=TRUE, col='#1515CB')
    text(0.61,0.15,'__', cex = 0.9,col = '#1515CB')
    text(0.18,0.15,auc5, cex = 0.9,col='#1515CB')
    
    plot(g6, add=TRUE, col='#15CBCB')
    text(0.61,0.05,'__', cex = 0.9,col = '#15CBCB')
    text(0.18,0.05,auc6, cex = 0.9,col='#15CBCB')
  
    
    
    
    
    plot(g1, col = 'black')
    text(0.68,0.55,'__', cex = 0.9, col = 'black')
    text(0.37,0.55,auc11, cex = 0.9, col = 'black')
    text(0.37,0.5,auc12, cex = 0.9, col = 'black')
    
    plot(g2, add=TRUE, col='#9D9D9D')
    text(0.68,0.45,'__', cex = 0.9,col = '#9D9D9D')
    text(0.31,0.45,auc21, cex = 0.9,col='#9D9D9D')
    text(0.37,0.4,auc22, cex = 0.9,col='#9D9D9D')
    
    plot(g3, add=TRUE, col='#CB1515')
    text(0.68,0.35,'__', cex = 0.9,col = '#CB1515')
    text(0.24,0.35,auc31, cex = 0.9,col='#CB1515')
    text(0.37,0.3,auc32, cex = 0.9,col='#CB1515')
    
    plot(g4, add=TRUE, col='#F38026')
    text(0.68,0.25,'__', cex = 0.9,col = '#F38026')
    text(0.18,0.25,auc41, cex = 0.9,col='#F38026')
    text(0.37,0.2,auc42, cex = 0.9,col='#F38026')
    
    plot(g5, add=TRUE, col='#1515CB')
    text(0.68,0.15,'__', cex = 0.9,col = '#1515CB')
    text(0.27,0.15,auc51, cex = 0.9,col='#1515CB')
    text(0.37,0.1,auc52, cex = 0.9,col='#1515CB')
    
    plot(g6, add=TRUE, col='#15CBCB')
    text(0.68,0.05,'__', cex = 0.9,col = '#15CBCB')
    text(0.21,0.05,auc61, cex = 0.9,col='#15CBCB')
    text(0.37,0.0,auc62, cex = 0.9,col='#15CBCB')
    
    roc.test(g1, g2)
    roc.test(g1, g3)
    roc.test(g1, g4)
    roc.test(g1, g5)
    roc.test(g1, g6)
    roc.test(g2, g3)
    roc.test(g2, g4)
    roc.test(g2, g5)
    roc.test(g2, g6)
    roc.test(g3, g4)
    roc.test(g3, g5)
    roc.test(g3, g6)
    roc.test(g4, g5)
    roc.test(g4, g6)
    roc.test(g5, g6)
    
    Indexed.ADC.Median + indexVol2total        
    
    
    
    fit2 <- glm(status ~ PSA_grouped + MRI_EPE + MRI_SVI + EPEscore_grouped +  SbxISUP_grouped, data=df, family = binomial)
    df$prob2=predict(fit2)
    g2 <- roc(status~prob2, data=df, ci=TRUE)
    auc2 = paste("Pre-operative: Systematic biopsy AUC = ",toString(round(g2$auc, digits = 3)))
    
    fit3 <- glm(status ~ PSA_grouped + MRI_EPE + MRI_SVI + EPEscore_grouped + FbxISUP_grouped, data=df, family = binomial)
    df$prob3=predict(fit3)
    g3 <- roc(status~prob3, data=data2, ci=TRUE)
    auc3 = paste("Pre-operative: MRI-guided biopsy AUC = ",toString(round(g3$auc, digits = 3)))
    
    fit4 <- glm(status ~ PSA_grouped + MRI_EPE + MRI_SVI + EPEscore_grouped + FbxISUP_grouped + SbxISUP_grouped, data=df, family = binomial)
    df$prob4=predict(fit4)
    g4 <- roc(status~prob4, data=df, ci=TRUE)
    auc4 = paste("Pre-operative: Combined biopsy AUC = ",toString(round(g4$auc, digits = 3)))
    
    fit5 <- glm(status ~ PSA_grouped + MRI_EPE + MRI_SVI + EPEscore_grouped + FbxISUP_grouped + SbxISUP_grouped + RpISUP_grouped + path_EPE + path_SVI + margins, data=df, family = binomial)
    df$prob5=predict(fit5)
    g5 <- roc(status~prob5, data=df, ci=TRUE)
    auc5 = paste("Combination pre- and post-operative AUC = ",toString(round(g5$auc, digits = 3)))
    
    fit6 <- glm(status ~ PSA_grouped +  SbxISUP_grouped, data=df, family = binomial)
    df$prob6=predict(fit6)
    g6 <- roc(status~prob6, data=df, ci=TRUE)
    auc6 = paste("Pre-operative: PSA + systematic biopsy AUC = ",toString(round(g6$auc, digits = 3)))
    
    fit7 <- glm(status ~ PSA_grouped + FbxISUP_grouped, data=df, family = binomial)
    df$prob7=predict(fit7)
    g7 <- roc(status~prob7, data=df, ci=TRUE)
    auc7 = paste("Pre-operative: PSA + systematic biopsy AUC = ",toString(round(g7$auc, digits = 3)))
    
    fit8 <- glm(status ~ PSA_grouped + FbxISUP_grouped + SbxISUP_grouped, data=df, family = binomial)
    df$prob8 = predict(fit8)
    g8 <- roc(status~prob8, data=df, ci=TRUE)
    auc8 = paste("Pre-operative: PSA + MRI-guided biopsy AUC = ",toString(round(g8$auc, digits = 3)))
    
    
    
    
    
    plot(g1, col = 'black')
    text(0.82,0.15,'--', col = 'black')
    text(0.46,0.15,auc1, col = 'black')
    
    plot(g2, add=TRUE, col='red2')
    text(0.82,0.1,'--', col = 'red2')
    text(0.33,0.1,auc2, col='red2')
    
    plot(g5, add=TRUE, col='#2166AC')
    text(0.82,0.05,'--', col = '#2166AC')
    text(0.39,0.05,auc5, col='#2166AC')
    
    
    
    
    
    
    
    
    fit9 <- glm(status ~ PSA + FbxISUP + MRI_EPE + MRI_SVI, data=df, family = binomial)
    prob9=predict(fit9)
    df$prob9 = prob9
    g9 <- roc(status~prob9, data=df, ci=TRUE)
    auc9 = paste("surgical margins only = ",toString(round(g9$auc, digits = 3)))
    
    # single feature
    plot(g9, col='#00FF33')
    text(0.75,0.25,'--', col='#00FF33')
    text(0.45,0.25,auc9, col='#00FF33')
    
    
    
    palette_check(c("#2166AC", "#92C5DE","#A50F15","#EF3B2C"),plot=TRUE)
    display.brewer.pal(n = 8, name = 'RdBu')
    brewer.pal(n=8, name='RdBu')
    
    #"purple", "#2166AC", "orange","#FF3300","#B2182B"
    
    # 1  Postoperative        
    plot(g1, col = 'black')
    text(0.82,0.20,'--', col = 'black')
    text(0.55,0.20,auc1, col = 'black')
    
    #2 PSA + staging + Sbx
    plot(g2, add=TRUE, col='#B2182B')
    text(0.82,0.15,'--', col='#B2182B')
    text(0.425,0.15,auc2, col='#B2182B')
    
    #3 PSA + staging + Fbx
    plot(g3, add=TRUE, col='#D6604D')
    text(0.82,0.10,'--', col='#D6604D')
    text(0.425,0.10,auc3, col='#D6604D')
    
    #4 PSA + staging + combo
    plot(g4, add=TRUE, col='#92C5DE')
    text(0.82,0.05,'--', col='#92C5DE')
    text(0.432,0.05,auc4, col='#92C5DE')
    
    #5 All features combined 
    plot(g5, add=TRUE, col='#2166AC')
    text(0.82,0.0,'--', col = '#2166AC')
    text(0.405,0.0,auc5, col='#2166AC')        
    
    #6 PSA and Sbx
    plot(g6, add = TRUE, col = '#2166AC')
    text(0.75,0.15,'--', col = '#2166AC')
    text(0.292,0.15,auc6, col = '#2166AC')
    
    #7 PSA and Fbx
    plot(g7, add=TRUE, col='black')
    text(0.75,0.2,'--', col='black')
    text(0.285,0.2,auc7, col='black')
    
    #8 PSA and Combo
    plot(g8, add=TRUE, col='black')
    text(0.75,0.2,'--', col='black')
    text(0.285,0.2,auc8, col='black')
    
    # 1  Postoperative
    # 2  PSA + staging + Sbx
    # 3  PSA + staging + Fbx
    # 4  PSA + staging + combo
    # 5  Combination Pre- and post-operative
    # 6  PSA and Sbx
    # 7  PSA and Fbx
    # 8  PSA and combo
    # 9  single feature
    auc(g1)
    roc.test(g1, g2)
    
    
    
    
    summary(fit1)
    g1$auc
    g1$ci
    summary(fit2)
    g2$auc
    g2$ci
    roc.test(g1,g2)
    summary(fit3)
    g3$auc
    g3$ci
    roc.test(g1,g3)
    summary(fit4)
    g4$auc
    g4$ci
    roc.test(g1,g4)
    summary(fit5)
    g5$auc
    g5$ci
    roc.test(g1,g5)
    summary(fit6)
    g6$auc
    g6$ci
    roc.test(g1,g2)
    summary(fit7)
    g7$auc
    g7$ci
    roc.test(g1,g2)
    summary(fit8)
    g8$auc
    g8$ci
    
    
    
    
    
    
    
    
    
 