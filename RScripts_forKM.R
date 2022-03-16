#####################################################
######### EXAMPLE DATAFRAME SET UP ##################
#####################################################
# for the purpose of the below code lets say I have a dataframe called "inData" with the following variables:
# Time = time to death/censor
# Censor = censor=0 if data is censored, = 1 if observed
# Variable1 = continous imaging variable
# Variable2 = continous imaging variable
# Variable3 = continous imaging variable
# Variable4 = categorical imaging variable


##########################################################################
############## EXAMPLE AUC AND YOUDEN INDEX CALCULATIONS #################
##########################################################################
  library(ROCR)
  
  #### write a function to calculate optimal cut-point based on Youden Index
  opt.cut = function(in.cut, in.sens, in.spec){
    d = in.sens + in.spec -1
    ind = which(d == max(d))
    ind = ind[1]
    opt.out <- c(cutoff = in.cut[[ind]])
    return(opt.out)
  }
  
  
  # remember that censor = 0 is not observed (no recurrence) and censor=1 is obersved (recurrence)
  pred <- prediction(inData$Variable1, inData$Censor) #prediction outcome from ADC data
  perf2 <- performance(pred, "auc") #calculate AUC from Variable1
  est.auc[i,1] = perf2@y.values[[1]]
  
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



######################################################
####### CATEGORICAL KM and COX-PH ANALYSIS ###########
######################################################
  # in my silly example dataframe, Variable4 is the only categorical variable
  # variable 4 kaplan meier
  surv.object <- Surv(inData$Time, inData$Outcome)
  myfits.v4 <- survfit(surv.object~Variable4, data = inData)
  myfits.v4.LR = survdiff(surv.object~Variable4,data=inData)
  myfits.v4.p = 1-pchisq(myfits.v4.LR$chisq, 1)
  plot(myfits.v4, col=c(1,2,4), lwd=2,main="This plots Variable 4 KM", xlab="Time", ylab="Progression Free Survival", cex.lab=2,cex.axis=2, mark.time = TRUE, mark =3, cex=2, cex.main = 2)
  # to add a legend we need to know the category names
  # in our fictious example it has 2 groups
  legend("topright", c(paste("v4 low, N=", myfits.v4$n[1], sep=""),paste("v4 high, N=", myfits.v4$n[2], sep="")),lwd=2, col=c(1,2,4), cex=1.75)
  #lets put p-value on the plot
  text(5,0.05,paste("P =",signif(as.numeric(myfits.v4.p),2)), cex=1.2)
  
  #variable 4 cox proportional-hazard
  v4_cph = coxph(surv.object ~ Variable4,data=inData)
  summary(v4_cph)



######################################################
######## CONTINUOUS VARIABLES COX-PH ANALYSIS ########
######################################################
  # as a short-cut, here is a neat way to print all continous variable data into a table
  
  # example for selecting variables by name
  # I only select continuous variables for this
  x_cox <- inData[, names(inData) %in% c("Time","Outcome","Variable1","Variable2","Variable3")]
  
  surv.object <- Surv(x_cox$Time, x_cox$Outcome)
  
  #do univariate regression on each variable individually
  cox.univ <- vector(mode="list", length=(dim(x_cox)[2]-1))
  cox.coef <- array(-1, dim=c(dim(x_cox)[2],8))
  colnames(cox.coef) <- c("name", "p_val", "coef", "exp(coef)","SE_coef","low_95CI","up_95CI", "rsq");
  
  # since my first two variables are Time and Outcome, I start at 3
  for (i in 3:dim(x_cox)[2]) {
    # I scale the data so that the HR are comparable
    cox.univ[[i]] <- coxph(surv.object ~ scale(x_cox[,i]))
    cox.coef[i,1] <- colnames(x_cox)[i] #variable name
    cox.coef[i,2] <- summary(cox.univ[[i]])$coefficients[1,5] #p-val
    cox.coef[i,3] <- summary(cox.univ[[i]])$coefficients[1,1] #the coefficient
    cox.coef[i,4] <- summary(cox.univ[[i]])$coefficients[1,2] #the exponential of the coefficient
    cox.coef[i,5] <- summary(cox.univ[[i]])$coefficients[1,3] #SE of coefficient
    cox.coef[i,6] <- summary(cox.univ[[i]])$conf.int[1,3] #lower 95% CI
    cox.coef[i,7] <- summary(cox.univ[[i]])$conf.int[1,4] #upper 95% CI
    cox.coef[i,8] <- summary(cox.univ[[i]])$rsq[1]
  }
  
  write.table(cox.coef, file = "clipboard", sep="\t", row.names = FALSE, col.names = FALSE) #copy to clipboard


#########################################################
########## NUMERICAL FIND OPTIMAL CUTOFF ################
#########################################################
  #finding optimal cutoff for imaging data
  # you can either make a new dataframe with only significant variables or you can just run them all
  #example for selecting variables by name
  sigVars <- x_cox[, names(x_cox) %in% c("Time","Outcome","Variable1","Variable2")]
  surv.object <- Surv(sigVars$Time, sigVars$Outcome)
  # PLOTTING CATEGORICAL DATA, make sure this folder exists
  # the for loop will save plots for all variables
  directory2save = "C:/Users/harmonsa/Documents/KM_curves/"
  
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




######################################################
######## MULTIVARIABLE COX-PH ANALYSIS ###############
######################################################
  library("survival", "KMsurv", "Hmisc")
  
  #after you find how each variable does in univariate assessment (both continuous and categorical), now we look at how they all do together
  # "forward" selection here is based on our univariate assessment, you can either pick (1) all variables, (2) variables with p<0.1 or p<0.2 in univariate, or (3) known variables
  multiX <- inData[, names(inData) %in% c("Time","Outcome","Variable1","Variable2","Variable 4")]
  # we only fit multivariable models in patients that have all variables observed, i.e. no "NA" values
  m_cox <- multiX[complete.cases(multiX), ]
  
  #fit the initial model
  cox.multi.1 <- coxph(Surv(m_cox$Time, m_cox$Outcome) ~  Variable1 + Variable2 + Variable4, data = m_cox)
  summary(cox.multi.1)
  
  #perform AIC selection
  cox.multi.AIC <- step(cox.multi.1, direction = "backward")
  
  #perform BIC selection
  n <- dim(m_cox)[1]
  cox.multi.BIC <- step(cox.multi.1, k=log(n), direction = "backward")
  
  # BIC tends to produce smaller models (fewer variables)
  summary(cox.multi.BIC)
  summary(cox.multi.AIC)    
  
  # please note this multivariable model kept the continuous variables as continuous, if you want to make them categorical based on optimal cut-point create a new variable



