library(precrec)

# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Project1_Rstudio_ADCvolume.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

# Logistic Regression Models
library(nnet)
#model 1
mymodel1 <- multinom(Status~SVR, data = data1)
title1 = "ROC Curve for Minimum ADC Mean"
#model 2
mymodel2 <- multinom(Status~PIRADS4, data = data1)
title2 = "ROC Curve for Minimum ADC Med"
#model 3
mymodel3 <- multinom(Status~PIRADS3, data = data1)
title3 = "ROC Curve for Minimum ADC 10 Percent"
#model4
mymodel4 <- multinom(Status~TotalVol, data = data1)
title4 = "ROC Curve for Maximum ADC Entropy"
#model5
mymodel5 <- multinom(Status~PSAdensity, data = data1)
title5 = "ROC Curve for Minimum ADC Entropy"



# Model Performance Evaluation
library(ROCR)
pred1 <- predict(mymodel1, data1, type = 'prob')
pred2 <- predict(mymodel2, data1, type = 'prob')
pred3 <- predict(mymodel3, data1, type = 'prob')
pred4 <- predict(mymodel4, data1, type = 'prob')
pred5 <- predict(mymodel5, data1, type = 'prob')

# Calculate ROC and Precision-Recall curves
sscurves2 <- evalmod(scores = pred1, labels = data1$Status)

# Show ROC and Precision-Recall plots
plot(sscurves2)

# Show a Precision-Recall plot
plot(sscurves2, "PRC")


# ### Plot ROC and PR curves
# The ggplot2 package is required 
library(ggplot2)

# Show both ROC and Precision-Recall plots
autoplot(sscurves2)

# Show a Precision-Recall plot only
autoplot(sscurves2, "PRC")

# Show AUC scores in a table
aucs <- auc(sscurves2)
knitr::kable(aucs)



#The join_scores function combines multiple score datasets.
scores1 <- join_scores(pred1, pred2, pred3, pred4, pred5)
# Join two label vectors
labels1 <- join_labels(data1$Status, data1$Status, data1$Status, data1$Status, data1$Status)

# Specify model names and dataset IDs
mmmdat <- mmdata(scores1, labels1, modnames = c("PIRADS5", "PIRADS4", "PIRADS3", "TotalVol", "PSAdensity"), dsids = c(1, 2, 3, 4, 5))
 

# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(mmmdat)

# Show ROC and Precision-Recall curves with the ggplot2 package
autoplot(mscurves)


# Show average Precision-Recall curves with the 95% confidence bounds
autoplot(mscurves, "PRC", show_cb = TRUE)


# Get a data frame with AUC scores
aucs2 <- auc(mscurves)

# Use knitr::kable to display the result in a table format
knitr::kable(aucs2)

#Get AUCs of Precision-Recall
aucs_roc2 <- subset(aucs2, curvetypes == "ROC")
knitr::kable(aucs_roc2)

#Get AUCs of Precision-Recall
aucs_prc2 <- subset(aucs2, curvetypes == "PRC")
knitr::kable(aucs_prc2)
