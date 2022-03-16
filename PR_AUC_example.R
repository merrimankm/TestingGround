library(precrec)

# Load a test dataset
data(P10N10)

# Calculate ROC and Precision-Recall curves
sscurves <- evalmod(scores = P10N10$scores, labels = P10N10$labels)

# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Project1_Rstudio_ADCvolume.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

#view structure of data
str(data1)

# Logistic Regression Model
library(nnet)
mymodel <- multinom(Status~ADC10, data = data1)
title = "ROC Curve for Minimum ADC 10 Percent"

# Misclassification Rate
p <- predict(mymodel, data1)
tab <- table(p, data1$Status)
tab
#    compare to actual percentages 
#    (is model accuracy better than simply stating overall odds of BCR?)
sum(diag(tab))/sum(tab)
1-sum(diag(tab))/sum(tab)
table(data1$Status)

# Model Performance Evaluation
library(ROCR)
pred1 <- predict(mymodel, data1, type = 'prob')

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
autoplot(sscurves)

# Show a Precision-Recall plot only
autoplot(sscurves, "PRC")

# Show AUC scores in a table
aucs <- auc(sscurves)
knitr::kable(aucs)






# ##### Combine multiple data
# Create second Logistic Regression Model
mymodel2 <- multinom(Status~ADCmed, data = data1)
title2 = "ROC Curve for Minimum ADC Median"

# Misclassification Rate
p2 <- predict(mymodel2, data1)
pred2 <- predict(mymodel2, data1, type = 'prob')

#The join_scores function combines multiple score datasets.
scores1 <- join_scores(pred1, pred2)
# Join two label vectors
labels1 <- join_labels(data1$Status, data1$Status)

# Specify model names and dataset IDs
mmmdat <- mmdata(scores1, labels1, modnames = c("ADC10", "ADCmed"), dsids = c(1, 2))
 

# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(mmmdat)

# Show ROC and Precision-Recall curves with the ggplot2 package
autoplot(mscurves)


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