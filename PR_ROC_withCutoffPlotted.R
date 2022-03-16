#https://cran.r-project.org/web/packages/precrec/vignettes/introduction.html

library(precrec)

# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\ContinuousData.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

# Logistic Regression Models
library(nnet)
#model 1
mymodel1 <- multinom(Status~Vol_P5, data = data1)
title1 = "Relative Total PIRADS 5 Volume"
#model 2
mymodel2 <- multinom(Status~PSA_density, data = data1)
title2 = "PSA Density"
#model 3
mymodel3 <- multinom(Status~ADCent_wp, data = data1)
title3 = "Whole Prostate ADC Entropy"
#model4
mymodel4 <- multinom(Status~Vol_indexed, data = data1)
title4 = "Indexed Lesion Relative Volume"
#model5
mymodel5 <- multinom(Status~SVR, data = data1)
title5 = "Indexed Lesion Surface to Volume Ratio"
#model6
mymodel6 <- multinom(Status~X3Ddiameter, data = data1)
title6 = "Indexed Lesion Maximum 3D Diamete"


# Model Performance Evaluation
library(ROCR)
pred1 <- predict(mymodel1, data1, type = 'prob')
pred2 <- predict(mymodel2, data1, type = 'prob')
pred3 <- predict(mymodel3, data1, type = 'prob')
pred4 <- predict(mymodel4, data1, type = 'prob')
pred5 <- predict(mymodel5, data1, type = 'prob')
pred6 <- predict(mymodel6, data1, type = 'prob')

# Calculate ROC and Precision-Recall curves
sscurves1 <- evalmod(scores = pred1, labels = data1$Status)
sscurves2 <- evalmod(scores = pred2, labels = data1$Status)
sscurves3 <- evalmod(scores = pred3, labels = data1$Status)
sscurves4 <- evalmod(scores = pred4, labels = data1$Status)
sscurves5 <- evalmod(scores = pred5, labels = data1$Status)
sscurves6 <- evalmod(scores = pred6, labels = data1$Status)

# Show ROC and Precision-Recall plots
plot(sscurves1)
plot(sscurves2)
plot(sscurves3)
plot(sscurves4)
plot(sscurves5)
plot(sscurves6)

# Show AUC scores in a table
aucs1 <- auc(sscurves1)
knitr::kable(aucs1, modnames = "test")
aucs2 <- auc(sscurves2)
knitr::kable(aucs2)
aucs3 <- auc(sscurves3)
knitr::kable(aucs3)
aucs4 <- auc(sscurves4)
knitr::kable(aucs4)
aucs5 <- auc(sscurves5)
knitr::kable(aucs5)
aucs6 <- auc(sscurves6)
knitr::kable(aucs6)



#The join_scores function combines multiple score datasets.
#scores1 <- join_scores(pred1, pred2, pred3)
#scores1 <- join_scores(pred1, pred2, pred3, pred4, pred5)
scores1 <- join_scores(pred1, pred2, pred3, pred4, pred5, pred6)
# Join two label vectors
#labels1 <- join_labels(data1$Status, data1$Status, data1$Status)
#labels1 <- join_labels(data1$Status, data1$Status, data1$Status, data1$Status, data1$Status)
labels1 <- join_labels(data1$Status, data1$Status, data1$Status, data1$Status, data1$Status, data1$Status)

# Specify model names and dataset IDs
#mmmdat <- mmdata(scores1, labels1, modnames = c(title1, title2, title3), dsids = c(1, 2, 3))
#mmmdat <- mmdata(scores1, labels1, modnames = c(title1, title2, title3, title4, title5), dsids = c(1, 2, 3, 4, 5))
mmmdat <- mmdata(scores1, labels1, modnames = c(title1, title2, title3, title4, title5, title6), dsids = c(1, 2, 3, 4, 5, 6))

# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(mmmdat)

# Show ROC and Precision-Recall curves with the ggplot2 package
autoplot(mscurves)


# Show average Precision-Recall curves with the 95% confidence bounds
autoplot(mscurves, "PRC", show_cb = TRUE)


# Get a data frame with AUC scores
multiaucs <- auc(mscurves)

# Use knitr::kable to display the result in a table format
knitr::kable(multiaucs)

#Get AUCs of Precision-Recall
aucs_roc2 <- subset(multiaucs, curvetypes == "ROC")
knitr::kable(aucs_roc2)

#Get AUCs of Precision-Recall
aucs_prc2 <- subset(multiaucs, curvetypes == "PRC")
knitr::kable(aucs_prc2)



# Show plots individually along with threshold point (calculated using PR_Thresholds.R)

plot(sscurves2, "ROC")
x1 = c(0.42391)
y1 = c(0.719512)
points(x1,y1,pch=19,cex = 1.5, col="red")
text(x1,y1, labels="Threshold = 0.0263",pos = 2, cex=0.9, font=2)

x2=c(0.7317)
y2=c(0.2353)
plot(sscurves2, "PRC")
points(x2,y2,pch=19,cex = 1.5, col="red")
text(x2,y2, labels="Threshold = 0.0263",pos = 3, cex=0.9, font=2)


