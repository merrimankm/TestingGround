#https://cran.r-project.org/web/packages/precrec/vignettes/introduction.html

library(precrec)

# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)


# Logistic Regression Models
library(nnet)
#model 1
mymodel1 <- multinom(Status~ADCmin_10, data = data1)
title1 = "Relative Total PIRADS 5 Volume"
#model 2
mymodel2 <- multinom(Status~invSVR, data = data1)
title2 = "PSA Density"
#model 3
mymodel3 <- multinom(Status~minusSVR, data = data1)
title3 = "Indexed Lesion Relative Volume"
#model4
mymodel4 <- multinom(Status~Vol_P5+PSA_density+Vol_indexed+SVR+X3Ddiameter, data = data1)
title4 = "Top Performers Combined"




# Model Performance Evaluation
#library(ROCR)
#library(ggplot2)
pred1 <- predict(mymodel1, data1, type = 'prob')
pred2 <- predict(mymodel2, data1, type = 'prob')
pred3 <- predict(mymodel3, data1, type = 'prob')
pred4 <- predict(mymodel4, data1, type = 'prob')

# Calculate ROC and Precision-Recall curves
sscurves1 <- evalmod(scores = pred1, labels = data1$Status)
sscurves2 <- evalmod(scores = pred2, labels = data1$Status)
sscurves3 <- evalmod(scores = pred3, labels = data1$Status)
sscurves4 <- evalmod(scores = pred4, labels = data1$Status)

# Show individual ROC and Precision-Recall plots
plot(sscurves1)
plot(sscurves2)
plot(sscurves3)
plot(sscurves4)

# Show individual AUC scores in a table
aucs1 <- precrec::auc(sscurves1)
knitr::kable(aucs1)
aucs2 <- auc(sscurves2)
knitr::kable(aucs2)
aucs3 <- auc(sscurves3)
knitr::kable(aucs3)
aucs4 <- auc(sscurves4)
knitr::kable(aucs4)



#The join_scores function combines multiple score datasets.
scores1 <- join_scores(pred1, pred2, pred3)
#scores1 <- join_scores(pred1, pred2, pred3, pred4, pred5)
#scores1 <- join_scores(pred1, pred2, pred3, pred4, pred5, pred6)
# Join two label vectors
labels1 <- join_labels(data1$Status, data1$Status, data1$Status)
#labels1 <- join_labels(data1$Status, data1$Status, data1$Status, data1$Status, data1$Status)
#labels1 <- join_labels(data1$Status, data1$Status, data1$Status, data1$Status, data1$Status, data1$Status)

# Specify model names and dataset IDs
mmmdat <- mmdata(scores1, labels1, modnames = c(title1, title2, title3), dsids = c(1, 2, 3))
#mmmdat <- mmdata(scores1, labels1, modnames = c(title1, title2, title3, title4, title5), dsids = c(1, 2, 3, 4, 5))
#mmmdat <- mmdata(scores1, labels1, modnames = c(title1, title2, title3, title4, title5, title6), dsids = c(1, 2, 3, 4, 5, 6))

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




# write prediction scores into column for easy addition to excel doc
for (i in 1:length(pred6)) {
  cox.coef[i-2,8] <- summary(cox.univ[[i]])$rsq[1]
}
dft<-as.data.frame(t(pred6))
write.table(df, file = "clipboard", sep="\t", row.names = FALSE, col.names = FALSE) #copy to clipboard

