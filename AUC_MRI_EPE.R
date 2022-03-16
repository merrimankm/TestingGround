
# https://www.youtube.com/watch?v=ypO1DPEKYFo


# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\EPE_SVI.csv"


data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

#view structure of data
str(data1)

# Logistic Regression Model
library(nnet)
mymodel <- multinom(Status~MRI_EPE, data = data1)

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
pred <- predict(mymodel, data1, type = 'prob')
head(pred)
head(data1)

hist(pred)
pred1 <- prediction(pred, data1$Status)
eval <- performance(pred1, "acc")
plot(eval)
abline(h=0.71, v=0.45)


# Identify Best Values
max <- which.max(slot(eval, "y.values")[[1]])
max
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
print(c(Accuracy = acc, Cutoff = cut))

# ROC
roc <- performance(pred1, "tpr", "fpr")
plot(roc,
     colorize=T,
     main = "ROC Curve",
     ylab = "Sensitivity",
     xlab = "1-Specificity")
# color in plot above will correspond to cutoff
abline(a=0, b=1)

# Area Under Curve (AUC)
auc <- performance(pred1, "auc")
auc <- unlist(slot(auc, "y.values"))
auc <- round(auc, 4)
legend(.6, .2, auc, title = "AUC", cex = 0.8)


# video on handling imbalanced data:
# https://www.youtube.com/watch?v=Ho2Klvzjegg
