# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

g1 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_neg.csv"
data_neg <- read.csv(g1, header=TRUE, stringsAsFactors=FALSE)

g2 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_pos.csv"
data_pos <- read.csv(g2, header=TRUE, stringsAsFactors=FALSE)


library(foreign)
library(nnet)
library(ggplot2)
library(reshape2)
library(precrec)
library(ROCR)
library(survival)
library(PRROC)
library(cutpointr)



opt.cut = function(in.cut, in.sens, in.spec){
  d = in.sens + in.spec -1
  ind = which(d == max(d))
  ind = ind[1]
  opt.out <- c(cutoff = in.cut[[ind]])
  return(opt.out)
}




####### VERSION 1 - LOGISTIC REGRESSION ########

# create training and validation sets
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(data1$Status, p = .7,
                                  list = FALSE,
                                  times = 1)
Train <- data1[ trainIndex,]
Valid <- data1[-trainIndex,]




###!!!!!UPDATE
test <- glm(Status ~ ADCmin_mean, family = binomial(), data = Train)

summary(test)

# perform 2 tailed p test
z <- summary(test)$coefficients/summary(test)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

# relative risk
exp(coef(test))


# predictions using multinom with "predict"
pred <- predict(test, Valid, type = 'response')
sscurves <- evalmod(scores = pred, labels = Valid$Status)
plot(sscurves)
auc <- precrec::auc(sscurves)
knitr::kable(auc)



######### VERSION 2 ###############

# predictions using variable with "prediction"



###!!!! UPDATE
pred2 <- prediction(data1$ADCmin_mean, data1$Status)

#prediction outcome from ADC data
perf2 <- performance(pred2, "auc") #calculate AUC from Variable1
est.auc = perf2@y.values[[1]]

# get acc, tpr, fpr plot information so we can see how performance changes at each Variable1 value
perf1 <- performance(pred2, "tpr", "fpr")
perf3 <- performance(pred2, "acc")
# sensitivity and specificity values
sensitivity = perf1@y.values[[1]]
specificity = 1-perf1@x.values[[1]]
# cutoffs (i.e. Variable1 values)
cutoffs = perf3@x.values[[1]]
# calculate optimal cut point based on our function
est.cut = opt.cut(cutoffs, sensitivity, specificity)

est.auc
est.cut





########### VERSION 3 ##############

###!!!! UPDATE
fg1 = data_pos$Vol_P5
bg1 = data_neg$Vol_P5

roc<-roc.curve(scores.class0=fg1, scores.class1=bg1, curve=TRUE)
pr<-pr.curve(scores.class0=fg1, scores.class1=bg1, curve=TRUE)

roc
pr

plot(roc)
curve.points1<-roc$curve




###!!!! UPDATE
opt_cut <- cutpointr(data1,Vol_P5, Status, direction = ">=", pos_class = 1,
                     neg_class = 0, method = maximize_metric, metric = youden)
opt_cut$optimal_cutpoint



plot(pr)
curve.points2<-pr$curve
f1 <- vector()
#f1 <- array(-1, 1,dim=c(dim(curve.points2)[1]))
for (i in 1:dim(curve.points2)[1]) {
  f1[i]=2/(1/curve.points2[i,1]+1/curve.points2[i,2])
}
max(f1)
which.max(f1)
pr$curve[which.max(f1),1:3]  


