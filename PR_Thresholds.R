#https://cran.r-project.org/web/packages/PRROC/vignettes/PRROC.pdf
# https://cran.r-project.org/web/packages/cutpointr/vignettes/cutpointr.html


# write prediction scores into column for easy addition to excel doc
for (i in 1:length(pred6)) {
  cox.coef[i-2,8] <- summary(cox.univ[[i]])$rsq[1]
}
dft<-as.data.frame(t(pred6))
write.table(df, file = "clipboard", sep="\t", row.names = FALSE, col.names = FALSE) #copy to clipboard



library(precrec)
library(PRROC)


# Read Data File
f <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllData_outlierRemoved_neg.csv"
data_neg <- read.csv(f, header=TRUE, stringsAsFactors=FALSE)

g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllData_outlierRemoved_pos.csv"
data_pos <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

h <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllData_outlierRemoved.csv"
data_full <- read.csv(h, header=TRUE, stringsAsFactors=FALSE)


feat = "Vol_P4"
fg1 = data_pos$Vol_P4
bg1 = data_neg$Vol_P4

roc<-roc.curve(scores.class0=fg1, scores.class1=bg1, curve=TRUE)
pr<-pr.curve(scores.class0=fg1, scores.class1=bg1, curve=TRUE)

roc
pr

plot(roc)

curve.points1<-roc$curve

library(cutpointr)

opt_cut <- cutpointr(data_full, Vol_P4, Status, direction = ">=", pos_class = 1,
                     neg_class = 0, method = maximize_metric, metric = youden)

points(roc$curve[253],pch="+",col="black")
points(0.424,pch="+",col="black")
points(0.2,pch="+",col="purple")
points(0.3,pch="+",col="green")
points(0.45,pch="o",col="blue")
annotate("point", x = 0.42391, y = 0.719512, colour = "red")
annotate("point", x = 0.6, y = 0.4, colour = "purple",pch="+")








text(dist ~speed, labels=rownames(cars),data=cars, cex=0.9, font=2)

plot(pr)

curve.points2<-pr$curve

f1 <- vector()
#f1 <- array(-1, 1,dim=c(dim(curve.points2)[1]))
for (i in 1:dim(curve.points2)[1]) {
  f1[i]=2/(1/curve.points2[i,1]+1/curve.points2[i,2])
}

feat
roc$auc
pr$auc.integral
opt_cut$youden
opt_cut$optimal_cutpoint
max(f1)
which.max(f1)
curve.points2[which.max(f1),1:3]
