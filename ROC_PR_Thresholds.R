#https://cran.r-project.org/web/packages/PRROC/vignettes/PRROC.pdf
# https://cran.r-project.org/web/packages/cutpointr/vignettes/cutpointr.html


library(precrec)
library(PRROC)


# Read Data File
f <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_neg.csv"
data_neg <- read.csv(f, header=TRUE, stringsAsFactors=FALSE)

g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_pos.csv"
data_pos <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

h <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData.csv"
data_full <- read.csv(h, header=TRUE, stringsAsFactors=FALSE)


fg1 = data_pos$SVR
fg2 = data_pos$invSVR
fg3 = data_pos$minusSVR
bg1 = data_neg$SVR
bg2 = data_neg$invSVR
bg3 = data_neg$minusSVR

roc1<-roc.curve(scores.class0=fg1, scores.class1=bg1, curve=TRUE)
pr1<-pr.curve(scores.class0=fg1, scores.class1=bg1, curve=TRUE)
roc2<-roc.curve(scores.class0=fg2, scores.class1=bg2, curve=TRUE)
pr2<-pr.curve(scores.class0=fg2, scores.class1=bg2, curve=TRUE)
roc3<-roc.curve(scores.class0=fg3, scores.class1=bg3, curve=TRUE)
pr3<-pr.curve(scores.class0=fg3, scores.class1=bg3, curve=TRUE)

roc1
pr1

plot(roc1)
plot(roc2)
plot(roc3)

curve.points1<-roc$curve

library(cutpointr)

opt_cut <- cutpointr(data_full, invSVR, Status, direction = ">=", pos_class = 1,
                     neg_class = 0, method = maximize_metric, metric = youden)

points(roc$curve[253,3],pch="+",col="black")



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


