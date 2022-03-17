library(glmnet)
library(survival)

g1 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_x.csv"
data_x <- read.csv(g1, header=TRUE, stringsAsFactors=FALSE)
x1<-data.matrix(data_x)

g2 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_y.csv"
data_y <- read.csv(g2, header=TRUE, stringsAsFactors=FALSE)
y1<-data.matrix(data_y)

fit1 <- glmnet(x1,y1,family="cox")
plot(fit1)
coef(fit1, s=0.04)

#cross-validation
set.seed(1)
cvfit1 <- cv.glmnet(x1,y1,family="cox", type.measure="C")
plot(cvfit1)

#find optimal lambda (where CV-error curve hits its minimum)
lam<-cvfit1$lambda.min
lam
coef(fit1, s=lam)

#find lambda for regularized model with CV-error within 1 st.dev of min
lam_se<-cvfit1$lambda.1se
lam_se
coef(fit1, s=lam_se)

