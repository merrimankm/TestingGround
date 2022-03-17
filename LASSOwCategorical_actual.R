library(glmnet)
library(survival)

g1 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_x.csv"
data_x <- read.csv(g1, header=TRUE, stringsAsFactors=FALSE)
x1<-data.matrix(data_x)

g2 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData_y.csv"
data_y <- read.csv(g2, header=TRUE, stringsAsFactors=FALSE)
y1<-data.matrix(data_y)




#####  COX Analysis ##############

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




####  BINOMIAL ANALYSIS - all data #####
# Create a y vector containing ONLY status data
ybin = y1[1:555,2]

# fit binomial logistic regression
fitbin <-glmnet(x1,ybin, family = "binomial")
plot(fitbin)

# cross-validation code using misclassification error as measurement criterion
cvfitbin <- cv.glmnet(x1, ybin, family = "binomial", type.measure = "class")
plot(cvfitbin)


#find optimal lambda (where CV-error curve hits its minimum)
lam_bin<-cvfitbin$lambda.min
lam_bin
coef(fitbin, s=lam_bin)

#find lambda for regularized model with CV-error within 1 st.dev of min
lam_bin_se<-cvfitbin$lambda.1se
lam_bin_se
coef(fitbin, s=lam_bin_se)






####  BINOMIAL ANALYSIS - 3 year #####

g3 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\x_3yr.csv"
data_x3 <- read.csv(g3, header=TRUE, stringsAsFactors=FALSE)
x3<-data.matrix(data_x3)

g4 <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\Y_3yr.csv"
data_y3 <- read.csv(g4, header=TRUE, stringsAsFactors=FALSE)
y3<-data.matrix(data_y3)


# fit binomial logistic regression
fit3 <-glmnet(x3,y3, family = "binomial")
plot(fit3)

# cross-validation code using misclassification error as measurement criterion
cvfit3 <- cv.glmnet(x3, y3, family = "binomial", type.measure = "class")
plot(cvfit3)


#find optimal lambda (where CV-error curve hits its minimum)
lam_bin3<-cvfit3$lambda.min
lam_bin3
coef(fit3, s=lam_bin3)

#find lambda for regularized model with CV-error within 1 st.dev of min
lam_bin_se3<-cvfit3$lambda.1se
lam_bin_se3
coef(fit3, s=lam_bin_se3)


