library(glmnet)
library(survival)


#### EXAMPLE  ####
data(CoxExample)
x<-CoxExample$x
y<-CoxExample$y

fit <- glmnet(x,y,family="cox")
plot(fit)
coef(fit, s=0.05)

set.seed(1)
cvfit <- cv.glmnet(x,y,family="cox", type.measure="C")
plot(cvfit)

#### ACTUAL DATA #####

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



####### Example for stop-start data #####
#I don't think this is actually needed as my data isn't left censored?
set.seed(2)
nobs <- 100; nvars <- 15
xvec <- rnorm(nobs * nvars)
xvec[sample.int(nobs * nvars, size = 0.4 * nobs * nvars)] <- 0
x <- matrix(xvec, nrow = nobs) # non-sparse x
x_sparse <- Matrix::Matrix(xvec, nrow = nobs, sparse = TRUE) # sparse x
# create start-stop data response
beta <- rnorm(5)
fx <- x_sparse[, 1:5] %*% beta / 3
ty <- rexp(nobs, drop(exp(fx)))
tcens <- rbinom(n = nobs, prob = 0.3, size = 1)
starty <- runif(nobs)
yss <- Surv(starty, starty + ty, tcens)