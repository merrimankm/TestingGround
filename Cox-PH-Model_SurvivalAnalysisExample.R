



library(survival)

## Get Stanford Heart Transplant Dataset
mismatchlevel<-cut(stanford2$t5, breaks=3, labels = c(0, 1, 2))
over40<-cut(stanford2$age, breaks=c(-Inf,40,Inf), labels=c(0,1))
stantemp<-cbind(stanford2, mismatchlevel, over40)
# remove any rows that have missing values (this part important for anova)
numbermissing<-apply(stantemp, 1, function(x) sum(is.na(x)))
nonemissing<-numbermissing==0
stan<-stantemp[nonemissing,]
rm(stantemp)
# if data wasn't already imported as a factor, would need to use
# stan$over40 <- as.factor(stan$over40)

summary(stan)

# look at KM model of data using only mismatch level
stan.KM <- survfit( Surv(stan$time, stan$status) ~ stan$mismatchlevel, type="kaplan-meier")
plot(stan.KM, conf.int=F, xlab="Time (days)", ylab="Survival = S(t)", main="KM-Model", col=c("red", "blue", "green"), las=1, lwd=2, mark.time=TRUE)
legend(18,0.95, legend=c("0","1","2"), lty=1, lwd=2, col=c("red","blue", "green"), bty="", cex=0.6)

# now fit the Cox-PH model, and use both Over 40 and MM level
cox.mod <-coxph( Surv(stan$time, stan$status) ~ stan$over40 + stan$mismatchlevel )
# and ask for a summary of the model
summary(cox.mod)
# note: there is no coef for "incercept" because COX-PH model does NOT estimate "baseline hazard"
# but still estimates hazard ratios.  REMINDER: THis cannot be used to predict survival
# How to interpret: 
#       Note 1.7249 in exp(coef) column for stan$over401. 
#       This is called the "exponentiated coefficient. This eans that at any given instant in time,
#       someone who is over 40 is 1.72 times as likely to die as someone who is under 40 adjusting for mismatch level.
#       This is the hazard ratio.
#       (Adjusting for mismatch level means if you take two individuals at SAME mismatch level, the one over 40 is 1.72 times as likely to die)
#       Subtracting 1 from the hazard ratio gives us percent MORE likely to die, so in this case 72% more likely to die as someone under 40
#
#       exp(-coef) is the reciprical of the hazard ratio and just changes which group we're talking about.  
#       So exp(-coef) of 0.5797 means that individuals UNDER 40 are 0.5797 times as likely to die as over 40
#
#       Bottom stats (3 tests) are testing overall model significance 
#         - testing null hypothesis that beta 1 = beta2 = ... = 0 vs alternative that at least 1 != 0
#       Concordance gives goodness of fit analysis. Also called c statistic. Equivalent to area under the curve.
#           - tells how well predictions match observations
#           - 0.5 is flip of coin for binomial statistics



# we can use Likelihood Ratio Test to compare nested models 
# - helps us determine if dropping or adding a variable changes the result significantly
# Can use AIC or BIC to test non-nested models
cox.mod <-coxph( Surv(stan$time, stan$status) ~ stan$over40 + stan$mismatchlevel )
cox.mod2 <-coxph( Surv(stan$time, stan$status) ~ stan$over40)

anova(cox.mod2, cox.mod, test="LRT")
# p value is quite large so we can drop mismatch level without significant loss of predictive power
# (not a significant difference between two models)




## We can also include numeric, rather than categoric, X's
cox.num <- coxph( Surv(stan$time, stan$status) ~ stan$age + stan$t5)
summary(cox.num)
# interp the HR for age: HR = 1.03
# interp the HR for t5 score: HR = 1.186




## We can also include numeric, rather than categoric, X's
cox.num <- coxph( Surv(stan$time, stan$status) ~ stan$over40 + stan$t5)
summary(cox.num)
# interp the HR for age: HR = 1.03
# interp the HR for t5 score: HR = 1.186
anova(cox.mod2, cox.num, test="LRT")



#### Check linearity
# Checking linearity (for the model that used NUM X's) using MARTINGALE residuals
plot(predict(cox.num), residuals(cox.num, type="martingale"), xlab="fitted values",ylab="Martingale residuals",main="Residual Plot", las=1)
# add a line ax y=residual=0
abline(h=0)
# fit a smoother through the points
lines(smooth.spline(predict(cox.num),residuals(cox.num, type="martingale")), col="red")

## checking linearity using DEVIANCE residuals
plot(predict(cox.num),residuals(cox.num, type="deviance"), xlab="fitted values",ylab="Deviance residuals",main="Residual Plot", las=1)
#add a line ax y=residual=0
abline(h=0)
# fit a smoother through the points
lines(smooth.spline(predict(cox.num),residuals(cox.num, type="deviance")), col="red")
# appears to be a non-linearity.  We can address this include categorizing, polynomials, etc




### CHECKING PROPORTIONAL HAZARDS ASSUMPTION
#Test for prop hazards using Schoenfeld test for PH
# H0: Hazards are proportional
# Ha: Hazards are NOT proportional
# will return test for each X, and for overall model
cox.zph(cox.num)
#tests if coef for variable(X) changes over time...
# if it changes over time -> non-proportional hazard



## test Hazard assumption with plots
# if coefficient does not change over time, expect to see beta of 0
# first, can look at plots for both ratios in one screen
par(mfrow=c(2,1))
plot(cox.zph(cox.num))
# if you get error "figure margins too large, just expand the plot window

# then can look at each plot closer
par(mfrow=c(1,1))
plot(cox.zph(cox.num)[1])
abline(h=0, col=2)
#note that line at 0 not always within 95% CI.  May not be able to trust that proportional hazard assumption is met

par(mfrow=c(1,1))
plot(cox.zph(cox.num)[2])
abline(h=0, col=2)
# here the CI always includes line at zero.  Appears to truly meet proportional hazard assumption


# If proportional hazard assumption not met, can stratify on the variable that the prop. assumpt is not met for
# also can do time dependent coefficient models aka time dependent parameter models