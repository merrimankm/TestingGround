

library(survival)
library(survminer)
library("RColorBrewer")
library("viridis")

# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)


# look at KM model of data using only one variable
data.KM1 <- survfit( Surv(data1$Time, data1$Status) ~ data1$TopPerfCat, type="kaplan-meier")
plot(data.KM1, conf.int=F, xlab="Time (days)", ylab="Survival = S(t)", col=c("blue", "cyan", "aquamarine", "seagreen", "darkseagreen", "green", "forestgreen", "red"), main="KM-Model", las=1, lwd=2, mark.time=TRUE)
legend(95,0.3, legend=c("0","1","2","3","4","5"), lty=1, lwd=2, col=c("blue", "cyan", "aquamarine", "seagreen", "darkseagreen", "green", "forestgreen", "red"), bty="", cex=0.6)





#### Cox-PH models

## Model 1

cox.num1 <- coxph( Surv(inData$Time, inData$Status) ~ inData$Vol_P5)
summary(cox.num1)


# Check linearity (if model uses NUM X's) using MARTINGALE residuals
plot(predict(cox.num1), residuals(cox.num1, type="martingale"), xlab="fitted values",ylab="Martingale residuals",main="Residual Plot", las=1)
abline(h=0)
lines(smooth.spline(predict(cox.num1),residuals(cox.num1, type="martingale")), col="red")

# Check proportional hazard assumption
cox.zph(cox.num1)
ggcoxzph(cox.zph(cox.num1))




## Model 2

cox.num2 <- coxph( Surv(data1$Time, data1$Status) ~ data1$Path_T3)
summary(cox.num2)


# Check linearity (if model uses NUM X's) using MARTINGALE residuals
plot(predict(cox.num2), residuals(cox.num2, type="martingale"), xlab="fitted values",ylab="Martingale residuals",main="Residual Plot", las=1)
abline(h=0)
lines(smooth.spline(predict(cox.num2),residuals(cox.num2, type="martingale")), col="red")

# Check proportional hazard assumption
cox.zph(cox.num2)
ggcoxzph(cox.zph(cox.num2))


## Model 3

cox.num3 <- coxph( Surv(data1$Time, data1$Status) ~ data1$P5_cat + data1$MRI_T3)
summary(cox.num3)


# Check linearity (if model uses NUM X's) using MARTINGALE residuals
plot(predict(cox.num3), residuals(cox.num3, type="martingale"), xlab="fitted values",ylab="Martingale residuals",main="Residual Plot", las=1)
abline(h=0)
lines(smooth.spline(predict(cox.num3),residuals(cox.num3, type="martingale")), col="red")

# Check proportional hazard assumption
cox.zph(cox.num3)
ggcoxzph(cox.zph(cox.num3))



## Model 4

cox.num4 <- coxph( Surv(data1$Time, data1$Status) ~ data1$max3Ddiameter_cat + data1$MRI_T3 + data1$Margins+data1$Path_T3 + data1$RP_time + data1$Age)
summary(cox.num4)


# Check linearity (if model uses NUM X's) using MARTINGALE residuals
plot(predict(cox.num4), residuals(cox.num4, type="martingale"), xlab="fitted values",ylab="Martingale residuals",main="Residual Plot", las=1)
abline(h=0)
lines(smooth.spline(predict(cox.num4),residuals(cox.num4, type="martingale")), col="red")

# Check proportional hazard assumption
cox.zph(cox.num4)
ggcoxzph(cox.zph(cox.num4))



#### Compare models with univariate and multivariate models

anova(cox.num1, cox.num3, test="LRT")
# if p value is large we can drop additional variable without significant loss of predictive power
# (not a significant difference between two models)









cox.num2 <- coxph( Surv(data1$Time, data1$Status) ~ data1$TopPerformersCombined)
summary(cox.num2)

cox.num3 <- coxph( Surv(data1$Time, data1$Status) ~ data1$Vol_prost)
summary(cox.num3)


cox.num4 <- coxph( Surv(data1$Time, data1$Status) ~ data1$Vol_indexed)
summary(cox.num4)


cox.num5 <- coxph( Surv(data1$Time, data1$Status) ~ data1$Vol_total)
summary(cox.num5)

cox.num6 <- coxph( Surv(data1$Time, data1$Status) ~ data1$MRI_EPE)
summary(cox.num6)



cox.comb <- coxph( Surv(data1$Time, data1$Status) ~ data1$TopPerformersCombined + data1$MRI_T3)
summary(cox.comb)
anova(cox.num2, cox.comb, test="LRT")








#### Create Histogram to look at status vs feature
library(ggplot2)
# reverse the order of status (put 1 before 0) so that 1s will stack on top of 0s
status = factor(data1$Status,levels=c(1,0))
ggplot(data1)+geom_histogram(aes(x=Vol_P5, fill=factor(status)))


ggplot(data1)+geom_histogram(aes(x=Vol_prost, fill=factor(status)))
ggplot(data1)+geom_histogram(aes(x=Vol_total, fill=factor(status)))
ggplot(data1)+geom_histogram(aes(x=TopPerformersCombined, fill=factor(status)))

                                 