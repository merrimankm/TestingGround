# Example from https://www.youtube.com/watch?v=6_AF9mMuk9E&list=PLqzoL9-eJTNDdnKvep_YHIwk2AMqHhuJ0&index=5




# Read in AIDSsurv.txt dataset
#AIDSsurv <- read.table(file="~/Documents/StatisticsClasses/KMexample.txt")
#names(AIDSsurv)
#attach(AIDSsurv)
# this data just looks at survival of a group... no x variables at all

## If I had the text file in the video, this would print out [1] "Time" "Death"
## Since I don't have the text file, the following three lines of code are used to create the dataset
Time<- c(2,3,6,6,7,10,15,15,16,27,30,32)
Death<-c(1,0,1,1,1,0,1,1,1,1,1,1)
AIDSsurv<-data.frame(cbind(Time, Death))

# NOTE: Death=1 means event occurred, Death=0 means no event (they're censored)
# When fitting a model, make sure you KNOW what 0 and 1 code for



# to fit survival models, we will need to use survival library
#this is built into R so no need to install, just load
library(survival)
km.model <- survfit(Surv(AIDSsurv$Time,AIDSsurv$Death) ~ 1 ,
                    type = "kaplan-meier")
# The ~1 is there since we have NO x-variables
# Kaplan-meier is default so we wouldn't have to specify that

#ask for some summaries of the model
km.model
# this gives the MEDIAN survival and a CI for it

#ask for a model summary
summary(km.model)
# this shows the table we produced by hand (in OneNote)
# there are no y-intercept values for this model as it's just a non-parametric step function

#we can look at a plot of this model
#set the conf.int=T if you want CI returned
#set mark.time=TRUE if you want censor marks
plot(km.model, conf.int=F, xlab="Time (months)", ylab="%Alive = S(t)", main="KM-Model", las=1, mark.time=TRUE)
#When plotting with CI, can see why CI was (7,inf
abline(h=0.5,col="red")


#####  NOW, we'll work with data that separates individuals based on whether they're over 40 years old
# remove old dataset
rm(AIDSsurv)

# create new dataset
Time<-c(2,3,6,6,7,10,15,15,16,27,30,32,1,1,1,1,2,3,3,9,22)
Death<-c(1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0)
Over40<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
daataa<-data.frame(cbind(Time,Death,Over40))

#create model
km.model2<-survfit(Surv(daataa$Time, daataa$Death) ~ Over40,
                   type="kaplan-meier")

# plot and look at summaries
plot(km.model2, conf.int=F, xlab="Time (months)", ylab="%Alive = S(t)", main="KM-Model", las=1, mark.time=TRUE)
summary(km.model2)
km.model2

#now let's add colors and a legend
#lwd is line width. 2 means 2x the width
#las=1 rotates the values on the y axis (just rotates labels so they're easier to read)
plot(km.model2, conf.int=F, xlab="Time (months)", ylab="%Alive = S(t)", main="KM-Model", col=c("red", "blue"), las=1, lwd=2, mark.time=TRUE)
# Place upper left corner at x= 18, y = 0.95
# lty is line type
# lwd is line width
# bty is box type.  Setting to N creates no box.  Leaving blank creates box
# cex sets font to percent of default, in this case 60%
legend(18,0.95, legend=c("under40","Over40"), lty=1, lwd=2, col=c("red","blue"), bty="", cex=0.6)


# do LOG-RANK-TEST to see if significant difference between groups
#     Ho: survival in two groups is same
#     Ha: surv not same between two groups
survdiff( Surv(daataa$Time,daataa$Death)~Over40 )
# refect Ho... the survival for the over40 and under 40 groups differ!


#also known as mantel-haenszel or cochran-mantel-haenszel method
# different option is wilcoxon test