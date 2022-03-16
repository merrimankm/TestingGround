# Example from https://www.youtube.com/watch?v=6_AF9mMuk9E&list=PLqzoL9-eJTNDdnKvep_YHIwk2AMqHhuJ0&index=5



# Read Data File
g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\AllLesionData.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)


# Create K-M curve
library(survival)
km.model <- survfit(Surv(data1$Time,data1$Status) ~ data1$TopPerfCat,type = "kaplan-meier")

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

#now let's add colors and a legend
#lwd is line width. 2 means 2x the width
#las=1 rotates the values on the y axis (just rotates labels so they're easier to read)
plot(km.model, conf.int=F, xlab="Time (months)", ylab="%Alive = S(t)", main="KM-Model", col=c("red", "blue"), las=1, lwd=2, mark.time=TRUE)
# Place upper left corner at x= 18, y = 0.95
# lty is line type
# lwd is line width
# bty is box type.  Setting to N creates no box.  Leaving blank creates box
# cex sets font to percent of default, in this case 60%
legend(18,0.95, legend=c("under40","Over40"), lty=1, lwd=2, col=c("red","blue"), bty="", cex=0.6)


# do LOG-RANK-TEST to see if significant difference between groups
#     Ho: survival in two groups is same
#     Ha: surv not same between two groups
survdiff( Surv(data1$Time,data1$Status)~data1$TopPerfCat )
# refect Ho... the survival for the over40 and under 40 groups differ!


#also known as mantel-haenszel or cochran-mantel-haenszel method
# different option is wilcoxon test