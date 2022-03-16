library(survival)
library(survminer)

g <- "T:\\MIP\\Katie_Merriman\\Project1Data\\Rstudio_data\\Pirads_score.csv"
data1 <- read.csv(g, header=TRUE, stringsAsFactors=FALSE)

#view structure of data
#str(data1)

fit <- survfit(Surv(Time, Status) ~ PIRADS_score, data = data1)
print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = "safe",
           title="MRI EPE Score - No Margins", xlab="Time", ylab="BCR Free Survival")

