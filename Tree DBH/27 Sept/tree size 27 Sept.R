# Load data
library(rstudioapi)
library(pROC)
library(lme4)

setwd(dirname(getActiveDocumentContext()$path))

df <- read.csv('../dataset/all_plots.csv')

colnames(df)
df <- df[,c(1,2,3,4,6,7)]

df[df$Large=='', 'Large'] <- 'No'
df <- df[!is.na(df$DBH.year.6),]
df <- df[!is.na(df$DBH.year.0),]

df$Treatment <- factor(df$Treatment)
df$Plot <- factor(df$Plot)
df$Large <- factor(df$Large)

df$log.DBH.year.6 <- log(df$DBH.year.6)
df$log.DBH.year.0 <- log(df$DBH.year.0)

rownames(df) <- 1:length(df[,1])


model1 <- glm(log.DBH.year.6 ~ (log.DBH.year.0+Treatment+Large)^2, data=df)
summary(model1)
plot(model1)
# Step function - all terms remain with AIC and BIC

outliers_model1<- model1$residuals[abs(model1$residuals)>0.15]
outliers_model1
df[names(outliers_model1),]


# Try mixed effect model - random indercept at the plot level
model2 <- lmer(log.DBH.year.6 ~ (log.DBH.year.0 + Treatment + Large)^2 +
               (1|Plot), data=df)
summary(model2)


df$Gap <- FALSE
df[(df$Treatment=='T2: Radial') & (df$Large=='Yes'), 'Gap'] <- TRUE

model3 <- lmer(log.DBH.year.6 ~ log.DBH.year.0 + Treatment + Gap + 
                 log.DBH.year.0:Treatment + log.DBH.year.0:Gap +
                 (1|Plot), data=df)
summary(model3)

model4 <- lmer(log.DBH.year.6 ~ log.DBH.year.0 + Treatment + Gap + 
                 + log.DBH.year.0:Gap + (1|Plot), data=df)
summary(model4)
anova(model4, model3) #Y0:Treatment can be dropped


model5 <- lmer(log.DBH.year.6 ~ log.DBH.year.0 + Treatment + Gap + 
                 log.DBH.year.0:Treatment + (1|Plot), data=df)
summary(model5)
anova(model5, model3) # Y0:Gap cannot be dropped

plot(model4)

qqnorm(residuals(model4))
qqline(residuals(model4))

outliers_model4 <- residuals(model4)[abs(residuals(model4))>0.15]
df[names(outliers_mixed),]

rss <- sum((df$log.DBH.year.6-fitted(model4))^2)
tss <- sum((df$log.DBH.year.6-mean(df$log.DBH.year.6))^2)
1-rss/tss
