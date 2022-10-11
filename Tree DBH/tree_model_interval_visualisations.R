library(lme4)
library(merTools)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


data = read.csv("./dataset/all_plots.csv", stringsAsFactors = FALSE)

DBHdata = data[1:9]
DBHdata = DBHdata[!is.na(DBHdata$DBH.year.6), ]
DBHdata = DBHdata[!is.na(DBHdata$DBH.year.0), ]
DBHdata$Treatment = factor(DBHdata$Treatment)
DBHdata$Large = factor(DBHdata$Large)

DBHdatarm = DBHdata[DBHdata$DBH.year.6/DBHdata$DBH.year.0 < 1.25 & DBHdata$DBH.year.6/DBHdata$DBH.year.0 > 0.9 & DBHdata$DBH.year.0 < 80, ]

DBHdatarm$InGap = (DBHdatarm$Treatment == "T2: Radial" & DBHdatarm$Large == "Yes")
DBHdatarm$InGap = factor(DBHdatarm$InGap)

mixed2.0 = lmer(log(DBH.year.6)~log(DBH.year.0)+Treatment+InGap+log(DBH.year.0):InGap+(1|Plot), data = DBHdatarm)



# Fitted values for Control
X <- data.frame('DBH.year.0'=1:60,
                'InGap'=FALSE,
                'Plot'=3,
                'Treatment'='Control')

PI <- predictInterval(merMod = mixed2.0, newdata = X,
                      level = 0.95, n.sims = 1000,
                      stat = "median", type="linear.prediction",
                      include.resid.var = TRUE)

# Plot Control vs. T1 (Not large trees)
X_T1 <- data.frame('DBH.year.0'=1:60,
                'InGap'=FALSE,
                'Plot'=3,
                'Treatment'='T1: Gap')

PI_T1 <- predictInterval(merMod = mixed2.0, newdata = X_T1,
                        level = 0.95, n.sims = 1000,
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)


plot(0, xlim=c(30,50), ylim=c(30,60), main = 'Control vs T1, T2',
     xlab= 'DBH at Year 0', ylab='Fitted value for DBH at Year 6')
lines(1:60, exp(PI$fit), type='l')
lines(1:60, exp(PI$upr), type='l', col='grey')
lines(1:60, exp(PI$lwr), type='l', col='grey')

lines(1:60, exp(PI_T1$fit), type='l', col='red', cex=4)
lines(1:60, exp(PI_T1$upr), type='l', lty='dashed', col='red')
lines(1:60, exp(PI_T1$lwr), type='l', lty='dashed', col='red')

X_T2 <- data.frame('DBH.year.0'=1:60,
                   'InGap'=FALSE,
                   'Plot'=3,
                   'Treatment'='T2: Radial')

PI_T2 <- predictInterval(merMod = mixed2.0, newdata = X_T1,
                        level = 0.95, n.sims = 1000,
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)

lines(1:60, exp(PI_T2$fit), type='l', col='green', cex=4)
lines(1:60, exp(PI_T2$upr), type='l', lty='dashed', col='green')
lines(1:60, exp(PI_T2$lwr), type='l', lty='dashed', col='green')


# Plot control vs. in-gap

X_ingap <- data.frame('DBH.year.0'=33:70,
                      'InGap'=TRUE,
                      'Plot'=3,
                      'Treatment'='T2: Radial')

PI_ingap <- predictInterval(merMod = mixed2.0, newdata = X_ingap,
                              level = 0.95, n.sims = 1000,
                              stat = "median", type="linear.prediction",
                              include.resid.var = TRUE)

plot(0, xlim=c(30,50), ylim=c(30,60), main='Control vs In Gap',
     
     xlab= 'DBH at Year 0', ylab='Fitted value for DBH at Year 6')
lines(1:60, exp(PI$fit), type='l')
lines(1:60, exp(PI$upr), type='l', col='grey')
lines(1:60, exp(PI$lwr), type='l', col='grey')

lines(33:70, exp(PI_ingap$fit)
      , type='l', col='red', cex=4)
lines(33:70, exp(PI_ingap$upr), type='l', lty='dashed', col='red')
lines(33:70, exp(PI_ingap$lwr), type='l', lty='dashed', col='red')
abline(v=33, lty='dashed')
      
summary(DBHdatarm[DBHdatarm$Large=='Yes', 'DBH.year.6'])
summary(DBHdatarm[DBHdatarm$Large=='No', 'DBH.year.6'])

alpha <- coef(mixed2.0)$Plot[,'(Intercept)']
summary(alpha)
exp(coef(mixed2.0)$Plot)

bootMer()
