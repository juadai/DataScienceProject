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

X <- data.frame('DBH.year.0'=1:100,
                'InGap'=FALSE,
                'Plot'=6,
                'Treatment'='Control')


PI <- predictInterval(merMod = mixed2.0, newdata = X,
                      level = 0.95, n.sims = 1000,
                      stat = "median", type="linear.prediction",
                      include.resid.var = TRUE)

plot(1:100, exp(PI$fit), type='l',
     xlim=c(20,60))
lines(1:100, exp(PI$upr), type='l', col='grey')
lines(1:100, exp(PI$lwr), type='l', col='grey')

X_ingap <- data.frame('DBH.year.0'=1:100,
                      'InGap'=TRUE,
                      'Plot'=6,
                      'Treatment'='T2: Radial')

PI_2 <- predictInterval(merMod = mixed2.0, newdata = X_ingap,
                              level = 0.95, n.sims = 1000,
                              stat = "median", type="linear.prediction",
                              include.resid.var = TRUE)


lines(1:100, exp(PI_2$fit), type='l', col='red')
lines(1:100, exp(PI_2$upr), type='l', lty='dashed', col='red')
lines(1:100, exp(PI_2$lwr), type='l', lty='dashed', col='red')
      
