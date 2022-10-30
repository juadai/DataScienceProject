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

final_model = lmer(log(DBH.year.6)~log(DBH.year.0)+Treatment+InGap+log(DBH.year.0):InGap+(1|Plot), data = DBHdatarm)
summary(final_model)



######################
# Hypothesis testing #
######################
library(car)

# H0:
# Treatment Radial = Treatment Gap

C <- c(0,0,1,-1,0,0)
linearHypothesis(final_model,C)


# H0:
# Treatment Radial + In Gap + In Gap:log(DBH Year 0) = 0

C <- c(0,0,0,1,1,log(50))
linearHypothesis(final_model,C)

C <- c(0,0,0,1,1,log(54))
linearHypothesis(final_model,C)

C <- c(0,0,0,1,1,log(55))
linearHypothesis(final_model,C)

C <- c(0,0,0,1,1,log(57))
linearHypothesis(final_model,C)

C <- c(0,0,0,1,1,log(59))
linearHypothesis(final_model,C)

C <- c(0,0,0,1,1,log(70))
linearHypothesis(final_model,C)


# H0:
# In Gap + In Gap:log(DBH Year 0) = 0

C <- c(0,0,0,0,1,log(40))
linearHypothesis(final_model,C)

C <- c(0,0,0,0,1,log(45))
linearHypothesis(final_model,C)

C <- c(0,0,0,0,1,log(42))
linearHypothesis(final_model,C)

C <- c(0,0,0,0,1,log(44))
linearHypothesis(final_model,C)

C <- c(0,0,0,0,1,log(46))
linearHypothesis(final_model,C)

C <- c(0,0,0,0,1,log(47))
linearHypothesis(final_model,C)


##################
# Visualisations #
##################



# Fitted values for Control
X <- data.frame('DBH.year.0'=1:60,
                'InGap'=FALSE,
                'Plot'=3,
                'Treatment'='Control')

PI <- predictInterval(merMod = final_model, newdata = X,
                      level = 0.95, n.sims = 1000,
                      stat = "median", type="linear.prediction",
                      include.resid.var = TRUE)

# Plot Control vs. T1 (Not large trees)
X_T1 <- data.frame('DBH.year.0'=1:60,
                'InGap'=FALSE,
                'Plot'=3,
                'Treatment'='T1: Gap')

PI_T1 <- predictInterval(merMod = final_model, newdata = X_T1,
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

PI_T2 <- predictInterval(merMod = final_model, newdata = X_T1,
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

PI_ingap <- predictInterval(merMod = final_model, newdata = X_ingap,
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

alpha <- coef(final_model)$Plot[,'(Intercept)']
summary(alpha)
exp(coef(final_model)$Plot)



# visualise predicted interval for the three treatments
treatments = c("Control", "T1: Gap", "T2: Radial")
shades = c(rgb(0,0,0, 25, maxColorValue = 255), rgb(126, 0, 222, 35, maxColorValue = 255), rgb(190,0,90, 35, maxColorValue = 255))
colours = c("black", "#7e00de", "#c7005a")

plot(1:100, exp(PI$fit), type='n', xlim=c(20,60), ylim = c(15,60), xlab = "DBH: Year 0", ylab = "Estimated DBH: Year 6",
     main='Tree Size Model: Treatments vs Control')
for (i in 1:3) {
  X <- data.frame('DBH.year.0'=1:100, 'InGap'=FALSE, 'Plot'=6, 'Treatment'=treatments[i])
  PI <- predictInterval(merMod = final_model, newdata = X, level = 0.95, n.sims = 1000, stat = "median", type="linear.prediction", include.resid.var = TRUE)
  
  polygon(c(1:100, 100:1), c(exp(PI$upr), rev(exp(PI$lwr))), col = shades[i], border = shades[i])
  lines(1:100, exp(PI$fit), type = "l", col = colours[i], lwd=1.5)
}
legend("topleft", legend = c("Control", "T1: Gap", "T2: Radial (Not in gap)", '95% CI'), 
       lty=c(1,1,1,NA),
       col= c(colours,NA), bty='n', lwd=1.5)




# visualise predicted interval for InGap VS Control
X_ingap <- data.frame('DBH.year.0'=30:100, 'InGap'=TRUE, 'Plot'=6, 'Treatment'='T2: Radial')
PI_2 <- predictInterval(merMod = final_model, newdata = X_ingap, level = 0.95, n.sims = 1000, stat = "median", type="linear.prediction", include.resid.var = TRUE)

X <- data.frame('DBH.year.0'=1:100, 'InGap'=FALSE, 'Plot'=6, 'Treatment'='Control')
PI <- predictInterval(merMod = final_model, newdata = X, level = 0.95, n.sims = 1000, stat = "median", type="linear.prediction", include.resid.var = TRUE)

plot(1
     :100, exp(PI$fit), type='n', xlim=c(20,60), ylim = c(15,60), xlab = "DBH: Year 0", ylab = "Estimated DBH: Year 6",
     main="Tree Size Model: 'In Gap' vs Control")

polygon(c(1:100, 100:1), c(exp(PI$upr), rev(exp(PI$lwr))), col = rgb(0,0,0, 25, maxColorValue = 255), border = "light grey")
lines(1:100, exp(PI$fit), type = "l", lwd=1.5)
polygon(c(30:100, 100:30), c(exp(PI_2$upr), rev(exp(PI_2$lwr))), col = rgb(190, 0, 90, 35, maxColorValue = 255), border = rgb(190, 0, 90, 60, maxColorValue = 255))
lines(30:100, exp(PI_2$fit), type = "l", col = "#c7005a", lwd=1.5)
legend("topleft", legend = c('Control',"T2: Radial (In gap)", '95% CI'),
       lty=c(1,1,NA),
       col= c('black','red','green',NA), bty='n', lwd=1.5)

