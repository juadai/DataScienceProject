library(lme4)
library(merTools)

data = read.csv("/Users/nickjolly/Desktop/tfn/DataScienceProject/Tree DBH/dataset/all_plots.csv", stringsAsFactors = FALSE)
DBHdata = data[1:9]
DBHdata = DBHdata[!is.na(DBHdata$DBH.year.6), ]
DBHdata = DBHdata[!is.na(DBHdata$DBH.year.0), ]
DBHdata$Treatment = factor(DBHdata$Treatment)
DBHdata$Large = factor(DBHdata$Large)
DBHdatarm = DBHdata[DBHdata$DBH.year.6/DBHdata$DBH.year.0 < 1.25 & DBHdata$DBH.year.6/DBHdata$DBH.year.0 > 0.9 & DBHdata$DBH.year.0 < 80, ]
DBHdatarm$InGap = (DBHdatarm$Treatment == "T2: Radial" & DBHdatarm$Large == "Yes")
DBHdatarm$InGap = factor(DBHdatarm$InGap)

mixed2.0 = lmer(log(DBH.year.6)~log(DBH.year.0)+Treatment+InGap+log(DBH.year.0):InGap+(1|Plot), data = DBHdatarm)

X_ingap <- data.frame('DBH.year.0'=30:100, 'InGap'=TRUE, 'Plot'=6, 'Treatment'='T2: Radial')
PI_2 <- predictInterval(merMod = mixed2.0, newdata = X_ingap, level = 0.95, n.sims = 1000, stat = "median", type="linear.prediction", include.resid.var = TRUE)

X <- data.frame('DBH.year.0'=1:100, 'InGap'=FALSE, 'Plot'=6, 'Treatment'='Control')
PI <- predictInterval(merMod = mixed2.0, newdata = X, level = 0.95, n.sims = 1000, stat = "median", type="linear.prediction", include.resid.var = TRUE)

plot(1:100, exp(PI$fit), type='n', xlim=c(20,60), ylim = c(15,60), xlab = "DBH: Year 0", ylab = "Estimated DBH: Year 6",
     main="Tree Size Model: 'In Gap' vs Control")

polygon(c(1:100, 100:1), c(exp(PI$upr), rev(exp(PI$lwr))), col = "light grey", border = "light grey")
lines(1:100, exp(PI$fit), type = "l")
polygon(c(30:100, 100:30), c(exp(PI_2$upr), rev(exp(PI_2$lwr))), col = rgb(255,0,0, 50, maxColorValue = 255), border = rgb(255,0,0, 50, maxColorValue = 255))
lines(30:100, exp(PI_2$fit), type = "l", col = "red")
legend("topleft", legend = c('Control',"T2: Radial (In gap)", '95% CI'),
       lty=c(1,1,NA),
       col= c('black','red','green',NA), bty='n')






for (t in c("Control", "T1: Gap", "T2: Radial")) {
  X <- data.frame('DBH.year.0'=1:100, 'InGap'=FALSE, 'Plot'=6, 'Treatment'=t)
  PI <- predictInterval(merMod = mixed2.0, newdata = X, level = 0.95, n.sims = 1000, stat = "median", type="linear.prediction", include.resid.var = TRUE)
  
  plot(1:100, exp(PI$fit), type='n', xlim=c(20,60), ylim = c(10,80), xlab = "DBH.year.0", ylab = "predicted DBH.year.6")
  polygon(c(1:100, 100:1), c(exp(PI$upr), rev(exp(PI$lwr))), col = "light grey", border = "light grey")
  lines(1:100, exp(PI$fit), type = "l")
  polygon(c(1:100, 100:1), c(exp(PI_2$upr), rev(exp(PI_2$lwr))), col = rgb(255,0,0, 50, maxColorValue = 255), border = rgb(255,0,0, 50, maxColorValue = 255))
  lines(1:100, exp(PI_2$fit), type = "l", col = "red")
  legend("topleft", legend = c("In Gap", t), fill = c("red", "black"))
}


?legend


par(bg="#fffdf3")
treatments = c("Control", "T1: Gap", "T2: Radial")
shades = c("light grey", rgb(255,0,0, 50, maxColorValue = 255), rgb(0,255,0, 50, maxColorValue = 255))
colours = c("black", "red", "green")

plot(1:100, exp(PI$fit), type='n', xlim=c(20,60), ylim = c(15,60), xlab = "DBH: Year 0", ylab = "Estimated DBH: Year 6",
     main='Tree Size Model: Treatments vs Control')
for (i in 1:3) {
  X <- data.frame('DBH.year.0'=1:100, 'InGap'=FALSE, 'Plot'=6, 'Treatment'=treatments[i])
  PI <- predictInterval(merMod = mixed2.0, newdata = X, level = 0.95, n.sims = 1000, stat = "median", type="linear.prediction", include.resid.var = TRUE)
  
  polygon(c(1:100, 100:1), c(exp(PI$upr), rev(exp(PI$lwr))), col = shades[i], border = shades[i])
  lines(1:100, exp(PI$fit), type = "l", col = colours[i])
}
legend("topleft", legend = c("Control", "T1: Gap", "T2: Radial (Not in gap)", '95% CI'), 
       lty=c(1,1,1,NA),
       col= c('black','red','green',NA), bty='n')


library(car)

summary(mixed2.0)

C <- c(0,0,1,0,0,0)
linearHypothesis(mixed2.0, C)

C <- c(0,0,0,0,1,1)
linearHypothesis(mixed2.0, C)
