setwd("~/Desktop/capstone/analysis/Juan")

CWDdata = read.csv("output.csv")

# model using gaussian
model = glm(difference ~ treatment, family = gaussian, data = CWDdata)
summary(model)

par(mfrow = c(2,2))
plot(model)


CWDdatarm3 = CWDdata[-c(3, 12),]
CWDdatarm6 = CWDdata[-c(6, 12),]
CWDdatarm36 = CWDdata[-c(3, 6, 12),]
CWDdata["abs_diff"] = abs(CWDdata["difference"] )

#guassian model + plot3 removed

modelrm = glm(difference ~ treatment , family = gaussian, data = CWDdatarm3)
summary(modelrm)

par(mfrow = c(2,2))
plot(modelrm)
# the plot here shows that the plot 3, 5, 9 are always outliers, but all of 
# them are from the radial treatment



# model using poisson + abs_diff
model0 = glm(abs_diff ~ treatment, family = poisson, data = CWDdata)
summary(model0)

par(mfrow = c(2,2))
plot(model0)

# model using poisson + difference
model1 = glm(difference ~ treatment, family = poisson, data = CWDdatarm6)
summary(model1)

par(mfrow = c(2,2))
plot(model1)


## quasipoisson + diff
model2 = glm(difference ~ treatment, family = quasipoisson, data = CWDdatarm6)
summary(model2)

par(mfrow = c(2,2))
plot(model2)


## Gamma + difference
model3 = glm(difference ~ treatment, family = Gamma(link = "inverse"), data = CWDdatarm6)
summary(model3)

par(mfrow = c(2,2))
plot(model3)



model4 = glm(difference ~ treatment, family = Gamma(link = "log"), data = CWDdatarm6)
summary(model4)
par(mfrow = c(2,2))
plot(model4)


## Gamma + abs_diff
model5 = glm(abs_diff ~ treatment, family = Gamma(link = "identity"), data = CWDdatarm6)
summary(model5)

par(mfrow = c(2,2))
plot(model5)


