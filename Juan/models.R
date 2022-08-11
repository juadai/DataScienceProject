setwd("~/Desktop/capstone/analysis/Juan")

CWDdata = read.csv("output.csv")

# model using gaussian
model = glm(difference ~ treatment, family = gaussian, data = CWDdata)
summary(model)

par(mfrow = c(1,1))
plot(model)

# log guassian model
CWDdata = CWDdata[-c(12)]
modellog = glm(log(difference) ~ treatment , family = gaussian, data = CWDdata)
summary(modellog)


par(mfrow = c(2,2))
plot(modellog)

# the plot here shows that the plot 3, 5, 9 are always outl, but all of 
# them are from the radial treatment

## quasipoisson + abs_diff
CWDdata["abs_diff"] = abs(CWDdata["difference"] )
model1 = glm(abs_diff ~ treatment, family = quasipoisson(link = "log"), data = CWDdata)
summary(model1)

par(mfrow = c(2,2))
plot(model1)

# model using poisson + abs_diff
model2 = glm(abs_diff ~ treatment, family = poisson, data = CWDdata)
summary(model2)

par(mfrow = c(2,2))
plot(model2)




## Gamma + abs_diff
model3 = glm(abs_diff ~ treatment, family = Gamma(link = "inverse"), data = CWDdata)
summary(model3)

par(mfrow = c(1,1))
plot(model3)

