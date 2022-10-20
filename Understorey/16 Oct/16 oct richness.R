library(rstudioapi)
library(pROC)

library(lme4)
library(lmerTest)


setwd(dirname(getActiveDocumentContext()$path))
richness_df <- read.csv('../dataset/excluding weeds/richness_quadrat.csv')

# Treatments, Plot numbers, Quadrat numbers, Life forms as factors
richness_df$Treatment <- factor(richness_df$Treatment)
richness_df$Plot.Number <- factor(richness_df$Plot.Number)
richness_df$Quadrat.Number <- factor(richness_df$Quadrat.Number)

richness_df$Gap <- richness_df$Gap=='True'
richness_df$Fenced <- richness_df$Fenced=='True'

# Remove null rows
richness_df <- na.omit(richness_df)

########################################
# Richness Modeling: Year 3: Y3 ~ Pois #
########################################

richness_Y3_01 <- glm(X3 ~ (X0 + Treatment + Gap + Fenced)^2, 
                      family='poisson', data=richness_df)
summary(richness_Y3_01)

richness_Y3_02 <- step(richness_Y3_01)
summary(richness_Y3_02)

# Excluding weeds - we lost the term associated with browsing

pchisq(263.32, 260, lower=F)
# Passes goodness-of-fit test

library(car)
C <- matrix(c(0,0,1,,1,1,0), nrow=1)
linearHypothesis(richness_Y3_02, C)
# No difference between (T1, In gap) and (Control, Not in gap)

C <- matrix(c(0,0,0,1,1,0,1), nrow=1)
linearHypothesis(richness_Y3_02, C)
# No difference between (T2, In gap) and (Control, Not in gap)

C <- matrix(c(0,0,1,-1,0,1,-1), nrow=1)
linearHypothesis(richness_Y3_02, C)
# No difference between treatments

########################################
# Richness Modeling: Year 3: Y6 ~ Pois #
########################################

richness_Y6_01 <- glm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, 
                      family='poisson', data=richness_df)
summary(richness_Y6_01)

richness_Y6_02 <- step(richness_Y6_01)
summary(richness_Y6_02)

pchisq(313.82, 255, lower=F)
# Does not pass goodness-of-fit test






richness_0 <- cbind(richness_df[1:5],richness_df[6])
richness_0['Year']=0
names(richness_0) <- c("Treatment", "Plot.Number", 'Quadrat.Number', "Fenced", "Gap", "X", 'Year')

richness_3 <- cbind(richness_df[1:5],richness_df[7])
richness_3['Year'] <- 3
names(richness_3) <- c("Treatment", "Plot.Number", 'Quadrat.Number', "Fenced", "Gap", "X", 'Year')

richness_6 <- cbind(richness_df[1:5],richness_df[8])
richness_6['Year'] <- 6
names(richness_6) <- c("Treatment", "Plot.Number", 'Quadrat.Number', "Fenced", "Gap", "X", 'Year')

richness <- rbind(richness_0, richness_3, richness_6)
richness$Year2 <- richness$Year^2

model1 <- lmer(X ~ (Treatment + Fenced + Gap + Year + Year2)^2 +
                 (1|Plot.Number), data=richness)
model2 <- get_model(step(model1))
summary(model2)

RSS <- sum((richness$X - fitted(model2, richness))^2)
TSS <- sum((richness$X - mean(richness$X))^2)
R2 <- 1-RSS/TSS
R2
# Low R^2 score, 0.25

