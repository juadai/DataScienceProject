# Load data
library(rstudioapi)
library(pROC)
library(lme4)
library(lmerTest)

setwd(dirname(getActiveDocumentContext()$path))
richness_df <- read.csv('../dataset/tables/species_richness_quadrats.csv')
richness_df_plot <- read.csv('../contingency tables/output/richness/species_richness_plots.csv')

# Treatments, Plot numbers, Quadrat numbers, Life forms as factors
richness_df$Treatment <- factor(richness_df$Treatment)
richness_df$Plot.Number <- factor(richness_df$Plot.Number)
richness_df$Quadrat.Number <- factor(richness_df$Quadrat.Number)

# Convert Fenced, Gap columns to boolean
richness_df$Gap <- richness_df$Gap=='True'
richness_df$Fenced <- richness_df$Fenced=='True'

richness_df_plot$Gap <- richness_df_plot$Gap=='True'
richness_df_plot$Fenced <- richness_df_plot$Fenced=='True'

# Remove null rows
richness_df <- na.omit(richness_df)
richness_df_plot <- na.omit(richness_df_plot)

########################################
# Richness Modeling: Year 3: Y3 ~ Pois #
########################################

str(richness_df)

richness_Y3_01 <- glm(X3 ~ (X0 + Treatment + Gap + Fenced)^2, 
                      family='poisson', data=richness_df)
summary(richness_Y3_01)

richness_Y3_02 <- step(richness_Y3_01)
summary(richness_Y3_02)
plot(richness_Y3_02)

anova(glm(formula = X3 ~ X0 + Treatment + Gap + Fenced + Treatment:Gap + 
            Treatment:Fenced, family = "poisson", data = richness_df),
      glm(formula = X3 ~ X0 + Treatment + Gap + Fenced + Treatment:Gap, family = "poisson", data = richness_df),
      test='Chi')

richness_Y3_03 <- glm(formula = X3 ~ X0 + Treatment + Gap + Fenced + Treatment:Gap, family = "poisson", data = richness_df)
summary(richness_Y3_03)
pchisq(265.356, 261, lower=F) # Pass


k = log(length(richness_df[,1]))
richness_Y3_04 <- step(richness_Y3_01, k=k)
summary(richness_Y3_04)

anova(richness_Y3_04, richness_Y3_02, test='Chi')
richness_Y3_02$call
richness_Y3_04$call
# Keep the larger model 2
summary(richness_Y3_02)

pchisq(263.32, 260, lower=F) # Pass


########################################
# Richness Modeling: Year 6: Y6 ~ Pois #
########################################

richness_Y6_01 <- glm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, 
                      family='poisson', data=richness_df)
summary(richness_Y6_01)

richness_Y6_02 <- step(richness_Y6_01)
summary(richness_Y6_02)

anova(richness_Y6_02, richness_Y6_01, test='Chi')

k = log(length(richness_df[,1]))
richness_Y6_03 <- step(richness_Y6_01, k=k)
summary(richness_Y6_03)

anova(richness_Y6_02, richness_Y6_03, test='Chi')
# Keep the larger model 2
summary(richness_Y6_02)

richness_Y6_01 <- glm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, 
                      family='poisson', data=richness_df)

pchisq(313.82, 225, lower=F) # Fail


modelY6 <- glmer(X6 ~ X0 + X3 + Treatment + Gap + Fenced + X0:X3 + 
                   X0:Fenced + X3:Treatment + X3:Fenced + Treatment:Gap +
                   (1 | Plot.Number), family = "poisson", 
                 data = richness_df)

summary(modelY6)$Plot.Number

# log(lambda) = x*beta
# lambda = exp(x*beta) 
exp(2.180927)
exp(1.689850)

coef(modelY6)
colnames(richness_df)

plot(modelY6)
plot(richness_Y6_02)

# Try negative binomial
# Try poisson without interaction


c(mean(richness_df$X0), var(richness_df$X0))
c(mean(richness_df$X3), var(richness_df$X3))
c(mean(richness_df$X6), var(richness_df$X6))

# Key assumption for Poisson distribution doesn't hold (mean = variance)

richness_df

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

# Expand the basis to be quadratic in YEAR
# Allows increase at Y3, decrease at Y6

m1 <- glm(X ~ (Treatment + Fenced + Gap + Year)^2, family='poisson', data=richness)
summary(m1)
m2 <- step(m1)
summary(m2)
plot(m2)


richness

model1 <- lmer(X ~ (Treatment + Fenced + Gap + Year + Year2)^2 +
                 (1|Plot.Number), data=richness)

model2 <- get_model(step(model1))
summary(model2)
ranova(model2)

RSS <- sum((richness$X - fitted(model2, richness))^2)
TSS <- sum((richness$X - mean(richness$X))^2)
R2 <- 1-RSS/TSS
R2

plot(model2)

# Not good - residuals are too large, very little correlation

