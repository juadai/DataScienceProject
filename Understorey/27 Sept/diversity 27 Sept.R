
# Load data
library(rstudioapi)
library(pROC)
library(lme4)
library(merTools)
library(lmerTest)

setwd(dirname(getActiveDocumentContext()$path))
diversity_df <- read.csv('../dataset/tables/species_divers_quadrats.csv')
diversity_df_plot <- read.csv('../contingency tables/output/diversity/species_divers_plots.csv')
diversity_df_treat <- read.csv('../contingency tables/output/diversity/species_divers_treatments.csv')


# Treatments, Plot numbers, Quadrat numbers, Life forms as factors
diversity_df$Treatment <- factor(diversity_df$Treatment)
diversity_df$Plot.Number <- factor(diversity_df$Plot.Number)
diversity_df$Quadrat.Number <- factor(diversity_df$Quadrat.Number)

# Convert Fenced, Gap columns to boolean
diversity_df$Gap <- diversity_df$Gap=='True'
diversity_df$Fenced <- diversity_df$Fenced=='True'

diversity_df_plot$Gap <- diversity_df_plot$Gap=='True'
diversity_df_plot$Fenced <- diversity_df_plot$Fenced=='True'

# Remove null rows
diversity_df <- na.omit(diversity_df)
diversity_df_plot <- na.omit(diversity_df_plot)



##############################
# Diversity Modeling: Y3 ~ N #
##############################

str(diversity_df)

diversity_Y3_01 <- glm(X3 ~ (X0 + Treatment + Gap + Fenced)^2, 
                       data=diversity_df)
summary(diversity_Y3_01)
diversity_Y3_02 <- step(diversity_Y3_01)
summary(diversity_Y3_02)

anova(diversity_Y3_02, diversity_Y3_01, test='Chi')
plot(diversity_Y3_02)

k = log(length(diversity_df[,1]))
diversity_Y3_03 <- step(diversity_Y3_01, k=k)
summary(diversity_Y3_03)

diversity_Y3_02$call
diversity_Y3_03$call
anova(diversity_Y3_03, diversity_Y3_02, test='Chi')
# Can't take the smaller model - Treatment, Treatment:Gap, Fenced are significant
summary(diversity_Y3_02)
plot(diversity_Y3_02)

RSS <- sum((diversity_df$X3 - fitted(diversity_Y3_02, diversity_df))^2)
TSS <- sum((diversity_df$X3 - mean(diversity_df$X3))^2)
R2 <- 1-RSS/TSS
R2


# TRYING TO SHOW THAT THE MODEL IS NOT ADEQUATE
X1 <- data.frame(
  'Treatment'='Control',
  'Gap'=FALSE,
  'Fenced'=FALSE,
  'X0'=seq(0,3,0.1)
)

point_est <- predict(diversity_Y3_02, newdata=X1, se.fit=T)
lower_ci = point_est$fit - 1.96*pr$se.fit
upper_ci = point_est$fit + 1.96*pr$se.fit


plot(seq(0,3,0.1),point_est$fit, type='l', ylim=c(0,3),
     main='Control Plot - Estimated Diversity',
     xlab='Diversity at Year 0', ylab='Diversity at Year 6')
lines(seq(0,3,0.1), lower_ci, lty='dashed')
lines(seq(0,3,0.1), upper_ci, lty='dashed')

x <- diversity_df[diversity_df$Treatment=='Control' &
                    diversity_df$Gap ==F &
                    diversity_df$Fenced==F,]$X0

y <- diversity_df[diversity_df$Treatment=='Control' &
                    diversity_df$Gap ==F &
                    diversity_df$Fenced==F,,]$X3
points(x,y, cex=0.5)




# Try incorporating mixed effects


mixed_Y3_01 <- lmer(X3 ~ X0 + Treatment + Gap + Fenced + X0:Gap + Treatment:Gap +
                      (1 | Plot.Number), data = diversity_df)
summary(mixed_Y3_01)
anova(mixed_Y3_01, diversity_Y3_02)
colnames(diversity_df)
coef(mixed_Y3_01)

plot(mixed_Y3_01)

mixed_Y3_02 <- get_model(step(mixed_Y3_01)

RSS <- sum((diversity_df$X3 - fitted(mixed_Y3_02, diversity_df))^2)
TSS <- sum((diversity_df$X3 - mean(diversity_df$X3))^2)
R2 <- 1-RSS/TSS
R2

##############################
# Diversity Modeling: Y6 ~ N #
##############################

str(diversity_df)

diversity_Y6_01 <- glm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, 
                       data=diversity_df)
summary(diversity_Y6_01)
diversity_Y6_02 <- step(diversity_Y6_01)
summary(diversity_Y6_02)

anova(diversity_Y6_02, diversity_Y6_01, test='Chi')
plot(diversity_Y6_02)

k = log(length(diversity_df[,1]))
diversity_Y6_03 <- step(diversity_Y6_01, k=k)
summary(diversity_Y6_03)

anova(diversity_Y6_03, diversity_Y6_02, test='Chi')
# Can't take the smaller model

summary(diversity_Y6_02)


RSS <- sum((diversity_df$X6 - fitted(diversity_Y6_02, diversity_df))^2)
TSS <- sum((diversity_df$X6 - mean(diversity_df$X6))^2)
R2 <- 1-RSS/TSS
R2


diversity_df

diversity_0 <- cbind(diversity_df[1:5],diversity_df[6])
diversity_0['Year']=0
names(diversity_0) <- c("Treatment", "Plot.Number", 'Quadrat.Number', "Fenced", "Gap", "X", 'Year')

diversity_3 <- cbind(diversity_df[1:5],diversity_df[7])
diversity_3['Year'] <- 3
names(diversity_3) <- c("Treatment", "Plot.Number", 'Quadrat.Number', "Fenced", "Gap", "X", 'Year')


diversity_6 <- cbind(diversity_df[1:5],diversity_df[8])
diversity_6['Year'] <- 6
names(diversity_6) <- c("Treatment", "Plot.Number", 'Quadrat.Number', "Fenced", "Gap", "X", 'Year')

diversity <- rbind(diversity_0, diversity_3, diversity_6)
diversity$Year2 <- diversity$Year^2

# Expand the basis to be quadratic in YEAR
# Allows increase at Y3, decrease at Y6


m1 <- lmer(X ~ (Treatment + Fenced + Gap + Year + Year2)^2 +
             (1|Plot.Number), data=diversity)
summary(m1)
m2 <- get_model(step(m1))
summary(m2)
plot(m2)

RSS <- sum((diversity$X - fitted(m2, diversity))^2)
TSS <- sum((diversity$X - mean(diversity$X))^2)
R2 <- 1-RSS/TSS
R2

# Not good - residuals are too large, very little correlation

