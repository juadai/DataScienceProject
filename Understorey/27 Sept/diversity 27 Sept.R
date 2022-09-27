
# Load data
library(rstudioapi)
library(pROC)
library(lme4)

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


mixed_Y3_01 <- lmer(X3 ~ X0 + Treatment + Gap + Fenced + X0:Gap + Treatment:Gap +
                      (1 | Plot.Number), data = diversity_df)
summary(mixed_Y3_01)
anova(mixed_Y3_01, diversity_Y3_02)
colnames(diversity_df)
coef(mixed_Y3_01)

plot(mixed_Y3_01)



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



