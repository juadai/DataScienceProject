library(rstudioapi)
library(pROC)
library(lme4)
library(merTools)
library(lmerTest)

setwd(dirname(getActiveDocumentContext()$path))
diversity_df <- read.csv('../dataset/excluding weeds/diversity_quadrat.csv')

diversity_df$Treatment <- factor(diversity_df$Treatment)
diversity_df$Plot <- factor(diversity_df$Plot)
diversity_df$Quadrat <- factor(diversity_df$Quadrat)

# Convert Fenced, Gap columns to boolean
diversity_df$Gap <- diversity_df$Gap=='True'
diversity_df$Fenced <- diversity_df$Fenced=='True'

diversity_df <- na.omit(diversity_df)


##############################
# Diversity Modeling: Y3 ~ N #
##############################

diversity_Y3_01 <- glm(X3 ~ (X0 + Treatment + Gap + Fenced)^2, 
                       data=diversity_df)
summary(diversity_Y3_01)
diversity_Y3_02 <- step(diversity_Y3_01)
summary(diversity_Y3_02)

plot(diversity_Y3_02)

RSS <- sum((diversity_df$X3 - fitted(diversity_Y3_02, diversity_df))^2)
TSS <- sum((diversity_df$X3 - mean(diversity_df$X3))^2)
R2 <- 1-RSS/TSS
R2
# Very low, 0.17


##############################
# Diversity Modeling: Y6 ~ N #
##############################

diversity_Y6_01 <- glm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, 
                       data=diversity_df)
summary(diversity_Y6_01)
diversity_Y6_02 <- step(diversity_Y6_01)
summary(diversity_Y6_02)

plot(diversity_Y6_02)

RSS <- sum((diversity_df$X6 - fitted(diversity_Y6_02, diversity_df))^2)
TSS <- sum((diversity_df$X6 - mean(diversity_df$X6))^2)
R2 <- 1-RSS/TSS
R2
# Very low, 0.18





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


m1 <- lmer(X ~ (Treatment + Fenced + Gap + Year + Year2)^2 +
             (1|Plot.Number), data=diversity)
m2 <- get_model(step(m1))
summary(m2)

RSS <- sum((diversity$X - fitted(m2, diversity))^2)
TSS <- sum((diversity$X - mean(diversity$X))^2)
R2 <- 1-RSS/TSS
R2
# Low - 0.22
