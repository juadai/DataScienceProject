library(rstudioapi)
library(pROC)
library(lme4)
library(merTools)
library(lmerTest)
library(dplyr)

setwd(dirname(getActiveDocumentContext()$path))
abundance_df <- read.csv('../dataset/excluding weeds/lf_abund_quadrat.csv')

grouped <- aggregate(abundance_df$Abundance, 
                     by=list(abundance_df$Treatment, abundance_df$Plot,
                             abundance_df$Quadrat, abundance_df$Fenced,
                             abundance_df$Gap, abundance_df$Year),
                     FUN=sum)
colnames(grouped) <- c('Treatment', 'Plot', 'Quadrat', 'Fenced', 'Gap', 'Year',
                       'Abundance')

grouped$Year2 <- grouped$Year^2
grouped$log_abund <- log(grouped$Abundance)



model1 <- lmer(Abundance~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                  (1|Plot), data=grouped)
model2 <- get_model(step(model1))
summary(model2)

plot(model2)
RSS <- sum((grouped$Abundance - fitted(model2, grouped))^2)
TSS <- sum((grouped$Abundance - mean(grouped$Abundance))^2)
R2 <- 1-RSS/TSS
R2
# 0.29


model3 <- lmer(log_abund ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                  (1|Plot), data=grouped)
model4 <- get_model(step(model3))
summary(model4)
plot(model4)

RSS <- sum((grouped$log_abund - fitted(model4, grouped))^2)
TSS <- sum((grouped$log_abund - mean(grouped$log_abund))^2)
R2 <- 1-RSS/TSS
R2
# 0.27



grouped <- grouped[grouped$log_abund>1,]

model3 <- lmer(log_abund ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                 (1|Plot), data=grouped)
model4 <- get_model(step(model3))
summary(model4)
plot(model4)

RSS <- sum((grouped$log_abund - fitted(model4, grouped))^2)
TSS <- sum((grouped$log_abund - mean(grouped$log_abund))^2)
R2 <- 1-RSS/TSS
R2
# 0.31