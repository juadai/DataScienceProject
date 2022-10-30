library(rstudioapi)
library(pROC)
library(lme4)
library(merTools)
library(lmerTest)

setwd(dirname(getActiveDocumentContext()$path))
df <- read.csv('../dataset/excluding weeds/understorey_cleaned.csv')


df$Plot <- factor(df$Plot)
df$Treatment <- factor(df$Treatment)
df$Fenced <- df$Fenced=="True"
df$Gap <- df$Gap=="True"
df$Shade.Tolerant <- df$Shade.Tolerant=="True"
df$Weed <- df$Weed=="True"
df$Year2 <- df$Year
df$log_abund <- log(df$Abundance)

#########################
# Model Large Forb/Herb #
#########################

# LF = Large Forb/herb
df_LF <- df[df$Life.Form=="Large Forb/Herb",]

LF_m1 <- lmer(Abundance~ (Treatment + Year + Year2 + Gap + Fenced + Weed)^2 +
                (1|Plot) + (1|Species.Name), data=df_LF)
LF_m2 <- get_model(step(LF_m1))
summary(LF_m2)


plot(LF_m2)
RSS <- sum((df_LF$Abundance - fitted(LF_m2, df_LF))^2)
TSS <- sum((df_LF$Abundance - mean(df_LF$Abundance))^2)
R2 <- 1-RSS/TSS
R2
# 0.6




LF_m3 <- lmer(log_abund ~ (Treatment + Year + Year2 + Gap + Fenced + Weed)^2 +
                (1|Plot) + (1|Species.Name), data=df_LF)
LF_m4 <- get_model(step(LF_m3))
summary(LF_m4)
plot(LF_m4)

RSS <- sum((df_LF$log_abund - fitted(LF_m4, df_LF))^2)
TSS <- sum((df_LF$log_abund - mean(df_LF$log_abund))^2)
R2 <- 1-RSS/TSS
R2
# 0.73




#############################
# Model Medium Tufted Grass #
#############################

df_MTG <- df[df$Life.Form=="Medium Tufted grass/sedge",]

MTG_m1 <- lmer(Abundance~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                (1|Plot) + (1|Species.Name), data=df_MTG)
MTG_m2 <- get_model(step(MTG_m1))
summary(MTG_m2)

plot(MTG_m2)
RSS <- sum((df_MTG$Abundance - fitted(MTG_m2, df_MTG))^2)
TSS <- sum((df_MTG$Abundance - mean(df_MTG$Abundance))^2)
R2 <- 1-RSS/TSS
R2
# 0.43

MTG_m3 <- lmer(log_abund ~ (Treatment + Year + Year2 + Gap + Fenced + Weed)^2 +
                (1|Plot) + (1|Species.Name), data=df_MTG)
MTG_m4 <- get_model(step(MTG_m3))
summary(MTG_m4)
plot(MTG_m4)

RSS <- sum((df_LF$log_abund - fitted(LF_m4, df_LF))^2)
TSS <- sum((df_LF$log_abund - mean(df_LF$log_abund))^2)
R2 <- 1-RSS/TSS
R2
# 0.73


#################################
# Model Medium Non-tufted Grass #
#################################

df_MNTG <- df[df$Life.Form=="Medium Non-tufted grass/sedge",]

MNTG_m1 <- lmer(Abundance~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                 (1|Plot) + (1|Species.Name), data=df_MNTG)
MNTG_m2 <- get_model(step(MNTG_m1))
summary(MNTG_m2)

plot(MNTG_m2)
RSS <- sum((df_MNTG$Abundance - fitted(MNTG_m2, df_MNTG))^2)
TSS <- sum((df_MNTG$Abundance - mean(df_MNTG$Abundance))^2)
R2 <- 1-RSS/TSS
R2
# 0.16


MNTG_m3 <- lmer(log_abund ~ (Treatment + Year + Year2 + Gap + Fenced + Weed)^2 +
                 (1|Plot) + (1|Species.Name), data=df_MNTG)
MNTG_m4 <- get_model(step(MNTG_m3))
summary(MNTG_m4)
plot(MNTG_m4)

RSS <- sum((df_MNTG$log_abund - fitted(MNTG_m4, df_MNTG))^2)
TSS <- sum((df_MNTG$log_abund - mean(df_MNTG$log_abund))^2)
R2 <- 1-RSS/TSS
R2
# 0.23




