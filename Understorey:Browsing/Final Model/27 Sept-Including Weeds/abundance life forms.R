
library(rstudioapi)
library(lme4)
library(lmerTest)

setwd(dirname(getActiveDocumentContext()$path))
df <- read.csv('../dataset/tables/abundance_species_for_analysis.csv')
unique(df['Life.Form'])
colnames(df)

#grouped <- aggregate(cbind(df$X0,df$X3,df$X6),
#                     list(df$Life.Form), FUN=sum)
#grouped <- grouped[order(grouped$V1, decreasing=T),]
#rownames(grouped) <- grouped$Group.1

df$Plot.Number <- factor(df$Plot.Number)
df$Treatment <- factor(df$Treatment)
df$Fenced <- df$Fenced=="True"
df$Gap <- df$Gap=="True"

# Remove rows with all abundance = 0
df <- df[!(df$X0 == 0 & df$X3 == 0 & df$X6 == 0),]

df_Y0 <- df[,c(1:7,8)]
df_Y0 <- cbind(df_Y0, 0)
colnames(df_Y0)[8:9] <- c("X", 'Year')

df_Y3 <- df[,c(1:7,9)]
df_Y3 <- cbind(df_Y3, 3)
colnames(df_Y3)[8:9] <- c("X", 'Year')

df_Y6 <- df[,c(1:7,10)]
df_Y6 <- cbind(df_Y6, 6)
colnames(df_Y6)[8:9] <- c("X", 'Year')

df_prime <- rbind(df_Y0,df_Y3,df_Y6)
df_prime$Year2 <- df_prime$Year^2

df_prime$log_X <- df_prime$X
df_prime[df_prime$X==0, 'log_X'] <- 0.1
df_prime$log_X <- log(df_prime$log_X)



#########################
# Model Large Forb/Herb #
#########################

# LF = Large Forb/herb
df_LF <- df_prime[df_prime$Life.Form=="Large Forb/Herb",]


LF_m1 <- lmer(X ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                (1|Plot.Number) + (1|Species.Name), data=df_LF)
LF_m2 <- get_model(step(LF_m1))
summary(LF_m2)

plot(LF_m2)

mean(df_LF[df_LF$Year==0,'X'])
var(df_LF[df_LF$Year==0,'X'])

mean(df_LF[df_LF$Year==3,'X'])
var(df_LF[df_LF$Year==3,'X'])

mean(df_LF[df_LF$Year==6,'X'])
var(df_LF[df_LF$Year==6,'X'])



# Exclude rows where abundance = 0
df_LF <- df_LF[df_LF$log_X >= 0,]

LF_m3 <- lmer(log_X ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                (1|Plot.Number) + (1|Species.Name), data=df_LF)
LF_m4 <- get_model(step(LF_m3))
summary(LF_m4)
plot(LF_m4)

RSS <- sum((df_LF$log_X - fitted(LF_m4, df_LF))^2)
TSS <- sum((df_LF$log_X - mean(df_LF$log_X))^2)
R2 <- 1-RSS/TSS
R2



#############################
# Model Medium Tufted Grass #
#############################

# MTG = Medium Tufted Grass
df_MTG <- df_prime[df_prime$Life.Form=="Medium Tufted grass/sedge",]

MTG_m1 <- lmer(X ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                (1|Plot.Number) + (1|Species.Name), data=df_MTG)
MTG_m2 <- get_model(step(MTG_m1))

summary(MTG_m2)
plot(MTG_m2)

# Exclude rows where abundance = 0
df_MTG <- df_MTG[df_MTG$log_X >= 0,]

MTG_m3 <- lmer(log_X ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                (1|Plot.Number) + (1|Species.Name), data=df_MTG)
MTG_m4 <- get_model(step(MTG_m3))
summary(MTG_m4)
plot(MTG_m4)

RSS <- sum((df_MTG$log_X - fitted(MTG_m4, df_MTG))^2)
TSS <- sum((df_MTG$log_X - mean(df_MTG$log_X))^2)
R2 <- 1-RSS/TSS
R2



#################################
# Model Medium Non-Tufted Grass #
#################################

# MNTG = Medium Non-Tufted Grass
df_MNTG <- df_prime[df_prime$Life.Form=="Medium Non-tufted grass/sedge",]

MNTG_m1 <- lmer(X ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                 (1|Plot.Number) + (1|Species.Name), data=df_MNTG)
MNTG_m2 <- get_model(step(MNTG_m1))
summary(MNTG_m2)
plot(MNTG_m2)


# Exclude rows where abundance = 0
df_MNTG <- df_MNTG[df_MNTG$log_X >= 0,]

MNTG_m3 <- lmer(log_X ~ (Treatment + Year + Year2 + Gap + Fenced)^2 +
                 (1|Plot.Number) + (1|Species.Name), data=df_MNTG)
MNTG_m4 <- get_model(step(MNTG_m3))
summary(MNTG_m4)
plot(MNTG_m4)

RSS <- sum((df_MNTG$log_X - fitted(MNTG_m4, df_MNTG))^2)
TSS <- sum((df_MNTG$log_X - mean(df_MNTG$log_X))^2)
R2 <- 1-RSS/TSS
R2



#########################
# Final selected models #
#########################

# Large Forb / Herb
summary(LF_m4)

# Medium Tufted Grass
summary(MTG_m4)


# Medium Non-Tufted Grass
summary(MNTG_m4)




