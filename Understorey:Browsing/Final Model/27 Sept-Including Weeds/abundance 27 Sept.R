# Load data
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# Load data
abundance_df <- read.csv('../dataset/tables/lf_abund_analysis.csv')
abundance_matrix <- read.csv('../contingency tables/output/abundance/lf_abund_quadrats.csv')

# Total abundance for each quadrat
n = length(abundance_matrix[,1])
X0 <-rowSums(abundance_matrix[,abundance_matrix[1,] == 0])[3:n]
X3 <-rowSums(abundance_matrix[,abundance_matrix[1,] == 3])[3:n]
X6 <-rowSums(abundance_matrix[,abundance_matrix[1,] == 6])[3:n]

abundance_totals <- abundance_matrix[,1:5]
colnames(abundance_totals) <- abundance_totals[2,]
abundance_totals <- cbind(abundance_totals[3:n,], X0, X3, X6)

abundance_totals <- abundance_totals[abundance_totals$X0>0,]
abundance_totals <- abundance_totals[abundance_totals$X3>0,]
abundance_totals <- abundance_totals[abundance_totals$X6>0,]

colnames(abundance_totals) <- c('Treatment', 'Plot.Number', 'Quadrat.Number',
                                'Fenced', 'Gap', 'X0', 'X3', 'X6')

# Treatments, Plot numbers, Quadrat numbers, Life forms as factors
abundance_df$Treatment <- factor(abundance_df$Treatment)
abundance_df$Plot.Number <- factor(abundance_df$Plot.Number)
abundance_df$Quadrat.Number <- factor(abundance_df$Quadrat.Number)
abundance_df$Life.Form <- factor(abundance_df$Life.Form)

abundance_totals$Treatment <- factor(abundance_totals$Treatment)
abundance_totals$Plot.Number <- factor(abundance_totals$Plot.Number)
abundance_totals$Quadrat.Number <- factor(abundance_totals$Quadrat.Number)

# Convert Fenced, Gap columns to boolean
abundance_df$Gap <- abundance_df$Gap=='True'
abundance_df$Fenced <- abundance_df$Fenced=='True'

abundance_totals$Gap <- abundance_totals$Gap=='True'
abundance_totals$Fenced <- abundance_totals$Fenced=='True'




################################
# Total Abundance: log(Y3) ~ N #
################################

abundance_totals$log_X0 <- log(abundance_totals$X0)
abundance_totals$log_X3 <- log(abundance_totals$X3)
abundance_totals$log_X6 <- log(abundance_totals$X6)

abundance_Y3_11 <- glm(log_X3 ~ (log_X0 + Treatment + Gap + Fenced)^2, data=abundance_totals)
abundance_Y3_12 <- step(abundance_Y3_11)
summary(abundance_Y3_12)

k <- log(length(abundance_totals[,1]))
abundance_Y3_13 <- step(abundance_Y3_11, k=k)
summary(abundance_Y3_13)
anova(abundance_Y3_13, abundance_Y3_12, test='Chi')
# Cannot take the smaller model

summary(abundance_Y3_12)
plot(abundance_Y3_12) # Some outliers to check

################################
# Total Abundance: log(Y6) ~ N #
################################

abundance_Y6_11 <- glm(log_X6 ~ (log_X0 + log_X3 + Treatment + Gap + Fenced)^2, 
                       data=abundance_totals)
abundance_Y6_12 <- step(abundance_Y6_11)
summary(abundance_Y6_12)

abundance_Y6_13 <- step(abundance_Y6_11, k=k)
summary(abundance_Y6_13)

anova(abundance_Y6_13, abundance_Y6_12, test='Chi')
# Fit of models are similar at 5% significance level

plot(abundance_Y6_13)

mean(abundance_totals$X6)




##########################
# Model Medium Forb/Herb #
##########################
#EVC benchmark: #spp: 10 (richness), %cover: 20 (abundance)

# Setup
library(lme4)
library(lmerTest)

df <- read.csv('../dataset/tables/abundance_species_for_analysis.csv')
unique(df['Life.Form'])
colnames(df)

grouped <- aggregate(cbind(df$X0,df$X3,df$X6),
                     list(df$Life.Form), FUN=sum)
grouped <- grouped[order(grouped$V1, decreasing=T),]
rownames(grouped) <- grouped$Group.1

par(mar=c(15,5,3,3))
plot(grouped$V1, cex=0.75, ylim=c(0,6000), xaxt='n', xlab='',
     ylab='Sum of abundance across trial', main='Total abundance of each life form')
axis(1, at = 1:19, labels = grouped$Group.1,las=2)
points(1:19, grouped$V2, col='red', cex=0.75)
points(1:19, grouped$V3, col='blue', cex=0.75)
legend('topright',legend=c('Year 0', 'Year 3', 'Year 6'),
       pch=1,
       col=c('black', 'red', 'blue'),
       cex=0.75)

df$Life.Form
df <- df[df$Life.Form=='Medium Forb/Herb',]

df$Plot <- factor(df$Plot)
df$Treatment <- factor(df$Treatment)
df$Fenced <- df$Fenced=="True"
df$Gap <- df$Gap=="True"


# Visualisations
grouped <- aggregate(cbind(df$X0,df$X3,df$X6), 
                     by=list(df$Treatment,df$Plot.Number,df$Quadrat.Number,df$Fenced,df$Gap,df$Life.Form), 
                     FUN=sum)

hist(grouped$V1, ylim=c(0,230), xlim=c(0,65),
     main='Medium Forb/Herb - Year 0', 
     xlab='Abundance (BB)', ylab='Number of quadrats')
abline(v=20, col='red')
hist(grouped$V2, ylim=c(0,230), xlim=c(0,65),
     main='Medium Forb/Herb - Year 0', 
     xlab='Abundance (BB)', ylab='Number of quadrats')
abline(v=20, col='red')
hist(grouped$V3, ylim=c(0,230), xlim=c(0,65),
     main='Medium Forb/Herb - Year 0', 
     xlab='Abundance (BB)', ylab='Number of quadrats')
abline(v=20, col='red')



# Modeling Y3
mixed_Y3_01 <- lmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 +
                         (1|Plot.Number) + (1 + Treatment|Species.Name),
                         data=df)
summary(mixed_Y3_01)
plot(mixed_Y3_01)

mixed_Y3_03 <- step(mixed_Y3_01)
mixed_Y3_03 <- get_model(mixed_Y3_03)
summary(mixed_Y3_03)
# Drops (Treatment|Species)

mixed_Y3_03

ranova(mixed_Y3_01)
# Other random effects highly significant
# (1|Plot Number), (1|Species.Name)

qqnorm(residuals(mixed_Y3_01))
qqline(residuals(mixed_Y3_01))

mixed_Y3_02 <- glmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 +
                          (1|Plot.Number) + (1 + Treatment + Gap | Species.Name),
                          family='poisson', data=df)
summary(mixed_Y3_02)
plot(mixed_Y3_02)


# Modeling Y6
mixed_Y6_01 <- lmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 +
                         (1|Plot.Number) + (1 + Treatment + Gap | Species.Name),
                         data=df)
summary(mixed_Y6_01)
plot(mixed_Y6_01)
qqnorm(residuals(mixed_Y6_01))
qqline(residuals(mixed_Y6_01))




############################################
# Model Prostrate or Mat-forming Forb/Herb #
############################################
#EVC benchmark: #spp: 5 (richness), %cover: 10 (abundance)


# Setup
PMFH_df <- read.csv('../dataset/tables/abundance_species_for_analysis.csv')

PMFH_df <- PMFH_df[PMFH_df$Life.Form=='Prostrate or Mat-forming Forb/Herb',]

PMFH_df$Plot <- factor(PMFH_df$Plot)
PMFH_df$Treatment <- factor(PMFH_df$Treatment)
PMFH_df$Fenced <- PMFH_df$Fenced=="True"
PMFH_df$Gap <- PMFH_df$Gap=="True"


# Visualisations
grouped <- aggregate(cbind(PMFH_df$X0,PMFH_df$X3,PMFH_df$X6), 
                     by=list(PMFH_df$Treatment,PMFH_df$Plot.Number,PMFH_df$Quadrat.Number,PMFH_df$Fenced,PMFH_df$Gap,PMFH_df$Life.Form), 
                     FUN=sum)

hist(grouped$V1, ylim=c(0,250), xlim=c(0,80),
     main='Prostrate or Mat-forming Forb/Herb - Year 0', 
     xlab='Abundance (BB)', ylab='Number of quadrats', nclass=20)
abline(v=10, col='red')
hist(grouped$V2, ylim=c(0,250), xlim=c(0,80),
     main='Prostrate or Mat-forming Forb/Herb - Year 0', 
     xlab='Abundance (BB)', ylab='Number of quadrats', nclass=20)
abline(v=10, col='red')
hist(grouped$V3, ylim=c(0,250), xlim=c(0,80),
     main='Prostrate or Mat-forming Forb/Herb - Year 0', 
     xlab='Abundance (BB)', ylab='Number of quadrats', nclass=10)
abline(v=10, col='red')


?hist
# Modeling Y3

# Fixed-effect models
PMFH_Y3_01 <- lm(X3 ~ (X0 + Treatment + Gap + Fenced)^2, data=PMFH_df)

PMFH_Y3_02 <- step(PMFH_Y3_01)
summary(PMFH_Y3_02)

anova(PMFH_Y3_02, PMFH_Y3_01) #choose smaller model


# Mixed-effect models
PMFH_mixed_Y3_01 <- lmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 + (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y3_01) #random effect of plot.number on intercept is not significant at all
coef(PMFH_mixed_Y3_01)
plot(PMFH_mixed_Y3_01)
AIC(PMFH_mixed_Y3_01)

qqnorm(residuals(PMFH_mixed_Y3_01))
qqline(residuals(PMFH_mixed_Y3_01))


PMFH_mixed_Y3_02 <- lmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 + (1 + Treatment + Gap|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y3_02)
coef(PMFH_mixed_Y3_02)
plot(PMFH_mixed_Y3_02)
AIC(PMFH_mixed_Y3_02) # lower than model 1

qqnorm(residuals(PMFH_mixed_Y3_02))
qqline(residuals(PMFH_mixed_Y3_02))


PMFH_mixed_Y3_03 <- glmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 +  (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y3_03)
plot(PMFH_mixed_Y3_03)
AIC(PMFH_mixed_Y3_03) #lower than linear mixed model


PMFH_mixed_Y3_04 <- glmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 +  (1 + Treatment + Gap | Plot.Number) + (1 + Treatment + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y3_04)
plot(PMFH_mixed_Y3_04) #lower residuals than model 3
AIC(PMFH_mixed_Y3_04)  #much lower than model 3


PMFH_mixed_Y3_05 <- lmer(X3 ~ X0 + Treatment + Gap + Fenced + X0:Treatment + X0:Fenced + Treatment:Gap 
                         + (1 + Treatment + Gap | Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y3_05)
plot(PMFH_mixed_Y3_05)
AIC(PMFH_mixed_Y3_05)

qqnorm(residuals(PMFH_mixed_Y3_04))
qqline(residuals(PMFH_mixed_Y3_04))


PMFH_mixed_Y3_06 <- glmer(X3 ~ X0 + Treatment + Gap + Fenced + X0:Treatment + X0:Fenced + Treatment:Gap 
                         + (1 + Treatment + Gap | Plot.Number) + (1 + Treatment + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y3_06)
plot(PMFH_mixed_Y3_06)
AIC(PMFH_mixed_Y3_06) #lower than model 4


PMFH_mixed_Y3_07 <- glmer(X3 ~ X0 + Treatment + Gap + Fenced + X0:Treatment + X0:Fenced 
                          + (1 + Treatment + Gap | Plot.Number) + (1 + Treatment + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y3_07)
plot(PMFH_mixed_Y3_07)
AIC(PMFH_mixed_Y3_07)


PMFH_mixed_Y3_08 <- glmer(X3 ~ X0 + Treatment + Gap + Fenced + X0:Fenced 
                          + (1 + Treatment + Gap | Plot.Number) + (1 + Treatment + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y3_08)
plot(PMFH_mixed_Y3_08)
AIC(PMFH_mixed_Y3_08)

anova(PMFH_mixed_Y3_08, PMFH_mixed_Y3_07)

PMFH_mixed_Y3_09 <- glmer(X3 ~ X0 + Gap + Fenced + X0:Fenced 
                          + (1 + Gap | Plot.Number) + (1 + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y3_09)
plot(PMFH_mixed_Y3_09)
AIC(PMFH_mixed_Y3_09)

#Mixed-effect model 8 is the best here 


# Modeling Y6

# Fixed-effect models
PMFH_Y6_01 <- lm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, data=PMFH_df)

PMFH_Y6_02 <- step(PMFH_Y6_01)
summary(PMFH_Y6_02)

anova(PMFH_Y6_02, PMFH_Y6_01) #choose smaller model


# Mixed-effect models
PMFH_mixed_Y6_01 <- lmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 + (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y6_01)
plot(PMFH_mixed_Y6_01)
AIC(PMFH_mixed_Y6_01)

qqnorm(residuals(PMFH_mixed_Y6_01))
qqline(residuals(PMFH_mixed_Y6_01)) #weird looking qq plot at the right tail

anova(PMFH_mixed_Y6_01)


PMFH_mixed_Y6_02 <- lmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 + (1 + Treatment + Gap|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y6_02)
plot(PMFH_mixed_Y6_02)
AIC(PMFH_mixed_Y6_02) # worse than model 1

qqnorm(residuals(PMFH_mixed_Y6_02))
qqline(residuals(PMFH_mixed_Y6_02))


PMFH_mixed_Y6_03 <- glmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 +  (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y6_03)
plot(PMFH_mixed_Y6_03) #larger residuals than model 1
AIC(PMFH_mixed_Y6_03) # much lower than model 1


PMFH_mixed_Y6_04 <- lmer(X6 ~ (X0 + X3 + Treatment + Fenced)^2 + (1|Plot.Number) + (1 + Treatment | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y6_04)
plot(PMFH_mixed_Y6_04)
AIC(PMFH_mixed_Y6_04) #lower than model 1 but nothing seems significant


PMFH_mixed_Y6_05 <- lmer(X6 ~ (X0 + X3 + Treatment + Fenced)^2 - X0:Fenced + (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y6_05)
plot(PMFH_mixed_Y6_05)
AIC(PMFH_mixed_Y6_05) #lower than model 4 but nothing seems significant


PMFH_mixed_Y6_06 <- lmer(X6 ~ (X3 + Treatment + Fenced)^2 + X0 + (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y6_06)
plot(PMFH_mixed_Y6_06)
AIC(PMFH_mixed_Y6_06) #lower than model 5 but nothing seems significant


PMFH_mixed_Y6_07 <- glmer(X6 ~ (X0 + X3 + Treatment + Fenced)^2 - X0:Fenced - X3:Fenced - X0:Treatment +  (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), family='poisson', data=PMFH_df)
summary(PMFH_mixed_Y6_07)
plot(PMFH_mixed_Y6_07)
AIC(PMFH_mixed_Y6_07)

#Try negative binomial model
PMFH_mixed_Y6_08 <- glmer.nb(X6 ~ (X0 + X3 + Treatment + Fenced)^2 - X0:Fenced - X3:Fenced - X0:Treatment +  (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=PMFH_df)
summary(PMFH_mixed_Y6_08)
plot(PMFH_mixed_Y6_08) #larger residuals than model 1
AIC(PMFH_mixed_Y6_08)
# Took too long to run 

#choose model 7 here (but the residuals are very high compared to the fitted values)




##################################
# Model Small Tufted grass/sedge #
##################################
#EVC benchmark: #spp: 9 (richness), % cover: 35 (abundance)


# Setup
STGS_df <- read.csv('../dataset/tables/abundance_species_for_analysis.csv')

STGS_df <- STGS_df[STGS_df$Life.Form=='Small Tufted grass/sedge',]

STGS_df$Plot <- factor(STGS_df$Plot)
STGS_df$Treatment <- factor(STGS_df$Treatment)
STGS_df$Fenced <- STGS_df$Fenced=="True"
STGS_df$Gap <- STGS_df$Gap=="True"


# Modeling Y3

# Fixed-effect models
STGS_Y3_01 <- lm(X3 ~ (X0 + Treatment + Gap + Fenced)^2, data=STGS_df)

STGS_Y3_02 <- step(PMFH_Y3_01)
summary(STGS_Y3_02)

anova(STGS_Y3_02, PMFH_Y3_01) #choose smaller model


# Mixed-effect models
STGS_mixed_Y3_01 <- lmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 + (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=STGS_df)
summary(STGS_mixed_Y3_01) 
plot(STGS_mixed_Y3_01)
AIC(STGS_mixed_Y3_01)

STGS_mixed_Y3_02 <- glmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 + (1|Plot.Number) + (1 + Treatment + Gap | Species.Name),family='poisson', data=STGS_df)
summary(STGS_mixed_Y3_02) 
plot(STGS_mixed_Y3_02) # yield lower residuals than gaussian GLM
AIC(STGS_mixed_Y3_02) #much lower than model 1

STGS_mixed_Y3_03 <- glmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 + (1 + Treatment + Gap |Plot.Number) + (1 + Treatment + Gap | Species.Name),family='poisson', data=STGS_df)
summary(STGS_mixed_Y3_03) 
plot(STGS_mixed_Y3_03) # yield a bit lower residuals than model 2
AIC(STGS_mixed_Y3_03)

STGS_mixed_Y3_04 <- glmer(X3 ~ (X0 + Treatment + Fenced)^2 + (1 + Treatment + Gap |Plot.Number) + (1 + Treatment + Gap | Species.Name),family='poisson', data=STGS_df)
summary(STGS_mixed_Y3_04) 
plot(STGS_mixed_Y3_04) 
AIC(STGS_mixed_Y3_04)

#choose model 4




# Modeling Y6

# Fixed-effect models
STGS_Y6_01 <- lm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, data=STGS_df)

STGS_Y6_02 <- step(STGS_Y6_01)
summary(STGS_Y6_02)

anova(STGS_Y6_02, STGS_Y6_01) #choose smaller model


# Mixed-effect models
STGS_mixed_Y6_01 <- lmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 + (1|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=STGS_df)
summary(STGS_mixed_Y6_01) # Plot number has no random effect on intercept
plot(STGS_mixed_Y6_01) # seems to have 2 outliers 
AIC(STGS_mixed_Y6_01)

qqnorm(residuals(STGS_mixed_Y6_01))
qqline(residuals(STGS_mixed_Y6_01)) #weird looking qq plot at the right tail

STGS_mixed_Y6_02 <- lmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 + (1 + Treatment + Gap|Plot.Number) + (1 + Treatment + Gap | Species.Name), data=STGS_df)
summary(STGS_mixed_Y6_02) #Random effects from Plot.Number are really small
plot(STGS_mixed_Y6_02)
AIC(STGS_mixed_Y6_02) # worse than model 1

qqnorm(residuals(STGS_mixed_Y6_02))
qqline(residuals(STGS_mixed_Y6_02))

STGS_mixed_Y6_03 <- glmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 + (1 + Treatment + Gap|Plot.Number) + (1 + Treatment + Gap | Species.Name),family='poisson', data=STGS_df)
summary(STGS_mixed_Y6_03) 
plot(STGS_mixed_Y6_03) # much lower residuals than the first 2 models
AIC(STGS_mixed_Y6_03) # much lower than the first 2 models

STGS_mixed_Y6_04 <- glmer(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2 - X3:Treatment + (1 + Treatment + Gap|Plot.Number) + (1 + Treatment + Gap | Species.Name),family='poisson', data=STGS_df)
summary(STGS_mixed_Y6_04) 
plot(STGS_mixed_Y6_04) 
AIC(STGS_mixed_Y6_04)

STGS_mixed_Y6_05 <- glmer(X6 ~  (X0 + X3 + Treatment + Gap + Fenced)^2 - X3:Treatment  - X0:Treatment + (1 + Treatment + Gap|Plot.Number) + (1 + Treatment + Gap | Species.Name),family='poisson', data=STGS_df)
summary(STGS_mixed_Y6_05) 
plot(STGS_mixed_Y6_05) 
AIC(STGS_mixed_Y6_05)

#model 5 is the best here, but the residuals are quite high compared to the fitted values
