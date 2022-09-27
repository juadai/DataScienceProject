
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

df <- read.csv('../dataset/tables/abundance_species_for_analysis.csv')
unique(df['Life.Form'])
colnames(df)

df <- df[df$Life.Form=='Medium Forb/Herb',]

df$Plot <- factor(df$Plot)
df$Treatment <- factor(df$Treatment)
df$Fenced <- df$Fenced=="True"
df$Gap <- df$Gap=="True"


# Modeling Y3
mixed_Y3_01 <- lmer(X3 ~ (X0 + Treatment + Gap + Fenced)^2 +
                         (1|Plot.Number) + (1 + Treatment + Gap | Species.Name),
                         data=df)
summary(mixed_Y3_01)
plot(mixed_Y3_01)

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

