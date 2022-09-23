# CONTENTS

# Setup
# Total Abundance Modeling: Year 3 ~ Normal
# Total Abundance Modeling: Year 6 ~ Normal
# Total Abundance Modeling: log(Year 3) ~ Normal
# Total Abundance Modeling: log(Year 6) ~ Normal
# Total Abundance Modeling: Year 3 ~ Pois
# Total Abundance Modeling: Year 3 ~ Pois
# Explore and visualise

#########
# Setup #
#########

# Load data
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# Load data
abundance_df <- read.csv('../dataset/tables/lf_abund_analysis.csv')
abundance_rel_df <- read.csv('../dataset/tables/lf_abund_rel_analysis.csv')
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

abundance_rel_df$Treatment <- factor(abundance_rel_df$Treatment)
abundance_rel_df$Plot.Number <- factor(abundance_rel_df$Plot.Number)
abundance_rel_df$Quadrat.Number <- factor(abundance_rel_df$Quadrat.Number)
abundance_rel_df$Life.Form <- factor(abundance_rel_df$Life.Form)

abundance_totals$Treatment <- factor(abundance_totals$Treatment)
abundance_totals$Plot.Number <- factor(abundance_totals$Plot.Number)
abundance_totals$Quadrat.Number <- factor(abundance_totals$Quadrat.Number)

# Convert Fenced, Gap columns to boolean
abundance_df$Gap <- abundance_df$Gap=='True'
abundance_df$Fenced <- abundance_df$Fenced=='True'

abundance_rel_df$Gap <- abundance_rel_df$Gap=='True'
abundance_rel_df$Fenced <- abundance_rel_df$Fenced=='True'

abundance_totals$Gap <- abundance_totals$Gap=='True'
abundance_totals$Fenced <- abundance_totals$Fenced=='True'

# Relative Abundance: convert from proportion to percentage, and round to integer
# Note: Trace observations (below 1%) rounded up to 1%
abundance_rel_df$response <- abundance_rel_df$X0*100
abundance_rel_df[abundance_rel_df$response < 1 & abundance_rel_df$response != 0, 'response'] <- 1
abundance_rel_df$response <- round(abundance_rel_df$response)




###########################
# Total Abundance: Y3 ~ N #
###########################

str(abundance_totals)
k <- log(length(abundance_totals[,1]))

abundance_Y3_01 <- glm(X3 ~ (X0 + Treatment + Gap + Fenced)^2, data=abundance_totals)
abundance_Y3_02 <- step(abundance_Y3_01)
summary(abundance_Y3_02)

anova(abundance_Y3_02, abundance_Y3_01, test='Chi')

abundance_Y3_03 <- step(abundance_Y3_01, k=k)
summary(abundance_Y3_03)

anova(abundance_Y3_03, abundance_Y3_02, test='Chi')
# Keep model 2 - Fenced and its interactions are significant


###########################
# Total Abundance: Y6 ~ N #
###########################

abundance_Y6_01 <- glm(X6 ~ (X0 + X3 + Treatment + Gap + Fenced)^2, data=abundance_totals)
abundance_Y6_02 <- step(abundance_Y6_01)
summary(abundance_Y6_02)

anova(abundance_Y6_02, abundance_Y6_01, test='Chi')
# No significant difference between the models

abundance_Y6_03 <- step(abundance_Y6_01, k=k)
summary(abundance_Y6_03)

anova(abundance_Y6_03, abundance_Y6_02, test='Chi')
# Cannot take the smaller model

plot(abundance_Y6_02)

################################
# Total Abundance: log(Y3) ~ N #
################################

abundance_totals$log_X0 <- log(abundance_totals$X0)
abundance_totals$log_X3 <- log(abundance_totals$X3)
abundance_totals$log_X6 <- log(abundance_totals$X6)

abundance_Y3_11 <- glm(log_X3 ~ (log_X0 + Treatment + Gap + Fenced)^2, data=abundance_totals)
abundance_Y3_12 <- step(abundance_Y3_11)
summary(abundance_Y3_12)

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

##############################
# Total Abundance: Y3 ~ Pois #
##############################

str(abundance_totals)
k <- log(length(abundance_totals[,1]))

abundance_totals$X0.Int <- round(abundance_totals$X0)
abundance_totals$X3.Int <- round(abundance_totals$X3)
abundance_totals$X6.Int <- round(abundance_totals$X6)

abundance_Y3_21 <- glm(X3.Int ~ (X0.Int + Treatment + Gap + Fenced)^2, 
                       family='poisson', data=abundance_totals)
abundance_Y3_22 <- step(abundance_Y3_21)
summary(abundance_Y3_21)

##############################
# Total Abundance: Y6 ~ Pois #
##############################

abundance_Y6_21 <- glm(X6.Int ~ (X0.Int + X3.Int + Treatment + Gap + Fenced)^2, data=abundance_totals)
abundance_Y6_22 <- step(abundance_Y6_21)
summary(abundance_Y6_22)

anova(abundance_Y6_22, abundance_Y6_21, test='Chi')
# No significant difference between the models

abundance_Y6_23 <- step(abundance_Y6_21, k=k)
summary(abundance_Y6_23)

anova(abundance_Y6_23, abundance_Y6_22, test='Chi')
# Cannot take the smaller model

summary(abundance_Y6_22)
plot(abundance_Y6_22) # Not looking good

#######################
# Explore & Visualise #
#######################

plot(abundance_totals$X0, ylab='Total Quadrat Abundance', main='Year 0 (Midpoint)')
lines(c(-50,300), rep(100,2), col='red')
plot(abundance_totals$X3, ylab='Total Quadrat Abundance', main='Year 3 (Midpoint)')
lines(c(-50,300), rep(100,2), col='red')
plot(abundance_totals$X6, ylab='Total Quadrat Abundance', main='Year 6 (Midpoint)')
lines(c(-50,300), rep(100,2), col='red')


cols1 <- c('springgreen3', 'orchid4', 'firebrick3')
cols2 <- c('mediumvioletred', 'midnightblue')

abundance_totals$Col1 <- cols1[1]
abundance_totals[abundance_totals$Treatment == 'Gap', 'Col1'] <- cols1[2]
abundance_totals[abundance_totals$Treatment == 'Radial', 'Col1'] <- cols1[3]

abundance_totals$Col2 <- cols2[1]
abundance_totals[abundance_totals$Fenced == FALSE, 'Col2'] <- cols2[2]
abundance_totals$Col3 <- cols2[1]
abundance_totals[abundance_totals$Gap == FALSE, 'Col3'] <- cols2[2]



line_plot <- function(df, label) {
  n <- length(df$Treatment)
  plot(0, type='l', xlim=c(0,6), ylim=c(-0.5,180),
       xlab='Year', ylab="Total abundance (%)",
       main=paste(label, 'Quadrats'), col='white')
  for (i in 1:n) {
    lines(c(df[i,]$X0, df[i,]$X3, df[i,]$X6) ~ c(0,3,6),
          col=df[i,]$Col1)
  }
}

# All quadrats
line_plot(abundance_totals, 'All')
# Control quadrats
line_plot(abundance_totals[abundance_totals$Treatment=='Control',], 'Control')
# Gap quadrats
line_plot(abundance_totals[abundance_totals$Treatment=='Gap',], 'Gap')
# Radial quadrats
line_plot(abundance_totals[abundance_totals$Treatment=='Radial',], 'Radial')


# In-Gap quadrats
line_plot(abundance_totals[abundance_totals$Gap=='True',], 'In-Gap')
# Not-In-Gap quadrats
line_plot(abundance_totals[abundance_totals$Gap=='False',], 'Not In-Gap')

# Fenced quadrats
line_plot(abundance_totals[abundance_totals$Fenced=='True',], 'Fenced')
# Unfenced quadrats
line_plot(abundance_totals[abundance_totals$Fenced=='False',], 'Unfenced')




