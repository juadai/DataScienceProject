# CONTENTS

# Setup
# Diversity Modeling: Year 0 to 3: P(Y3>Y0)
# Diversity Modeling: Year 3 to 6: P(Y6>Y3)
# Diversity Modeling: Year 0 to 6: P(Y6>Y0)
# Diversity Modeling: Year 3: Y3 ~ N
# Diversity Modeling: Year 6: Y6 ~ N
# Explore and visualise


#########
# Setup #
#########

# Load data
library(rstudioapi)
library(pROC)

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

# Define function to diagnose binomial response
confusion <- function(model, y, data, threshold) {
  
  y_hat <- predict(model, data, type='response') > threshold
  conf <- matrix(0, nrow=2, ncol=2)
  rownames(conf) <- c('Predict 1', 'Predict 0')
  colnames(conf) <- c('Observe 1', 'Observe 0')
  
  conf[1,1] <- sum(y==1 & y_hat==1)
  conf[1,2] <- sum(y==0 & y_hat==1)
  conf[2,1] <- sum(y==1 & y_hat==0)
  conf[2,2] <- sum(y==0 & y_hat==0)
  
  tp_rate <- conf[1,1] / sum(conf[1,]) # Sensitivity
  tn_rate <- conf[2,2] / sum(conf[2,]) # Specificity
  
  return(list(Confusion=conf, Sensitivity=tp_rate, Specificity=tn_rate))
}

#############################################
# Diversity Modeling: Year 0 to 3: P(Y3>Y0) #
#############################################

# Model the probability that diversity increased between years 0 and 3
diversity_df$Y0.Y3.Increase <- diversity_df$X3 > diversity_df$X0
diversity_df$Y0.Y3.Increase <- as.numeric(diversity_df$Y0.Y3.Increase)

# Model 1: Treatment, Fenced, Gap, first-order interaction terms
diversity_Y0_Y3_01 <- glm(Y0.Y3.Increase ~ (Treatment + Fenced + Gap)^2,
                          family='binomial', data=diversity_df)
summary(diversity_Y0_Y3_01)
# Nothing significant yet

# Model 2: Step-wise selection fom Model 1 (AIC)
diversity_Y0_Y3_02 <- step(diversity_Y0_Y3_01)
summary(diversity_Y0_Y3_02)
# Treatment, Gap, and interaction are kept. But not significant (Wald test)

anova(diversity_Y0_Y3_01, diversity_Y0_Y3_02, test='Chi')
# No significant difference between the models


# Model 3: Check a model without interaction terms
diversity_Y0_Y3_03 <- glm(Y0.Y3.Increase ~ Treatment + Fenced + Gap,
                          family='binomial', data=diversity_df)
summary(diversity_Y0_Y3_03)
# Treatment, Fenced, Gap are in the model


# Model 4: Try step-wise selection from Model 1 (BIC)
k = log(length(diversity_df[,1]))
diversity_Y0_Y3_04 <- step(diversity_Y0_Y3_01, k=k)
summary(diversity_Y0_Y3_04)
# Only intercept is kept

anova(diversity_Y0_Y3_04, diversity_Y0_Y3_02, test='Chi')
# At 95% significance level, no indication that the larger model is better
# p-value 0.54, it's close

###########
# Results #
###########

# Investigate model 2
summary(diversity_Y0_Y3_02)

# Model 2: Intercept + Treatment + Gap + Treatment:Gap

# Estimated odds that diversity increased in a Gap Treatment quadrat (vs Control)
exp(-0.04382) # close to 1 ... Not significant

# Estimated odds that diversity increased in a Radial Treatment quadrat (vs Control)
exp(-0.14870) # 0.86

# Estimated 95% CI for these odds
c(exp(-0.14870-1.96*0.51500), exp(-0.14870+1.96*0.51500))
# [0.31, 2.36] ... Not significant


# Estimated odds that diversity increases in in-gap quadrats
exp(1.32721) # 3.77

# Estimated 95% CI for these odds
c(exp(1.32721-1.96*1.07676), exp(1.32721+1.96*1.07676))
# Extremely wide, not significant

r <- roc(diversity_df$Y0.Y3.Increase, predict(diversity_Y0_Y3_02, diversity_df, type='response'))
r$auc #0.60
plot.roc(r)

c <- confusion(diversity_Y0_Y3_02, diversity_df$Y0.Y3.Increase, diversity_df, 0.5)

c$Confusion # Matrix
c$Sensitivity # True Positive Rate
c$Specificity # True Negative Rate

# This model wouldn't generalise, i.e. not really effective for prediction



#############################################
# Diversity Modeling: Year 3 to 6: P(Y6>Y3) #
#############################################

# Model the probability that diversity increased between years 3 and 6
diversity_df$Y3.Y6.Increase <- diversity_df$X6 > diversity_df$X3
diversity_df$Y3.Y6.Increase <- as.numeric(diversity_df$Y3.Y6.Increase)

# Model 1: Treatment, Fenced, Gap, with first-order interaction terms
diversity_Y3_Y6_01 <- glm(Y3.Y6.Increase ~ (Treatment + Fenced + Gap)^2,
                          family='binomial', data=diversity_df)
summary(diversity_Y3_Y6_01)
# Nothing significant yet

# Model 2: Step-wise selection from Model 1 (AIC)
diversity_Y3_Y6_02 <- step(diversity_Y3_Y6_01)
summary(diversity_Y3_Y6_02)
# Treatment and Gap are significant

anova(diversity_Y3_Y6_02, diversity_Y3_Y6_01, test='Chi')
# No significant evidence that the larger model is a better fit

# Model 3: Step-wise selection from Model 1 (BIC)
k = log(length(diversity_df[,1]))
diversity_Y3_Y6_03 <- step(diversity_Y3_Y6_01)
summary(diversity_Y3_Y6_03)
# Same outcome

###########
# Results #
###########

# Not sufficient evidence that Fenced is significant

# Model 3: Intercept + Treatment + Gap
summary(diversity_Y3_Y6_03)

# Estimated odds that diversity increased in a Gap Treatment quadrat (vs Control)
exp(0.0332) # 1 ... not significant

# Estimated odds that diversity increased in a Radial Treatment quadrat (vs Control)
exp(1.1277) # 3.1

# Estimated 95% CI for these odds
c(exp(1.1277-1.96*0.4094), exp(1.1277+1.96*0.4094))
# [1.4, 6.9] ... Significant

# Estimated odds that diversity in an In Gap quadrat will increase (vs Not In Gap)
exp(-0.6685) # 0.51

# Estimated 95% CI for these odds
c(exp(-0.6685-1.96*0.3522), exp(-0.6685+1.96*0.3522))
# [0.26, 1.02] ... Maybe significant!


r <- roc(diversity_df$Y3.Y6.Increase, predict(diversity_Y3_Y6_03, diversity_df, type='response'))
r$auc # 0.64
plot.roc(r)

c <- confusion(diversity_Y3_Y6_03, diversity_df$Y3.Y6.Increase, diversity_df, 0.5)

c$Confusion # Matrix
c$Sensitivity # True Positive Rate
c$Specificity # True Negative Rate



#############################################
# Diversity Modeling: Year 0 to 6: P(Y6>Y0) #
#############################################

# Model the probability that diversity increased between years 0 and 6
diversity_df$Y0.Y6.Increase <- diversity_df$X6 > diversity_df$X0
diversity_df$Y0.Y6.Increase <- as.numeric(diversity_df$Y0.Y6.Increase)

# Model 1: Treatment, Fenced, Gap, and first-order interaction terms
diversity_Y0_Y6_01 <- glm(Y0.Y6.Increase ~ (Treatment + Fenced + Gap)^2,
                          family='binomial', data=diversity_df)
summary(diversity_Y0_Y6_01)
# Only intercept is signficant


# Model 2: Step-wise selection from Model 1 (AIC)
diversity_Y0_Y6_02 <- step(diversity_Y0_Y6_01)
summary(diversity_Y0_Y6_02)
# Treatment and Gap are in the model
anova(diversity_Y0_Y6_02, diversity_Y0_Y6_01, test='Chi')
# No significant difference between the two models


# Model 3: Step-wise selection from Model 1 (BIC)
k = log(length(diversity_df[,1]))
diversity_Y0_Y6_03 <- step(diversity_Y0_Y6_01, k=k)
summary(diversity_Y0_Y6_03)
# Only intercept remains

anova(diversity_Y0_Y6_03, diversity_Y0_Y6_02, test='Chi')
# Cannot take smaller model


###########
# Results #
###########

# Model 2: Intercept + Fenced + Gap + Intercept:Gap
summary(diversity_Y0_Y6_02)

# Estimated odds that diversity will increase in a gap treatment
exp(0.9014) # 2.4

# Estimated 95% CI for this interval
c(exp(0.9014-1.96*0.3345), exp(0.9014+1.96*0.3345))
# [1.3, 4.7] ... Significant


# Estimated odds that diversity will increase in a gap treatment, In Gap
exp(0.9014-0.7351) # 1.2

# Estimated odds that diversity will increase in a radial treatment
exp(0.8151) # 2.3

# Estimated 95% CI for this interval
c(exp(0.8151-1.96*0.3753), exp(0.8151+1.96*0.3753))
# [1.08, 4.71] ... Significant

# Estimated odds that diversity will increase in a radial treatment, In Gap
exp(0.8151-0.7351) # 1.08 ... Not significant


# Estimated odds that diversity will increase in a quadrat that's In Gap (Control)
exp(-0.7351) # 0.48

# Estimated 95% CI for this interval
c(exp(-0.7351-1.96*0.3026), exp(-0.7351+1.96*0.3026))
# [0.26, 0.87] ... Significant





r <- roc(diversity_df$Y0.Y6.Increase, predict(diversity_Y0_Y6_02, diversity_df, type='response'))
r$auc # 0.64
plot.roc(r)

c <- confusion(diversity_Y0_Y6_02, diversity_df$Y0.Y6.Increase, diversity_df, 0.5)

c$Confusion # Matrix
c$Sensitivity # True Positive Rate
c$Specificity # True Negative Rate

# ROC curve for Model 2 is a bit better than null model




###########################################
# Diversity Modeling: Year 3 to 6: Y3 ~ N #
###########################################

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

##################################
# Diversity Modeling: Year 6 ~ N #
##################################

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



###########
# Results #
###########

diversity_Y3_02$call
diversity_Y6_02$call

summary(diversity_Y3_02)
summary(diversity_Y6_02)




#######################
# Explore & Visualise #
#######################

cols1 <- c('springgreen3', 'orchid4', 'firebrick3')
cols2 <- c('mediumvioletred', 'midnightblue')

diversity_df$Col1 <- cols1[1]
diversity_df[diversity_df$Treatment == 'Gap', 'Col1'] <- cols1[2]
diversity_df[diversity_df$Treatment == 'Radial', 'Col1'] <- cols1[3]

diversity_df$Col2 <- cols2[1]
diversity_df[diversity_df$Fenced == FALSE, 'Col2'] <- cols2[2]
diversity_df$Col3 <- cols2[1]
diversity_df[diversity_df$Gap == FALSE, 'Col3'] <- cols2[2]

diversity_df_plot$Col1 <- cols1[1]
diversity_df_plot[diversity_df_plot$Treatment == 'Gap', 'Col1'] <- cols1[2]
diversity_df_plot[diversity_df_plot$Treatment == 'Radial', 'Col1'] <- cols1[3]

diversity_df_plot$Col2 <- cols2[1]
diversity_df_plot[diversity_df_plot$Fenced == FALSE, 'Col12'] <- cols1[2]
diversity_df_plot$Col3 <- cols2[1]
diversity_df_plot[diversity_df_plot$Gap == FALSE, 'Col3'] <- cols2[2]



line_plot <- function(df, label) {
  n <- length(df$Treatment)
  plot(0, type='l', xlim=c(0,6), ylim=c(-0.5,3.25),
       xlab='Year', ylab="Shannon's Diversity Index",
       main=paste(label, 'Plots'), col='white')
  for (i in 1:n) {
    lines(c(df[i,]$X0, df[i,]$X3, df[i,]$X6) ~ c(0,3,6),
          col=df[i,]$Col1)
  }
}

# All plots
line_plot(diversity_df_plot, 'All')
# Control plots
line_plot(diversity_df_plot[diversity_df_plot$Treatment=='Control',], 'Control')
# Gap plots
line_plot(diversity_df_plot[diversity_df_plot$Treatment=='Gap',], 'Gap')
# Radial plots
line_plot(diversity_df_plot[diversity_df_plot$Treatment=='Radial',], 'Radial')




#####################
# Plot by Treatment #
#####################

# Change in SDI Y0 to Y6 (plots)
plot(diversity_df_plot$X6-diversity_df_plot$X0,
     col=diversity_df_plot$Col1, pch=16, 
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y6")
lines(-25:300, rep(mean(diversity_df_plot$X6-diversity_df_plot$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in SDI Y0 to Y3 (plots)
plot(diversity_df_plot$X3-diversity_df_plot$X0,
     col=diversity_df_plot$Col1, pch=16, 
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y3")
lines(-25:300, rep(mean(diversity_df_plot$X3-diversity_df_plot$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in SDI Y3 to Y6 (plots)
plot(diversity_df_plot$X6-diversity_df_plot$X3,
     col=diversity_df_plot$Col1, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y3 to Y6")
lines(-25:300, rep(mean(diversity_df_plot$X6-diversity_df_plot$X3),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)




# Change in SDI Y0 to Y6 (quadrats)
plot(diversity_df$X6-diversity_df$X0,
     col=diversity_df$Col1, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y6")
lines(-25:300, rep(mean(diversity_df$X6-diversity_df$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in SDI Y0 to Y3
plot(diversity_df$X3-diversity_df$X0,
     col=diversity_df$Col1, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y3")
lines(-25:300, rep(mean(diversity_df$X3-diversity_df$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in SDI Y3 to Y6
plot(diversity_df$X6-diversity_df$X3,
     col=diversity_df$Col1, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y3 to Y6")
lines(-25:300, rep(mean(diversity_df$X6-diversity_df$X3),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)



# Magnitude of change in SDI Y0 to Y6
plot(abs(diversity_df$X6-diversity_df$X0),
     col=diversity_df$Col1, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change',
     main="Change in diversity: Y0 to Y6")
lines(-25:300, rep(mean(abs(diversity_df$X6-diversity_df$X0)),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Magnitude of change in SDI Y0 to Y3
plot(abs(diversity_df$X3-diversity_df$X0),
     col=diversity_df$Col1, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change',
     main="Change in diversity: Y0 to Y3")
lines(-25:300, rep(mean(abs(diversity_df$X3-diversity_df$X0)),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Magnitude of change in SDI Y3 to Y6
plot(abs(diversity_df$X6-diversity_df$X3),
     col=diversity_df$Col1, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change over time',
     main="Change in diversity: Y3 to Y6")
lines(-25:300, rep(mean(abs(diversity_df$X6-diversity_df$X3)),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)





##################
# Plot by Fenced #
##################

# Change in SDI Y0 to Y6 (quadrats)
plot(diversity_df$X6-diversity_df$X0,
     col=diversity_df$Col2, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y6")
lines(-25:300, rep(mean(diversity_df$X6-diversity_df$X0),326),
      col='red')
legend('topleft', c('Fenced', 'Unfenced'), pch=16, bty='n',
       col=cols2)

# Change in SDI Y0 to Y3
plot(diversity_df$X3-diversity_df$X0,
     col=diversity_df$Col2, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y3")
lines(-25:300, rep(mean(diversity_df$X3-diversity_df$X0),326),
      col='red')
legend('topleft', c('Fenced', 'Unfenced'), pch=16, bty='n',
       col=cols2)

# Change in SDI Y3 to Y6
plot(diversity_df$X6-diversity_df$X3,
     col=diversity_df$Col2, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y3 to Y6")
lines(-25:300, rep(mean(diversity_df$X6-diversity_df$X3),326),
      col='red')
legend('topleft', c('Fenced', 'Unfenced'), pch=16, bty='n',
       col=cols2)



# Magnitude of change in SDI Y0 to Y6
plot(abs(diversity_df$X6-diversity_df$X0),
     col=diversity_df$Col2, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change',
     main="Change in diversity: Y0 to Y6")
lines(-25:300, rep(mean(abs(diversity_df$X6-diversity_df$X0)),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in SDI  Y0 to Y3
plot(abs(diversity_df$X3-diversity_df$X0),
     col=diversity_df$Col2, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change',
     main="Change in diversity: Y0 to Y3")
lines(-25:300, rep(mean(abs(diversity_df$X3-diversity_df$X0)),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in SDI  Y3 to Y6
plot(abs(diversity_df$X6-diversity_df$X3),
     col=diversity_df$Col2, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change over time',
     main="Change in diversity: Y3 to Y6")
lines(-25:300, rep(mean(abs(diversity_df$X6-diversity_df$X3)),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)







##################
# Plot by In Gap #
##################

# Change in SDI Y0 to Y6 (quadrats)
plot(diversity_df$X6-diversity_df$X0,
     col=diversity_df$Col3, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y6")
lines(-25:300, rep(mean(diversity_df$X6-diversity_df$X0),326),
      col='red')
legend('topleft', c('In Gap', 'Not in Gap'), pch=16, bty='n',
       col=cols2)

# Change in SDI Y0 to Y3
plot(diversity_df$X3-diversity_df$X0,
     col=diversity_df$Col3, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y0 to Y3")
lines(-25:300, rep(mean(diversity_df$X3-diversity_df$X0),326),
      col='red')
legend('topleft', c('In Gap', 'Not in Gap'), pch=16, bty='n',
       col=cols2)

# Change in SDI Y3 to Y6
plot(diversity_df$X6-diversity_df$X3,
     col=diversity_df$Col3, pch=16,
     ylim=c(-3,3), ylab='Change over time',
     main="Change in diversity: Y3 to Y6")
lines(-25:300, rep(mean(diversity_df$X6-diversity_df$X3),326),
      col='red')
legend('topleft', c('In Gap', 'Not in Gap'), pch=16, bty='n',
       col=cols2)



# Magnitude of change in SDI Y0 to Y6
plot(abs(diversity_df$X6-diversity_df$X0),
     col=diversity_df$Col3, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change',
     main="Change in diversity: Y0 to Y6")
lines(-25:300, rep(mean(abs(diversity_df$X6-diversity_df$X0)),326),
      col='red')
legend('topleft', c('In Gap', 'Not in Gap'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in SDI  Y0 to Y3
plot(abs(diversity_df$X3-diversity_df$X0),
     col=diversity_df$Col3, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change',
     main="Change in diversity: Y0 to Y3")
lines(-25:300, rep(mean(abs(diversity_df$X3-diversity_df$X0)),326),
      col='red')
legend('topleft', c('In Gap', 'Not in Gap'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in SDI Y3 to Y6
plot(abs(diversity_df$X6-diversity_df$X3),
     col=diversity_df$Col3, pch=16,
     ylim=c(-0,3), ylab='Magnitude of change over time',
     main="Change in diversity: Y3 to Y6")
lines(-25:300, rep(mean(abs(diversity_df$X6-diversity_df$X3)),326),
      col='red')
legend('topleft', c('In Gap', 'Not in Gap'), pch=16, bty='n',
       col=cols2)


