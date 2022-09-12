#########
# Setup #
#########

# Load data
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
richness_df <- read.csv('../dataset/tables/species_richness_quadrats.csv')
richness_df_plot <- read.csv('../contingency tables/output/richness/species_richness_plots.csv')

# Treatments, Plot numbers, Quadrat numbers, Life forms as factors
richness_df$Treatment <- factor(richness_df$Treatment)
richness_df$Plot.Number <- factor(richness_df$Plot.Number)
richness_df$Quadrat.Number <- factor(richness_df$Quadrat.Number)

# Convert Fenced, Gap columns to boolean
richness_df$Gap <- richness_df$Gap=='True'
richness_df$Fenced <- richness_df$Fenced=='True'

richness_df_plot$Gap <- richness_df_plot$Gap=='True'
richness_df_plot$Fenced <- richness_df_plot$Fenced=='True'

# Remove null rows
richness_df <- na.omit(richness_df)
richness_df_plot <- na.omit(richness_df_plot)

#######################################
# Richness Modeling: Year 0 to Year 3 #
#######################################

# Model the probability that richness increased between years 0 and 3
richness_df$Y0.Y3.Increase <- richness_df$X3 > richness_df$X0
richness_df$Y0.Y3.Increase <- as.numeric(richness_df$Y0.Y3.Increase)

# Model 1: Treatment, Fenced, Gap, all first order interaction terms
richness_Y0_Y3_01 <- glm(Y0.Y3.Increase ~ (Treatment + Fenced + Gap)^2,
                         family='binomial', data=richness_df)
summary(richness_Y0_Y3_01)
# In Gap is significant, as is Gap Treatment:In Gap interaction

# Model 2: Step-wise selection from Model 1 (AIC)
richness_Y0_Y3_02 <- step(richness_Y0_Y3_01)
summary(richness_Y0_Y3_02)
# Removed interaction between treatment and fenced

anova(richness_Y0_Y3_01, richness_Y0_Y3_02, test='Chi')
# Insufficient evidence to suggest that model 1 is better

# Model 3: Step-wise selection from Model 1 (BIC)
k = log(length(richness_df[,1]))
richness_Y0_Y3_03 <- step(richness_Y0_Y3_01, k=k)
summary(richness_Y0_Y3_03)
# Removed pretty much everything

# Compare BIC-selected model with full model
anova(richness_Y0_Y3_01, richness_Y0_Y3_03, test='Chi')
# Cannot reject the hypothesis that there's no difference in fit
# Keep larger model 1 over very simple model 3

# Compare BIC- and AIC-selected models
anova(richness_Y0_Y3_02, richness_Y0_Y3_03, test='Chi')
# Cannot reject the hypothesis that there's no difference in fit
# Keep model 2 (selected by AIC) over model 3 (selected by BIC)



# Model 4: Now without interaction
richness_Y0_Y3_04 <- glm(Y0.Y3.Increase ~ Treatment + Fenced + Gap,
                         family='binomial', data=richness_df)
summary(richness_Y0_Y3_04)
# Radial treatment and In-Gap are now significant


anova(richness_Y0_Y3_04, richness_Y0_Y3_02, test='Chi')
# At 95% significance level, can reject hypothesis that interaction is 
# important. But just barely (p-value 0.055)

# Model 5: Step-wise selection from Model 4 (AIC)
richness_Y0_Y3_05 <- step(richness_Y0_Y3_04)
summary(richness_Y0_Y3_05)
# Radial is significant, as is Gap=True


# Model 6: Try treatment, gap, and first order interactions (exclude fenced)
richness_Y0_Y3_06 <- glm(Y0.Y3.Increase ~ (Treatment + Gap)^2,
                         family='binomial', data=richness_df)
summary(richness_Y0_Y3_06)

# Model 7: Step-wise selection from Model 6 (AIC)
richness_Y0_Y3_07 <- step(richness_Y0_Y3_06)
# No change

# Test for significance of interaction, without the presence of Fenced
anova(richness_Y0_Y3_05, richness_Y0_Y3_06, test='Chi')
# At 95% significance level, may reject hypothesis that larger 
# model (with interaction) is better (p-value 0.072)

# Test the two best candidate models
anova(richness_Y0_Y3_02, richness_Y0_Y3_05, test='Chi')
# May reject the hypothesis that larger model has a better fit
# So we take the simpler model (model 5)

###########
# Results #
###########

# Model 3: Increase ~ Intercept + Treatment + In Gap
summary(richness_Y0_Y3_05)
# Model excludes interaction, and was selected through AIC

# However - should further investigate the impact of interaction


# Estimated odds of an In-Gap quadrat increasing in richness,
# between years 0 and 3 (vs Not-InGap quadrat)
exp(0.71121) # 2

# Estimated 95% CI for these odds:
c(exp(0.71121-1.96*0.29798), exp(0.71121+1.96*0.29798))
# [1.14, 3.65]

# Estimated odds of a Radial Treatment quadrat increasing in richness,
# between years 0 and 3 (vs Control quadrat)
exp(-0.84195) # 0.43

# Estimated 95% CI for these odds:
c(exp(-0.84195-1.96*0.36442), exp(-0.84195+1.96*0.36442))
# [0.21, 0.88]


# Estimated odds of a Gap Treatment quadrat increasing in richness,
# between years 0 and 3 (vs Control quadrat)
exp(0.10991) # 1.12

# Estimated 95% CI for these odds:
c(exp(0.10991-1.96*0.31672), exp(0.10991+1.96*0.31672))
# [0.6, 2.1] ... No significant difference Gap vs Control



##################################
# Richness Modeling: Year 3 to 6 #
##################################


# Model the probability that richness increased between years 3 and 6
richness_df$Y3.Y6.Increase <- richness_df$X6 > richness_df$X3
richness_df$Y3.Y6.Increase <- as.numeric(richness_df$Y3.Y6.Increase)

# Model 1: Treatment, Fenced, Gap, all first order interaction terms
richness_Y3_Y6_01 <- glm(Y3.Y6.Increase ~ (Treatment + Fenced + Gap)^2,
                         family='binomial', data=richness_df)
summary(richness_Y3_Y6_01)
# Only intercept is significant

# Model 2: Step-wise selection from Model 1 (AIC)
richness_Y3_Y6_02 <- step(richness_Y3_Y6_01)
summary(richness_Y3_Y6_02)

# Now fenced also a significant factor
# Some interaction, but the std errors are large and terms are not significant

# Model 3: Step-wise selection from Model 1 (BIC)
k = log(length(richness_df[,1]))
richness_Y3_Y6_03 <- step(richness_Y3_Y6_01, k=k)
summary(richness_Y3_Y6_03)

# Fenced is the only significant factor

anova(richness_Y3_Y6_03, richness_Y3_Y6_01, test='Chi')
# At 95% significant level, we may reject hypothesis that model 1 is better



# Model 4: Now without interaction
richness_Y3_Y6_04 <- glm(Y3.Y6.Increase ~ Treatment + Fenced + Gap,
                         family='binomial', data=richness_df)
summary(richness_Y3_Y6_04)

anova(richness_Y3_Y6_04, richness_Y3_Y6_01, test='Chi')
# Suggests interactions are important. But their std. errors are messed up

# Try a model with fencing and gap (exclude treatments)
richness_Y3_Y6_05 <- glm(Y3.Y6.Increase ~ Fenced + Gap,
                         family='binomial', data=richness_df)
summary(richness_Y3_Y6_05)

anova(richness_Y3_Y6_05, richness_Y3_Y6_04, test='Chi')
# At 95% significant level, we may reject hypothesis that model 1 is better
# So treatment is not significant

# Now test significance of Gap
anova(richness_Y3_Y6_03, richness_Y3_Y6_05, test='Chi')
# At 95% significant level, we may reject hypothesis that model 5 is better
# So gap is not significant


###########
# Results #
###########

# Model selected by BIC is the simplest
# There is insufficient evidence to suggest a more complex model is better

# Model 3: Increase ~ Intercept + Fenced=True
summary(richness_Y3_Y6_03)

# Estimated odds of a fenced quadrat increasing in richness,
# between years 3 and 6 (vs unfenced quadrat):
exp(-0.8276)
# 0.44

# Estimated 95% CI for these odds:
c(exp(-0.8276-1.96*0.3635), exp(-0.8276+1.96*0.3635))
# [0.21, 0.89]

# Estimated probability of increase in richness
fenced <- predict(richness_Y3_Y6_03, richness_df[130,])
1/(1+exp(-1*fenced))
# 0.32

unfenced <- predict(richness_Y3_Y6_03, richness_df[1,])
1/(1+exp(unfenced))
# 0.68



##################################
# Richness Modeling: Year 0 to 6 #
##################################

# Model the probability that richness increased between years 3 and 6
richness_df$Y0.Y6.Increase <- richness_df$X6 > richness_df$X0
richness_df$Y0.Y6.Increase <- as.numeric(richness_df$Y0.Y6.Increase)

# Model 1: Treatment, Fenced, Gap, all first order interaction terms
richness_Y0_Y6_01 <- glm(Y0.Y6.Increase ~ (Treatment + Fenced + Gap)^2,
                         family='binomial', data=richness_df)
summary(richness_Y0_Y6_01)
# Nothing significant

# Model 2: Step-wise selection from Model 1 (AIC)
richness_Y0_Y6_02 <- step(richness_Y0_Y6_01)
summary(richness_Y0_Y6_02)

# Now In Gap is a significant factor. As is Treatment Gap:In Gap
# Some interaction still in the model but maybe not significant


# Model 3: Step-wise selection from Model 1 (BIC)
k = log(length(richness_df[,1]))
richness_Y0_Y6_03 <- step(richness_Y0_Y6_01, k=k)
summary(richness_Y0_Y6_03)

# Now just In Gap is a significant factor

# Test different models
anova(richness_Y0_Y6_01, richness_Y0_Y6_02, test='Chi')
# Smaller model OK (p-value 0.3)
anova(richness_Y0_Y6_01, richness_Y0_Y6_03, test='Chi')
# Smaller model OK, but less confident (p-value 0.06)
anova(richness_Y0_Y6_02, richness_Y0_Y6_03, test='Chi')
# Keep larger model 2 (p-value 0.45)
summary(richness_Y0_Y6_02)



# Model 4: Now without interaction
richness_Y0_Y6_04 <- glm(Y0.Y6.Increase ~ Treatment + Fenced + Gap,
                         family='binomial', data=richness_df)
summary(richness_Y0_Y6_04)
# Intercept and In Gap are significant


anova(richness_Y0_Y6_02, richness_Y0_Y6_04, test='Chi')
# Interaction between Fenced and Treatment may be important
summary(richness_Y0_Y6_02) # Specifically Treatment Gap : Fenced
# We can reject at 95% significance level, but only just (p-value 0.054)

# Model 5: Step-wise selection from Model 4 (AIC)
richness_Y0_Y6_05 <- step(richness_Y0_Y6_04)
summary(richness_Y0_Y6_05)
# Removed Fenced

# Model 6: Step-wise selection from Model 4 (BIC)
k = log(length(richness_df[,1]))
richness_Y0_Y6_06 <- step(richness_Y0_Y6_04)
summary(richness_Y0_Y6_06)
#Same result

# Test candidate models
anova(richness_Y0_Y6_05, richness_Y0_Y6_02, test='Chi')
# Can reject the hypothesis that the larger model is better (p-value 0.13)


###########
# Results #
###########

# Should be aware of the possibility that interaction is significant
# But based on this, we can exlcude interaction

# Model 5: Increase ~ Intercept + Treatment + In Gap
summary(richness_Y0_Y6_05)

# Estimated odds of an In Gap quadrat increasing in richness from Y0 to Y6
# (Compared to a Not In Gap quadrat)
exp(-0.6956) # 0.5

# 95% CI for these odds
c(exp(-0.6956-1.96*0.3066), exp(-0.6956+1.96*0.3066))
# [0.27, 0.91] # Significant
# In Gap less likely to increase in richness vs Not In Gap

# Estimated odds of a Gap Treatment quadrat increasing in richness from Y0 to Y6
# Compared to a Control quadrat
exp(0.4519) # 1.57

# 95% CI for these odds
c(exp(0.4519-1.96*0.3231), exp(0.4519+1.96*0.3231))
# [0.93, 3] ... Cannot say that this is significant

# Estimated odds of a Radial Treatment quadrat increasing in richness from Y0 to Y6
# Compared to a Control quadrat
exp(-0.3332) # 0.72

# 95% CI for these odds
c(exp(-0.3332-1.96*0.3821), exp(-0.3332+1.96*0.3821))
# [0.34, 1.5] ... Cannot say that this is significant



#######################
# Explore & Visualise #
#######################

cols1 <- c('springgreen3', 'orchid4', 'firebrick3')
cols2 <- c('mediumvioletred', 'midnightblue')

richness_df$Col1 <- cols1[1]
richness_df[richness_df$Treatment == 'Gap', 'Col1'] <- cols1[2]
richness_df[richness_df$Treatment == 'Radial', 'Col1'] <- cols1[3]

richness_df$Col2 <- cols2[1]
richness_df[richness_df$Fenced == FALSE, 'Col2'] <- cols2[2]
richness_df$Col3 <- cols2[1]
richness_df[richness_df$Gap == FALSE, 'Col3'] <- cols2[2]

richness_df_plot$Col1 <- cols1[1]
richness_df_plot[richness_df_plot$Treatment == 'Gap', 'Col1'] <- cols1[2]
richness_df_plot[richness_df_plot$Treatment == 'Radial', 'Col1'] <- cols1[3]

richness_df_plot$Col2 <- cols2[1]
richness_df_plot[richness_df_plot$Fenced == FALSE, 'Col12'] <- cols1[2]
richness_df_plot$Col3 <- cols2[1]
richness_df_plot[richness_df_plot$Gap == FALSE, 'Col3'] <- cols2[2]

line_plot <- function(df, label) {
  n <- length(df$Treatment)
  plot(0, type='l', xlim=c(0,6), ylim=c(-0.5,70.5),
       xlab='Year', ylab="Richness",
       main=paste(label, 'Plots'), col='white')
  for (i in 1:n) {
    lines(c(df[i,]$X0, df[i,]$X3, df[i,]$X6) ~ c(0,3,6),
          col=df[i,]$Col1)
  }
}



##############
# Line Plots #
##############

# All plots
line_plot(richness_df_plot, 'All')
# Control plots
line_plot(richness_df_plot[richness_df_plot$Treatment=='Control',], 'Control')
# Gap plots
line_plot(richness_df_plot[richness_df_plot$Treatment=='Gap',], 'Gap')
# Radial plots
line_plot(richness_df_plot[richness_df_plot$Treatment=='Radial',], 'Radial')


#####################
# Plot by Treatment #
#####################

# Change in Species Richness Y0 to Y6 (plots)
plot(richness_df_plot$X6-richness_df_plot$X0,
     col=richness_df_plot$Col1, pch=16, 
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y6")
lines(-25:300, rep(mean(richness_df_plot$X6-richness_df_plot$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in Species Richness Y0 to Y3 (plots)
plot(richness_df_plot$X3-richness_df_plot$X0,
     col=richness_df_plot$Col1, pch=16, 
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y3")
lines(-25:300, rep(mean(richness_df_plot$X3-richness_df_plot$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in Species Richness Y3 to Y6 (plots)
plot(richness_df_plot$X6-richness_df_plot$X3,
     col=richness_df_plot$Col1, pch=16,
      ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y3 to Y6")
lines(-25:300, rep(mean(richness_df_plot$X6-richness_df_plot$X3),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)




# Change in Species Richness Y0 to Y6 (quadrats)
plot(richness_df$X6-richness_df$X0,
     col=richness_df$Col1, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y6")
lines(-25:300, rep(mean(richness_df$X6-richness_df$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in Species Richness Y0 to Y3
plot(richness_df$X3-richness_df$X0,
     col=richness_df$Col1, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y3")
lines(-25:300, rep(mean(richness_df$X3-richness_df$X0),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Change in Species Richness Y3 to Y6
plot(richness_df$X6-richness_df$X3,
     col=richness_df$Col1, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y3 to Y6")
lines(-25:300, rep(mean(richness_df$X6-richness_df$X3),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)



# Magnitude of change in Species Richness Y0 to Y6
plot(abs(richness_df$X6-richness_df$X0),
     col=richness_df$Col1, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change',
     main="Change in Richness: Y0 to Y6")
lines(-25:300, rep(mean(abs(richness_df$X6-richness_df$X0)),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Magnitude of change in Species Richness  Y0 to Y3
plot(abs(richness_df$X3-richness_df$X0),
     col=richness_df$Col1, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change',
     main="Change in Richness: Y0 to Y3")
lines(-25:300, rep(mean(abs(richness_df$X3-richness_df$X0)),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)

# Magnitude of change in Species Richness  Y3 to Y6
plot(abs(richness_df$X6-richness_df$X3),
     col=richness_df$Col1, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change over time',
     main="Change in Richness: Y3 to Y6")
lines(-25:300, rep(mean(abs(richness_df$X6-richness_df$X3)),326),
      col='red')
legend('topleft', c('Control', 'Gap', 'Radial'), pch=16, bty='n',
       col=cols1)





##################
# Plot by Fenced #
##################

# Change in Species Richness Y0 to Y6 (quadrats)
plot(richness_df$X6-richness_df$X0,
     col=richness_df$Col2, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y6")
lines(-25:300, rep(mean(richness_df$X6-richness_df$X0),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)

# Change in Species Richness Y0 to Y3
plot(richness_df$X3-richness_df$X0,
     col=richness_df$Col2, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y3")
lines(-25:300, rep(mean(richness_df$X3-richness_df$X0),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)

# Change in Species Richness Y3 to Y6
plot(richness_df$X6-richness_df$X3,
     col=richness_df$Col2, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y3 to Y6")
lines(-25:300, rep(mean(richness_df$X6-richness_df$X3),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)



# Magnitude of change in Species Richness Y0 to Y6
plot(abs(richness_df$X6-richness_df$X0),
     col=richness_df$Col2, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change',
     main="Change in Richness: Y0 to Y6")
lines(-25:300, rep(mean(abs(richness_df$X6-richness_df$X0)),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in Species Richness  Y0 to Y3
plot(abs(richness_df$X3-richness_df$X0),
     col=richness_df$Col2, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change',
     main="Change in Richness: Y0 to Y3")
lines(-25:300, rep(mean(abs(richness_df$X3-richness_df$X0)),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in Species Richness  Y3 to Y6
plot(abs(richness_df$X6-richness_df$X3),
     col=richness_df$Col2, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change over time',
     main="Change in Richness: Y3 to Y6")
lines(-25:300, rep(mean(abs(richness_df$X6-richness_df$X3)),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)







##################
# Plot by In Gap #
##################

# Change in Species Richness Y0 to Y6 (quadrats)
plot(richness_df$X6-richness_df$X0,
     col=richness_df$Col3, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y6")
lines(-25:300, rep(mean(richness_df$X6-richness_df$X0),326),
      col='red')
legend('topleft', c('Unfenced', 'Fenced'), pch=16, bty='n',
       col=cols2)

# Change in Species Richness Y0 to Y3
plot(richness_df$X3-richness_df$X0,
     col=richness_df$Col3, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y0 to Y3")
lines(-25:300, rep(mean(richness_df$X3-richness_df$X0),326),
      col='red')
legend('topleft', c('Not in Gap', 'In Gap'), pch=16, bty='n',
       col=cols2)

# Change in Species Richness Y3 to Y6
plot(richness_df$X6-richness_df$X3,
     col=richness_df$Col3, pch=16,
     ylim=c(-20,30), ylab='Change over time',
     main="Change in Richness: Y3 to Y6")
lines(-25:300, rep(mean(richness_df$X6-richness_df$X3),326),
      col='red')
legend('topleft', c('Not in Gap', 'In Gap'), pch=16, bty='n',
       col=cols2)



# Magnitude of change in Species Richness Y0 to Y6
plot(abs(richness_df$X6-richness_df$X0),
     col=richness_df$Col3, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change',
     main="Change in Richness: Y0 to Y6")
lines(-25:300, rep(mean(abs(richness_df$X6-richness_df$X0)),326),
      col='red')
legend('topleft', c('Not in Gap', 'In Gap'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in Species Richness  Y0 to Y3
plot(abs(richness_df$X3-richness_df$X0),
     col=richness_df$Col3, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change',
     main="Change in Richness: Y0 to Y3")
lines(-25:300, rep(mean(abs(richness_df$X3-richness_df$X0)),326),
      col='red')
legend('topleft', c('Not in Gap', 'In Gap'), pch=16, bty='n',
       col=cols2)

# Magnitude of change in Species Richness  Y3 to Y6
plot(abs(richness_df$X6-richness_df$X3),
     col=richness_df$Col3, pch=16,
     ylim=c(-1,22), ylab='Magnitude of change over time',
     main="Change in Richness: Y3 to Y6")
lines(-25:300, rep(mean(abs(richness_df$X6-richness_df$X3)),326),
      col='red')
legend('topleft', c('Not in Gap', 'In Gap'), pch=16, bty='n',
       col=cols2)





