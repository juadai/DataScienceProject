#########
# Setup #
#########

# Load data
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# Load data
abundance_df <- read.csv('../dataset/tables/lf_abund_analysis.csv')
abundance_rel_df <- read.csv('../dataset/tables/lf_abund_rel_analysis.csv')

# Treatments, Plot numbers, Quadrat numbers, Life forms as factors
abundance_df$Treatment <- factor(abundance_df$Treatment)
abundance_df$Plot.Number <- factor(abundance_df$Plot.Number)
abundance_df$Quadrat.Number <- factor(abundance_df$Quadrat.Number)
abundance_df$Life.Form <- factor(abundance_df$Life.Form)

abundance_rel_df$Treatment <- factor(abundance_rel_df$Treatment)
abundance_rel_df$Plot.Number <- factor(abundance_rel_df$Plot.Number)
abundance_rel_df$Quadrat.Number <- factor(abundance_rel_df$Quadrat.Number)
abundance_rel_df$Life.Form <- factor(abundance_rel_df$Life.Form)

# Convert Fenced, Gap columns to boolean
abundance_df$Gap <- abundance_df$Gap=='True'
abundance_df$Fenced <- abundance_df$Fenced=='True'

abundance_rel_df$Gap <- abundance_rel_df$Gap=='True'
abundance_rel_df$Fenced <- abundance_rel_df$Fenced=='True'

# Relative Abundance: convert from proportion to percentage, and round to integer
# Note: Trace observations (below 1%) rounded up to 1%
abundance_rel_df$response <- abundance_rel_df$X0*100
abundance_rel_df[abundance_rel_df$response < 1 & abundance_rel_df$response != 0, 'response'] <- 1
abundance_rel_df$response <- round(abundance_rel_df$response)

# Weight for BIC
k <- log(length(abundance_rel_df$response))

###########################################
# Relative Abundance Modeling: Log-linear #
###########################################

model1 <- glm(response ~ (Treatment + Fenced + Gap + Life.Form + Year)^2, 
              family='poisson', data=abundance_rel_df)
summary(model1)

model2 <- step(model1)
summary(model2)
# Dropped Treatment:Gap, Fenced:Gap
# Retained 
#     Treatment:Fenced, Treatment:Life Form, Treatment:Year,
#     Fenced:Life Form, Fenced:Year, 
#     Gap:Life Form, Gap:Year,
#     Life Form:Year
anova(model1, model2, test='Chi')
# Can drop these 2 interaction terms

anova(glm(response ~ Treatment + Fenced + Gap + Life.Form + Year, 
          family='poisson', data=abundance_rel_df), model2, test='Chi')
# Can't drop the other interaction terms!


model3 <- glm(response ~ (Treatment + Fenced + Gap + Life.Form + factor(Year))^2, 
              family='poisson', data=abundance_rel_df)
model4 <- step(model3, k=k)
summary(model4)
# BIC selection with year as factor
#     Treatment, Fenced, Gap, Life Form Year
#     Treatment:Fenced, Treatment:Life Form, Treatment:Year
#     Fenced:Life Form, Fenced:Year
#     Gap:Life Form, Gap:Year
#     Life Form:Year
# No change from Model 2


# Is 3-way interaction at play here?
anova(
  glm(response ~ (Treatment + Fenced + Gap + Life.Form + Year)^3, 
      family='poisson', data=abundance_rel_df),
  glm(response ~ (Treatment + Fenced + Gap + Life.Form + Year)^2, 
      family='poisson', data=abundance_rel_df),
  test='Chi')
# Yes it is


# How about 4-way interaction terms?
anova(
  glm(response ~ (Treatment + Fenced + Gap + Life.Form + Year)^4, 
      family='poisson', data=abundance_rel_df),
  glm(response ~ (Treatment + Fenced + Gap + Life.Form + Year)^3, 
      family='poisson', data=abundance_rel_df),
  test='Chi')
# There should even be 4-way interaction in this model! 
# It's too complicated




# Big compute in this block
"
model5 <- step(
  glm(response ~ (Treatment + Fenced + Gap + Life.Form + Year)^3, 
    family='poisson', data=abundance_rel_df),
  k=k)
summary(model5)
# What 3-way interaction is in the model?
#     Treatment, Fenced, Gap, Life Form, Year

#     Treatment:Fenced, Treatment:Gap, Treatment:Life Form, Treatment:Year,
#     Fenced:Life Form, Fenced:Year,
#     Gap:Life Form, Gap:Year,
#     Life Form:Year

#     Treatment:Fenced:Life Form, Treatment:Fenced:Year
#     Treatment:Gap:Life Form             
#     Fenced:Life Form:Year
#     Gap:Life Form:Year

# (Can't drop Treatment:Gap as Treatment:Gap:Life Form is in the model)

model6 <- glm(formula = response ~ Treatment + Fenced + Gap + Life.Form + 
                Year + Treatment:Fenced + Treatment:Gap + Treatment:Life.Form + 
                Treatment:Year + Fenced:Life.Form + Fenced:Year + Gap:Life.Form + 
                Gap:Year + Life.Form:Year + Treatment:Fenced:Life.Form + 
                Treatment:Fenced:Year + Fenced:Life.Form:Year + 
                Gap:Life.Form:Year, family = 'poisson', data = abundance_rel_df)
anova(model6, model5, test='Chi')
# Treatment:Gap:Life Form is significant
"





# Just consider un-fenced quadrats, don't try to model browsing yet
abundance_rel_df <- abundance_rel_df[abundance_rel_df$Fenced==FALSE,]
k = log(length(abundance_rel_df$response))

abundance_model1 <- glm(response ~ (Treatment + Gap + Year + Life.Form)^2, 
                        family='poisson', data=abundance_rel_df)
abundance_model2 <- glm(response ~ (Treatment + Gap + Year + Life.Form)^3, 
                        family='poisson', data=abundance_rel_df)
summary(abundance_model2)
anova(abundance_model1, abundance_model2, test='Chi')

abundance_model3 <- step(abundance_model2, k=k)
summary(abundance_model3)
# There are still some third-order interaction terms
#       Treatment:Gap:Life.Form, Gap:Year:Life.Form

summary(glm(response ~ (Treatment + Gap + Year + Life.Form)^4, 
            family='poisson', data=abundance_rel_df))

anova(
  glm(response ~ (Treatment + Gap + Year + Life.Form)^4, 
      family='poisson', data=abundance_rel_df),
  glm(response ~ (Treatment + Gap + Year + Life.Form)^3, 
      family='poisson', data=abundance_rel_df),
  test='Chi'
)
# Still wants to keep the high order interaction terms

palette("R4")
plot(abundance_df$X0,
     col=abundance_df$Life.Form,
     pch=16)


# Can we get rid of some life forms 
# (i.e. simplify this factor, fewer parameters in modeling?)

lf <- unique(abundance_df$Life.Form)
n <- length(lf)
abund = rep(0,n)
for (i in 1:n) {
  abund[i] <- sum(abundance_df[abundance_df$Life.Form==lf[i], 'X0'])
}
names(abund) <- lf
plot(sort(abund, decreasing=TRUE))

sum(abund > 10)
sum(abund > 50)
sum(abund > 100)
sum(abund > 200)

rel_abund = rep(0,n)
for (i in 1:n) {
  rel_abund[i] <- sum(abundance_rel_df[abundance_rel_df$Life.Form==lf[i], 'X0'])
}
names(rel_abund) <- lf
round(rel_abund*100,4)
plot(sort(rel_abund, decreasing=TRUE))

sum(rel_abund > 1)
sum(rel_abund > 2.5)
sum(rel_abund > 5)
sum(rel_abund > 8)
sum(rel_abund > 10)

rel_abund <- rel_abund[rel_abund>10]
abundance_rel_df <- abundance_rel_df[abundance_rel_df$Life.Form %in% names(rel_abund),]


model1.1 <- glm(response ~ (Treatment + Fenced + Gap + Life.Form + Year)^3, 
                family='poisson', data=abundance_rel_df)
k=log(length(abundance_rel_df$response))
model2.1 <- step(model1.1, k=k)
summary(model2.1)

# It doesn't really help at all



