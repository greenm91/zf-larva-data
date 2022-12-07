## Load libraries
library(bannerCommenter) # for comments
library(car) # for ANOVA
library(MASS) # for box-cox transformation
library(lme4)# for modelling (glmer.nb)
library(RVAideMemoire) # to test overdispersion
library(performance) # to perform model diagnostics
library(glmmTMB) # for zero inflated nb glmm
#
## Use file.choose to get file path
filepath <- file.choose()
#
## Read CSV from file path, use file encoding to ensure correct format.
raw.data <- read.csv(filepath,header=TRUE, sep=",",fileEncoding="UTF-8-BOM")
#
## Plot distribution of response variables
par(mfrow=c(1,2))
hist(raw.data$dist_cm, breaks = 30)
hist(raw.data$mov, breaks = 30)
par(mfrow=c(1,1))
banner("Both responses are clearly very skewed")
#
## Formally check normality using shapiro-wilks test
shapiro.test(raw.data$dist_cm)
shapiro.test(raw.data$mov)
#
## Check variance of groups
par(mfrow=c(1,2))
boxplot(raw.data$dist_cm ~ raw.data$treat)
boxplot(raw.data$mov ~ raw.data$treat)
par(mfrow=c(1,1))
banner("Variances are very different")
#
## I am going to check the performance of the models using
## the performance package. Documentation on this model can 
## be found at https://easystats.github.io/performance/
##
banner("Movement Events", emph = TRUE)
#
##
banner("NB.GLMM")
#
## Fit a negative binomial GLMM
nb.mov <- glmer.nb(mov ~ treat + (1|p_id) + (1|m_id), data=raw.data)
#
## NB.GLMM summary
summary(nb.mov)
car::Anova(nb.mov)
#
## Visually check diagnostics
check_model(nb.mov)
#
##
banner("Total Distance", emph = TRUE)
#
##
banner("NB.GLMM")
#
## Fit a negative binomial GLMM
nb.dist <- glmer.nb(dist_cm ~ treat + (1|p_id/m_id), data=raw.data)
#
## NB.GLMM summary
summary(nb.dist)
car::Anova(nb.dist)
#
## Visually check diagnostics
check_model(nb.dist)
#