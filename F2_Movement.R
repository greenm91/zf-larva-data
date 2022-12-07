## Load libraries
library(bannerCommenter) # for comments
library(car) # for ANOVA
library(MASS) # for box-cox transformation
library(lme4)# for modelling (glmer.nb)
library(RVAideMemoire) # to test overdispersion
library(performance) # to perform model diagnostics
library(glmmTMB) # for zero inflated nb glmm
library(paletteer) # for colour palettes
#
filepath <- file.choose()
## Read CSV from file path, use file encoding to ensure correct format.
raw.data <- read.csv(filepath,header=TRUE)
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
boxplot(raw.data$dist_cm ~ raw.data$g_treat*raw.data$p_treat)
boxplot(raw.data$mov ~ raw.data$g_treat*raw.data$p_treat)
par(mfrow=c(1,1))
banner("Variances are different")
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
nb.mov <- glmer.nb(mov ~ g_treat*p_treat + (1|p_id) + (1|g_id) + (1|m_id), data=raw.data)
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
nb.dist <- glmer.nb(dist_cm ~ g_treat*p_treat + (1|p_id) + (1|g_id) + (1|m_id), data=raw.data)
#
## NB.GLMM summary
summary(nb.dist)
car::Anova(nb.dist)
#
## Visually check diagnostics
check_model(nb.dist)
#
## Set level order for plots
raw.data$gp_treat<-factor(raw.data$gp_treat, levels=c('EE','SE','ES','SS'))
#
## Plots
banner("Plots", emph = TRUE)
#
## Distance
banner("Total Distance")
dist.plot <- ggplot(raw.data, aes(x = gp_treat, y = dist_cm, colour = gp_treat)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = gp_treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 1, 
                                             dodge.width = 1), aes(shape = gp_treat, fill = gp_treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("EE" = "ENR-ENR", "ES" = "ENR-STD", "SE" = "STD-ENR", "SS" = "STD-STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  ylim(0,200) +
  ylab("Total Distance Moved (cm)") +
  xlab("") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 18, hjust = 1),
        axis.title = element_text(size = 18),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 3),
        legend.position="none",
        aspect.ratio=10/8,
        plot.margin = margin(1, 1, 1, 1, "cm"))
dist.plot
#
##
banner("Movement Events")
mov.plot <- ggplot(raw.data, aes(x = gp_treat, y = mov, colour = gp_treat)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = gp_treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 1, 
                                             dodge.width = 1), aes(shape = gp_treat, fill = gp_treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("EE" = "ENR-ENR", "ES" = "ENR-STD", "SE" = "STD-ENR", "SS" = "STD-STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  ylim(0,100) +
  ylab("Number of Movement Events") +
  xlab("") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 18, hjust = 1),
        axis.title = element_text(size = 18),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 3),
        legend.position="none",
        aspect.ratio=10/8,
        plot.margin = margin(1, 1, 1, 1, "cm"))
mov.plot

mov.plot+dist.plot

