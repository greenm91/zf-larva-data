library(MASS) # FOR GLM.NB
library(ggpubr) # for compare_means and plotting sig.
library(lme4) # for modelling
library(car) # functionality
library(sjPlot) # plotting model estimates
library(sjmisc) # addon for sjplot
library(reshape2) # for melt
library(Rmisc) # summarySE
library(ggplot2) # plotting
library(pscl) # for ZI poisson
library(RVAideMemoire) # test overdispersion
library(dplyr) # for piping
library(emmeans) # contrasts
library(patchwork)

#
# Open the file for analysis using file.choose
# Use file encoding to ensure correct format.
#
filepath <- file.choose()
raw.data <- read.csv(filepath,header=TRUE, sep=",",fileEncoding="UTF-8-BOM")
#
raw.data$comb <- paste(raw.data$parental, "-", raw.data$offspring)
#
cols <- c(1:4)
varlist <- names(raw.data)[cols]
datacols <- c(10:11)
extract_cols <- c(1:4,10:11)
#
## Full labels
#
full.labs <- c("Distance (cm)","Movement Events")
names(full.labs) <- c("dist_cm", "mov")
#
stack.data <- raw.data[extract_cols]
stack <- melt(stack.data, id = c("parental","offspring","clutch","id"))
#
ggplot(stack, aes(x = value)) +
  geom_histogram(color = "black",fill = "white",bins = 10,alpha = 0.5) +
  xlab("Value") + 
  ylab("Frequency") +
  facet_wrap( ~ variable, scale = "free", labeller = labeller(variable=full.labs))
#
# Calculate clutch means
#
clutch.means <- aggregate(raw.data[c("dist_cm","mov")], 
                        by = raw.data[c("parental","clutch")],
                        FUN=mean, na.rm=TRUE)
#
#
# Distance
#
## Try GLMER
#
dist.glmer <- glmer(dist_cm ~ parental*offspring + (1|clutch), data = raw.data, family = poisson(link=log))
overdisp.glmer(dist.glmer)
#
# This model is overdispersed. 
# Try glm.nb to account for this
#
nb.dist <- glmer.nb(dist_cm ~ parental*offspring + (1|clutch), data=raw.data)
overdisp.glmer(nb.dist)
#
## Overdispersion is now fine.
#
plot(residuals(nb.dist))
qqnorm(residuals(nb.dist), main = "QQ Plot - Distance");qqline(residuals(nb.dist))
#
## Residuals look good
#
# effects
#
summary(nb.dist)
car::Anova(nb.dist)
#
plot_model(nb.dist,type="pred", terms=c("parental","offspring"))
#
#
## Movement Events
#
## Same as before. Although 
## I'm expecting it to be over-
## dispersed, its worth checking
#
mov.glmer <- glmer(mov ~ parental*offspring + (1|clutch), data = raw.data, family = poisson(link=log))
overdisp.glmer(mov.glmer)
#
## Really overdispersed again.
#
nb.mov <- glmer.nb(mov ~ parental*offspring + (1|clutch), data = raw.data, family = poisson(link=log), verbose=TRUE)
overdisp.glmer(nb.mov)
#
# Looks better
#
plot(residuals(nb.mov))
qqnorm(residuals(nb.mov), main = "QQ Plot - Distance");qqline(residuals(nb.mov))
#
## Residuals look okay.
#
summary(nb.mov)
car::Anova(nb.mov)
#
#######################
#
#
# Plotting
#
#
#######################
#
## Generate summaries
dist.summ <- summarySE(raw.data,measurevar="dist_cm",groupvars = c("parental","offspring"))
dist.summ

mov.summ <- summarySE(raw.data,measurevar="mov",groupvars = c("parental","offspring"))
mov.summ

mov.plot <- ggplot(mov.summ, aes(x=offspring, y=mov, group=parental, shape=parental, colour = parental)) +
              geom_errorbar(aes(ymin = mov - se, ymax = mov + se), 
                            width=0.05,size = 0.4) +
              geom_point(fill="white", size = 2) +
              geom_line(size = 0.45) +
              ylab("Number of Movement Events \u00b1 SE") +
              ylim(0,40) +
              xlab("F1") +
              scale_x_discrete(labels=c("E" = "ENR", "S" = "STD")) +
              scale_color_manual(name = 'F0', labels=c("E" = "ENR", "S" = "STD"), values=c('#E69F00', '#56B4E9')) +
              scale_shape_manual(name = 'F0', labels=c("E" = "ENR", "S" = "STD"), values=c(5,15)) +
              theme_classic()
mov.plot
dist.plot <- ggplot(dist.summ, aes(x=offspring, y=dist_cm, group=parental, shape=parental, colour = parental)) +
              geom_errorbar(aes(ymin = dist_cm - se, ymax = dist_cm + se), 
                            width=0.05,size = 0.4) +
              geom_point(fill="white", size = 2) +
              geom_line(size = 0.45) +
              ylab("Total Distance (cm) \u00b1 SE") +
              ylim(0,60) +
              xlab("F1") +
              scale_x_discrete(labels=c("E" = "ENR", "S" = "STD")) +
              scale_color_manual(name = 'F0', labels=c("E" = "ENR", "S" = "STD"), values=c('#E69F00', '#56B4E9')) +
              scale_shape_manual(name = 'F0', labels=c("E" = "ENR", "S" = "STD"), values=c(5,15)) +
              theme_classic()


tiff("F2_Locomotion.tiff", units="in", width=4.5, height=7.5, res=600)
dist.plot/mov.plot
dev.off()


ggplot(mov.summ, aes(x=offspring, y=mov, group=parental, shape=parental, colour = parental)) +
  geom_errorbar(aes(ymin = mov - se, ymax = mov + se), 
                width=0.05,size = 0.4) +
  geom_point(fill="white", size = 2) +
  geom_line(size = 0.45) +
  ylab("Number of Movement Events \u00b1 SE") +
  ylim(0,40) +
  xlab("F1") +
  scale_x_discrete(labels=c("E" = "ENR", "S" = "STD")) +
  scale_color_manual(name = 'F0', labels=c("E" = "ENR", "S" = "STD"), values=c('#E69F00', '#56B4E9')) +
  scale_shape_manual(name = 'F0', labels=c("E" = "ENR", "S" = "STD"), values=c(5,15)) +
  theme_classic()




tiff("movF2.tiff",units = "in", width = 3, height = 4, res = 2160)
ggplot(mov.summ, aes(x = parental, y = mov, fill = offspring)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = mov - se, ymax = mov + se),
                width=.2, position=position_dodge(.9)) +
  xlab("GE") + 
  ylab("Movement Events  ± s.e.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_discrete(name = "PE") +
  theme_minimal()
dev.off()

dist_plot <- ggplot(dist.summ, aes(x = parental, y = dist_cm, color = offspring)) +
  geom_point(aes(shape = offspring)) +
  geom_line(aes(group = offspring), size = 1) +
  geom_errorbar(aes(ymax=dist_cm+se, ymin=dist_cm-se), width=.05, size = 1) +
  xlab("") + 
  ylab("Total Distance (cm) +- s.e") + 
  theme(axis.text=element_text(size=12),legend.text=element_text(size=12),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.line = element_blank())
tiff("distF2.tiff",units = "in", width = 6.5, height = 4, res = 2160)
gg.gap(dist_plot, ylim = c(0,60), segments = c(1,5), tick_width = c(10), c(0.2,0,0.8))
add.legend(plot = dist_plot, margin = c(top=150,right=1,bottom=250,left=350))
dev.off()

mov_plot <- ggplot(mov.summ, aes(x = parental, y = mov, color = offspring)) +
  geom_point(aes(shape = offspring)) +
  geom_line(aes(group = offspring), size = 1) +
  geom_errorbar(aes(ymax=mov+se, ymin=mov-se), width=.05, size = 1) +
  xlab("") + 
  ylab("No. Movement Events +- s.e") + 
  theme(axis.text=element_text(size=12),legend.text=element_text(size=12),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        axis.line = element_blank())
tiff("movF2.tiff",units = "in", width = 6.5, height = 4, res = 2160)
gg.gap(mov_plot, ylim = c(0,45), segments = c(1,5), tick_width = c(10), c(0.2,0,0.8))
add.legend(plot = mov_plot, margin = c(top=150,right=1,bottom=250,left=350))
dev.off()
