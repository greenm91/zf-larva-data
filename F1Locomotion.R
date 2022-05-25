library(MASS) # FOR GLM.NB
library(ggpubr) # for compare_means and plotting sig.
library(lme4) # for modelling
library(car) # functionality
library(reshape2) # for melt
library(Rmisc) # summarySE
library(pscl)
library(ggplot2) # plotting
library(RVAideMemoire) # test overdispersion
library(dplyr) # for piping
library(patchwork) # for ggplot figures
#
# Open the file for analysis using file.choose
# Use file encoding to ensure correct format.
#
filepath <- file.choose()
raw.data <- read.csv(filepath,header=TRUE, sep=",",fileEncoding="UTF-8-BOM")
#
## Define data columns
## and factor columns
cols <- c(9,10)
varlist <- names(raw.data)[cols]
datacols <- c(1:3,9:10)
raw.data <- raw.data[datacols]
#
## Create full labels
full.labs <- c("Distance (cm)","Movement Events")
names(full.labs) <- c("dist_cm", "mov")
#
## Create stack data for melting
stack.data <- raw.data
## Melt stack data
stack <- melt(stack.data, id = c("parental","clutch","id"))
#
## Plot distributions with
## melted data
ggplot(stack, aes(x = value)) +
  geom_histogram(color = "black",fill = "white",bins = 30, alpha = 0.5) +
  xlab("Value") + 
  ylab("Frequency") +
  facet_wrap( ~ variable, scale = "free", labeller = labeller(variable=full.labs))
#
## Calculate clutch means
clutch.means <- aggregate(raw.data[c("dist_cm","mov")], 
                        by = raw.data[c("parental","clutch")],
                        FUN=mean, na.rm=TRUE)
#
## FIT MODELS
## Total Distance
nb.dist <- glmer.nb(dist_cm ~ parental + (1|clutch), data=raw.data)
overdisp.glmer(nb.dist)
#
## Overdispersion is fine.
#
plot(residuals(nb.dist))
qqnorm(residuals(nb.dist), main = "QQ Plot - Distance");qqline(residuals(nb.dist))
#
## Residuals look good!
#
# effects
#
summary(nb.dist)
car::Anova(nb.dist)
#
#
## Movement Events
nb.mov <- glmer.nb(mov ~ parental + (1|clutch), data=raw.data, verbose=TRUE)
overdisp.glmer(nb.mov)
#
# Overdispersion is fine
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
## Movement events
mov.plot <- ggplot(raw.data, aes(x=parental, y=mov, 
                             colour = parental)) +
              geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = parental),
                           size = 1.4) +
              geom_jitter(width = 0.25, aes(fill = parental), colour = "black",
                          pch = 21, size = 2.5) +
              stat_compare_means(comparisons = f1_comparison, label.y = 820, size = 5)+
              scale_fill_manual(values=c('#E69F00', '#56B4E9')) +
              scale_colour_manual(values=c('#E69F00', '#56B4E9')) +
              scale_x_discrete(labels=c("S" = "STD", "E" = "ENR"), expand=c(0.6, 0)) +
              ylab("Number of Movement Events") +
              ylim(0,1000) +
              xlab("Paternal Experience") +
              theme_classic() + 
              theme(text = element_text(size = 20),
                    axis.text = element_text(size = 18),
                    axis.title = element_text(size = 20),
                    axis.text.x = element_text(vjust = -1),
                    axis.title.x = element_text(vjust = -3),
                    axis.text.y = element_text(hjust = 1),
                    axis.title.y = element_text(vjust = 3),
                    plot.margin = margin(1, 1, 1, 1, "cm"),
                    legend.position = "none")
# Total distance swam
dist.plot <- ggplot(raw.data, aes(x=parental, y=dist_cm, 
                                 colour = parental)) +
              geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = parental),
                           size = 1.4) +
              geom_jitter(width = 0.25, aes(fill = parental), colour = "black",
                          pch = 21, size = 2.5) +
              scale_fill_manual(values=c('#E69F00', '#56B4E9')) +
              scale_colour_manual(values=c('#E69F00', '#56B4E9')) +
              scale_x_discrete(labels=c("S" = "STD", "E" = "ENR"), expand=c(0.6, 0)) +
              ylab("Total Distance Swam (cm)") +
              ylim(0,500) +
              xlab("Paternal Experience") +
              theme_classic() + 
              theme(text = element_text(size = 20),
                    axis.text = element_text(size = 18),
                    axis.title = element_text(size = 20),
                    axis.text.x = element_text(vjust = -1),
                    axis.title.x = element_text(vjust = -3),
                    axis.text.y = element_text(hjust = 1),
                    axis.title.y = element_text(vjust = 3),
                    plot.margin = margin(1, 1, 1, 1, "cm"),
                    legend.position = "none")
#
## 2- figure plot
tiff("F1_figure.tiff", units="in", width=8, height=6.5, res=600)
mov.plot + dist.plot
dev.off()
