# Load packages

library(bannerCommenter)
library(geomorph)
library(lme4)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

banner("Load landmark data and warp scores")

# Landmark data
F1landmarkdata <- readland.tps("C:\\Users\\green\\OneDrive\\JEZB2022\\.tps\\F1_List.tps", readcurves = TRUE)
F2landmarkdata <- readland.tps("C:\\Users\\green\\OneDrive\\JEZB2022\\.tps\\F2_list.TPS", readcurves = TRUE)

# Partial warp scores
f1_pws <- read.csv("C:\\Users\\green\\Desktop\\F1_PartialWarpScores.csv",header=TRUE)
f2_pws <- read.csv("C:\\Users\\green\\Desktop\\F2_PartialWarpScores.csv",header=TRUE)

banner("Create classifier variables")

# F1
f1_classifier <- c(rep("S", 95),rep("E",94))
f1_classifierfactor <- as.factor(f1_classifier)

# F2
f2_classifier <- c(rep("SS", 37),rep("SE", 40), rep("ES", 39), rep("EE", 40))
f2_classifierfactor <- as.factor(f2_classifier)

banner("Define links and warp grid parameters")

# LINK MATRIX FOR WIREFRAME
links = matrix(c(1,2,3,4,5,6,6,8,9,10,
                 2,3,4,5,6,7,8,9,10,1),
               nrow = 10, ncol = 2)

# GRID PARAMETERS
GPar <- gridPar(
  pt.bg = "gray",  pt.size = 1.5,  link.col = "gray",  link.lwd = 2,
  link.lty = 1,  out.col = "gray",  out.cex = 0.1,  tar.pt.bg = "black",
  tar.pt.size = 0.75,  tar.link.col = "black",  tar.link.lwd = 1,
  tar.link.lty = 1,  tar.out.col = "black",  tar.out.cex = 0.1,
  n.col.cell = 25,  grid.col = "darkgreen",  grid.lwd = 1,  grid.lty = 2,
  txt.adj = NULL,  txt.pos = 1,  txt.cex = 1,  txt.col = "black")



banner("Section 1:", "F1 Geometric Morphometrics", emph = TRUE)

banner("Perform GPA and plot consensus")

# GPA
f1_gpa <- gpagen(F1landmarkdata)

#CONSENSUS
f1_ref <- mshape(f1_gpa$consensus)

# PLOT REFERENCE SHAPE
plotRefToTarget(f1_ref,f1_ref, mag = 1.2, links = links,  
                gridPars = GPar) +
  title("Consensus Shape", line = -1, cex.main = 2)

banner("Perform PCA")

# ANALYSIS
F1_PCA <- gm.prcomp(f1_gpa$coords)

#SUMMARY
summary(F1_PCA)

# COMPONENT LOADINGS
# head(F1_PCA$rotation)
# COMPONENT SCORES
# head(F1_PCA$x)

banner("Plot PCA Scores")

# CREATE DATA FRAME FOR PLOTTING WARP SCORES
F1_PCA_DATA <- as.data.frame(cbind(f1_classifierfactor,F1_PCA$x))

# PLOT RW SCORES
f1_pca_plot <- plot(F1_PCA_DATA$Comp1, F1_PCA_DATA$Comp2 , col = F1_PCA_DATA$f1_classifierfactor, 
     xlab = "PW1 Score", ylab = "PW2 Score",
     ylim = c(-0.05,0.05), xlim = c(-0.05,0.05),
     axes = F)
axis(1, c(-0.04,-0.02, 0, 0.02,0.04), 
     col = NA, col.ticks = 1,cex.axis = 0.75)
axis(2, c(-0.04,-0.02, 0, 0.02,0.04), 
     col = NA, col.ticks = 1,cex.axis = 0.75, las = 1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = 1, equilogs = TRUE)
abline(v=0, h=0, col="black")

banner("Plot shape deformations")

# COMP 1 MIN
plotRefToTarget(f1_ref,F1_PCA$shapes$shapes.comp1$min, mag = 1.2, links = links,
                gridPars = GPar) +
  title("Warp 1 Min", line = -1, cex.main = 2)

# COMP 1 MAX
plotRefToTarget(f1_ref,F1_PCA$shapes$shapes.comp1$max, mag = 1.5, links = links, 
                gridPars = GPar) +
  title("Warp 1 Max", line = -1, cex.main = 2)

# COMP 2 MIN
plotRefToTarget(f1_ref,F1_PCA$shapes$shapes.comp2$min, mag = 1.5, links = links,
                gridPars = GPar) +
  title("Warp 2 Min", line = -1, cex.main = 2)
# COMP 2 MAX
plotRefToTarget(f1_ref,F1_PCA$shapes$shapes.comp2$max, mag = 1.5, links = links, 
                gridPars = GPar) +
  title("Warp 2 Max", line = -1, cex.main = 2)

# COMP 3 MIN
plotRefToTarget(f1_ref,F1_PCA$shapes$shapes.comp3$min, mag = 2, links = links,
                gridPars = GPar) +
  title("Warp 3 Min", line = -1, cex.main = 2)

# COMP 3 MAX
plotRefToTarget(f1_ref,F1_PCA$shapes$shapes.comp3$max, mag = 2, links = links, 
                gridPars = GPar) +
  title("Warp 3 Max", line = -1, cex.main = 2)

banner("Fit models to first 3 warps")

# MODEL 1
f1_model1 <- lmer(PW1 ~ treat + (1|p_id) + (1|m_id), data = f1_pws)
summary(f1_model1)
car::Anova(f1_model1)
block(paste("No main effects",
            collapse = " "), fold = TRUE)

# MODEL 2
f1_model2 <- lmer(PW2 ~ treat + (1|p_id) + (1|m_id), data = f1_pws)
summary(f1_model2)
car::Anova(f1_model2)
block(paste("No main effectss",
            collapse = " "), fold = TRUE)

# MODEL 3
f1_model3 <- lmer(PW3 ~ treat + (1|p_id) + (1|m_id), data = f1_pws)
summary(f1_model3)
car::Anova(f1_model3)
block(paste("MAIN EFFECT",
            collapse = " "), fold = TRUE)

## Plots
banner("Plots", emph = TRUE)
#
## F1-PW1
banner("F1-PW1")
f1_pw1.plot <- ggplot(f1_pws, aes(x = treat, y = PW1, colour = treat)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, 
                                             dodge.width = 1), aes(shape = treat, fill = treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("E" = "ENR", "S" = "STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  geom_hline(yintercept = 0) +
  ylim(-0.03,0.03) +
  ylab("PW1 Score") +
  xlab("") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_text(vjust = -3),
        axis.text.y = element_text(hjust = 1),
        axis.title.y = element_text(vjust = 3),
        legend.position="none",
        aspect.ratio=10/8,
        plot.margin = margin(1, 1, 1, 1, "cm"))
#
## F1-PW2
banner("F1-PW2")
f1_pw2.plot <- ggplot(f1_pws, aes(x = treat, y = PW2, colour = treat)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, 
                                             dodge.width = 1), aes(shape = treat, fill = treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("E" = "ENR", "S" = "STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  geom_hline(yintercept = 0) +
  ylim(-0.03,0.03) +
  ylab("PW2 Score") +
  xlab("") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_text(vjust = -3),
        axis.text.y = element_text(hjust = 1),
        axis.title.y = element_text(vjust = 3),
        legend.position="none",
        aspect.ratio=10/8,
        plot.margin = margin(1, 1, 1, 1, "cm"))
#
## F1-PW3
banner("F1-PW3")
f1_pw3.plot <- ggplot(f1_pws, aes(x = treat, y = PW3, colour = treat)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, 
                                             dodge.width = 1), aes(shape = treat, fill = treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("E" = "ENR", "S" = "STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  geom_hline(yintercept = 0) +
  ylim(-0.04,0.04) +
  ylab("PW3 Score") +
  xlab("") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_text(vjust = -3),
        axis.text.y = element_text(hjust = 1),
        axis.title.y = element_text(vjust = 3),
        legend.position="none",
        aspect.ratio=10/8,
        plot.margin = margin(1, 1, 1, 1, "cm"))
#
f1_pw1.plot+f1_pw2.plot+f1_pw3.plot
#
##
banner("Section 2:", "F2 Geometric Morphometrics", emph = TRUE)

# PERFORM GPA
f2_gpa <- gpagen(F2landmarkdata)

# DEFINE CONSENSUS
f2_ref <- mshape(f2_gpa$consensus)

# PLOT REFERENCE SHAPE
plotRefToTarget(f2_ref,f2_ref, mag = 1.2, links = links,  
                gridPars = GPar) +
                title("Consensus Shape", line = -1, cex.main = 2)

# GEOMETRIC PCA ON GPA COORDS
F2_PCA <- gm.prcomp(f2_gpa$coords)

# PCA VARIANCE EXPLAINED
summary(F2_PCA)

# MODEL 1
f2_model1 <- lmer(PW1 ~ g_treat*p_treat + (1|p_id/g_id) + (1|m_id), data = f2_pws)
summary(f2_model1)
car::Anova(f2_model1)
block(paste("Interaction, no main effects",
            collapse = " "), fold = TRUE)
# MODEL 2
f2_model2 <- lmer(PW2 ~ g_treat*p_treat + (1|p_id/g_id) + (1|m_id), data = f2_pws)
summary(f2_model2)
car::Anova(f2_model2)
block(paste("Main effect of paternal condition",
            collapse = " "), fold = TRUE)
# MODEL 3
f2_model3 <- lmer(PW3 ~ g_treat*p_treat + (1|p_id/g_id) + (1|m_id), data = f2_pws)
summary(f2_model3)
car::Anova(f2_model3)
block(paste("Main effect of paternal condition",
            collapse = " "), fold = TRUE)

# PLOT SHAPE DEFORMATIONS

# COMP 1 MIN
plotRefToTarget(f2_ref,F2_PCA$shapes$shapes.comp1$min, mag = 1.0, links = links,
                gridPars = GPar) +
  title("Warp 1 Min", line = 1, cex.main = 2)

# COMP 1 MAX
plotRefToTarget(f2_ref,F2_PCA$shapes$shapes.comp1$max, mag = 1.0, links = links, 
                gridPars = GPar) +
  title("Warp 1 Max", line = 1, cex.main = 2)

# COMP 2 MIN
plotRefToTarget(f2_ref,F2_PCA$shapes$shapes.comp2$min, mag = 1.0, links = links,
                gridPars = GPar) +
  title("Warp 2 Min", line = 1, cex.main = 2)

# COMP 2 MAX
plotRefToTarget(f2_ref,F2_PCA$shapes$shapes.comp2$max, mag = 1.0, links = links, 
                gridPars = GPar) +
  title("Warp 2 Max", line = 1, cex.main = 2)

# COMP 3 MIN
plotRefToTarget(f2_ref,F2_PCA$shapes$shapes.comp3$min, mag = 1.0, links = links,
                gridPars = GPar) +
  title("Warp 3 Min", line = 1, cex.main = 2)

# COMP 3 MAX
plotRefToTarget(f2_ref,F2_PCA$shapes$shapes.comp3$max, mag = 1.0, links = links, 
                gridPars = GPar) +
  title("Warp 3 Max", line = 1, cex.main = 2)
#
## Set level order for plots
f2_pws$gp_treat<-factor(f2_pws$gp_treat, levels=c('EE','SE','ES','SS'))
#
## Plots
banner("Plots", emph = TRUE)
#
## F2-PW1
banner("F2-PW3")
f2_pw1.plot <- ggplot(f2_pws, aes(x = gp_treat, y = PW1, colour = gp_treat)) +
  geom_hline(yintercept = 0, linetype = "twodash") +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = gp_treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 1, 
                                             dodge.width = 1), aes(shape = gp_treat, fill = gp_treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("EE" = "ENR-ENR", "ES" = "ENR-STD", "SE" = "STD-ENR", "SS" = "STD-STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  ylab("F2-PW1 Score") +
  ylim(-0.03,0.03) +
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
f2_pw1.plot
#
## F2-PW2
banner("F2-PW2")
f2_pw2.plot <- ggplot(f2_pws, aes(x = gp_treat, y = PW2, colour = gp_treat)) +
  geom_hline(yintercept = 0, linetype = "twodash") +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = gp_treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 1, 
                                             dodge.width = 1), aes(shape = gp_treat, fill = gp_treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("EE" = "ENR-ENR", "ES" = "ENR-STD", "SE" = "STD-ENR", "SS" = "STD-STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  ylab("F2-PW2 Score") +
  ylim(-0.03,0.02) +
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
f2_pw2.plot
#
## F2-PW3
banner("F2-PW3")
f2_pw3.plot <- ggplot(f2_pws, aes(x = gp_treat, y = PW3, colour = gp_treat)) +
  geom_hline(yintercept = 0, linetype = "twodash") +
  geom_boxplot(outlier.shape = NA, width = 0.8, aes(colour = gp_treat),
               position = position_dodge(width = 1), size = 1.2) +
  geom_point(position = position_jitterdodge(jitter.width = 1, 
                                             dodge.width = 1), aes(shape = gp_treat, fill = gp_treat), size = 2.5, pch = 21,
             colour = "black") +
  scale_x_discrete(labels=c("EE" = "ENR-ENR", "ES" = "ENR-STD", "SE" = "STD-ENR", "SS" = "STD-STD"), expand=c(0, 0)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_shape_manual(values=c(3, 16, 17, 4)) +
  ylab("F2-PW3 Score") +
  ylim(-0.06,0.04) +
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
f2_pw3.plot


f2_pw1.plot+f2_pw2.plot+f2_pw3.plot
