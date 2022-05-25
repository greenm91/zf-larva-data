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
F1landmarkdata <- readland.tps("E:\\Data\\F1-LARVAE\\IMAGES\\F1\\F1.list.tps", readcurves = TRUE)
F2landmarkdata <- readland.tps("E:\\Data\\F2-LARVAE\\IMAGES\\F2.list.TPS", readcurves = TRUE)

# Partial warp scores
f1_pws <- read.csv("E:\\Data\\F1-LARVAE\\IMAGES\\F1\\p.csv",header=TRUE)
f2_pws <- read.csv("E:\\Data\\F2-LARVAE\\IMAGES\\p.csv",header=TRUE)

banner("Create classifier variables")

# F1
f1_classifier <- c(rep("c", 95),rep("e",94))
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

#
F1_PCA$x

#SUMMARY
summary(F1_PCA)

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
f1_model1 <- lmer(PW1 ~ P + (1|ID), data = f1_pws)
summary(f1_model1)
car::Anova(f1_model1)
block(paste("No main effects",
            collapse = " "), fold = TRUE)
# MODEL 2
f1_model2 <- lmer(PW2 ~ P + (1|ID), data = f1_pws)
summary(f1_model2)
car::Anova(f1_model2)
block(paste("No main effectss",
            collapse = " "), fold = TRUE)
# MODEL 3
f1_model3 <- lmer(PW3 ~ P + (1|ID), data = f1_pws)
summary(f1_model3)
car::Anova(f1_model3)
block(paste("MAIN EFFECT",
            collapse = " "), fold = TRUE)

banner("Create data frame for plotting individual warps")

# DEFINE COMPARISON
f1_comparison <- list(c("c","e"))
# CREATE DATA FRAME
f1_pws_data <- as.data.frame(f1_pws)
# BIND CLASSIFIER COLUMN
f1_pws_data <- cbind(f1_classifier,f1_pws_data)
# AS TIBBLE
f1_pws_tibble <- as_tibble(f1_pws_data)

banner("Calculate summary statistics")

# GENERATE SUMMARY STATS
f1_PW3_summary <- f1_pws_tibble %>%
  group_by(f1_classifier) %>%
  summarise(count = n(),
    mean_PW3 = mean(PW3),
    sd_PW3 = sd(PW3),
    n_PW3 = n(),
    SE_PW3 = sd(PW3)/sqrt(n()))
names(f1_PW3_summary)[names(f1_PW3_summary) == "mean_PW3"] <- "PW3"

banner("Plot individual warps")

# PERFORM STAT TEST FOR PLOT
pw3.stat.test <- compare_means(
  PW3 ~ P, data = f1_pws,
  method = "wilcox.test")%>%
  mutate(y.position = c(0.0035))

# PLOTS
## F1-PW3
p <- ggplot(f1_PW3_summary, aes(x=f1_classifier, y=PW3, group=1)) +
      geom_errorbar(aes(ymin = PW3 - SE_PW3, ymax = PW3 + SE_PW3), 
                    width=0.05,size = 0.4) +
      geom_point(fill="white", size = 2) +
      geom_line(size = 0.45) +
      geom_hline(yintercept = 0, linetype="twodash") + 
      stat_pvalue_manual(data = pw3.stat.test, label = "p.signif",
                         bracket.size = 0.4, size = 5) +
      ylab("PW3 Score \u00b1 SE") +
      xlab("Paternal Experience") +
      scale_x_discrete(labels=c("c" = "STD", "e" = "ENR")) +
      scale_color_manual(values=c('#E69F00', '#56B4E9')) +
      theme_classic()

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
f2_model1 <- lmer(PW1 ~ F0*F1 + (1|parental.id), data = f2_pws)
summary(f2_model1)
car::Anova(f2_model1)
block(paste("Interaction, no main effects",
            collapse = " "), fold = TRUE)
# MODEL 2
f2_model2 <- lmer(PW2 ~ F0*F1 + (1|parental.id), data = f2_pws)
summary(f2_model2)
car::Anova(f2_model2)
block(paste("Main effect of paternal condition",
            collapse = " "), fold = TRUE)
# MODEL 3
f2_model3 <- lmer(PW3 ~ F0*F1 + (1|parental.id), data = f2_pws)
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

# DEFINE COMPARISON
f2_comparison <- list(c("SS","SE","ES","EE"))
# CREATE DATA FRAME
f2_pws_data <- as.data.frame(f2_pws)
# BIND COMBINED TREATMENT COLUMN
f2_pws<- cbind(f2_classifier,f2_pws)
# AS TIBBLE
f2_pws_tibble <- as_tibble(f2_pws_data)

# SUMMARIES
# PW 1
f2_PW1_summary <- f2_pws_tibble %>%
  group_by(F0,F1) %>%
  summarise(count = n(),
            mean_PW1 = mean(PW1),
            sd_PW1 = sd(PW1),
            n_PW1 = n(),
            SE_PW1 = sd(PW1)/sqrt(n()))
names(f2_PW1_summary)[names(f2_PW1_summary) == "mean_PW1"] <- "PW1"
# PW 2
f2_PW2_summary <- f2_pws_tibble %>%
  group_by(F1,F0) %>%
  summarise(count = n(),
            mean_PW2 = mean(PW2),
            sd_PW2 = sd(PW2),
            n_PW2 = n(),
            SE_PW2 = sd(PW2)/sqrt(n()))
names(f2_PW2_summary)[names(f2_PW2_summary) == "mean_PW2"] <- "PW2"
# PW 3
f2_PW3_summary <- f2_pws_tibble %>%
  group_by(F1,F0) %>%
  summarise(count = n(),
            mean_PW3 = mean(PW3),
            sd_PW3 = sd(PW3),
            n_PW3 = n(),
            SE_PW3 = sd(PW3)/sqrt(n()))
names(f2_PW3_summary)[names(f2_PW3_summary) == "mean_PW3"] <- "PW3"

# PLOTS
## F2-PW1
p <- ggplot(f2_PW1_summary, aes(x=F1, y=PW1, group=F0, shape=F0, colour = F0)) +
  geom_errorbar(aes(ymin = PW1 - SE_PW1, ymax = PW1 + SE_PW1), 
                width=0.05,size = 0.4) +
  geom_point(fill="white", size = 2) +
  geom_line(size = 0.45) +
  geom_hline(yintercept = 0, linetype="twodash") + 
  ylab("PW1 Score \u00b1 SE") +
  xlab("F1") +
  scale_x_discrete(labels=c("E" = "ENR", "S" = "STD")) +
  scale_color_manual(labels=c("E" = "ENR", "S" = "STD"), values=c('#E69F00', '#56B4E9')) +
  scale_shape_manual(labels=c("E" = "ENR", "S" = "STD"), values=c(5,15)) +
  theme_classic()
##  F2-PW2
p2 <- ggplot(f2_PW2_summary, aes(x=F1, y=PW2, group=F0, shape=F0, colour = F0)) +
  geom_errorbar(aes(ymin = PW2 - SE_PW2, ymax = PW2 + SE_PW2), 
                width=0.05,size = 0.4) +
  geom_point(fill="white", size = 2) +
  geom_line(size = 0.45) +
  geom_hline(yintercept = 0, linetype="twodash") + 
  ylab("PW2 Score \u00b1 SE") +
  xlab("F1") +
  scale_x_discrete(labels=c("E" = "ENR", "S" = "STD")) +
  scale_color_manual(labels=c("E" = "ENR", "S" = "STD"), values=c('#E69F00', '#56B4E9')) +
  scale_shape_manual(labels=c("E" = "ENR", "S" = "STD"), values=c(5,15)) +
  theme_classic()
## F2-PW3
p3 <- ggplot(f2_PW3_summary, aes(x=F1, y=PW3, group=F0, shape=F0, colour = F0)) +
  geom_errorbar(aes(ymin = PW3 - SE_PW3, ymax = PW3 + SE_PW3), 
                width=0.05,size = 0.4) +
  geom_point(fill="white", size = 2) +
  geom_line(size = 0.45) +
  geom_hline(yintercept = 0, linetype="twodash") + 
  ylab("PW3 Score \u00b1 SE") +
  xlab("F1") +
  scale_x_discrete(labels=c("E" = "ENR", "S" = "STD")) +
  scale_color_manual(labels=c("E" = "ENR", "S" = "STD"), values=c('#E69F00', '#56B4E9')) +
  scale_shape_manual(labels=c("E" = "ENR", "S" = "STD"), values=c(5,15)) +
  theme_classic()
## 3-FIGURE PLOT
p/p2/p3


