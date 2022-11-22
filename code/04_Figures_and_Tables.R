#04_Figures_and_Tables

#Create figures for manuscript, as well as
# tables of information to include in manuscript

path <- "/Users/curculion/Documents/GitHub"
source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

# 1. Species Niches #####
plot(0)
ellipse_broad_avg <- data.frame(mixtools::ellipse(mu_avg, sd_broad, alpha=0.32, npoints=250), draw=F, newplot=F)

sampled_broad_avg   <- data.frame(mixtools::ellipse(c(estSp$broad_avg$X_bar2[1001, 1], estSp$broad_avg$X_bar2[1001, 5]), diag(c(estSp$broad_avg$X_bar2[1001, 2], estSp$broad_avg$X_bar2[1001,6])), alpha=0.32, npoints=250), draw=F, newplot=F)

ellipse_broad_ext <- data.frame(mixtools::ellipse(mu_ext, sd_broad, alpha=0.32, npoints=250), draw=F, newplot=F)

sampled_broad_ext   <- data.frame(mixtools::ellipse(c(estSp$broad_ext$X_bar2[1001, 1], estSp$broad_ext$X_bar2[1001, 5]), diag(c(estSp$broad_ext$X_bar2[1001, 2], estSp$broad_ext$X_bar2[1001,6])), alpha=0.32, npoints=250), draw=F, newplot=F)


ellipse_narrow_avg <- data.frame(mixtools::ellipse(mu_avg, sd_narrow, alpha=0.32, npoints=250), draw=F, newplot=F)

sampled_narrow_avg   <- data.frame(mixtools::ellipse(c(estSp$narrow_avg$X_bar2[1001, 1], estSp$narrow_avg$X_bar2[1001, 5]), diag(c(estSp$narrow_avg$X_bar2[1001, 2], estSp$narrow_avg$X_bar2[1001,6])), alpha=0.32, npoints=250), draw=F, newplot=F)



ellipse_narrow_ext <- data.frame(mixtools::ellipse(mu_ext, sd_narrow, alpha=0.32, npoints=250), draw=F, newplot=F)

sampled_narrow_ext   <- data.frame(mixtools::ellipse(c(estSp$narrow_ext$X_bar2[1001, 1], estSp$narrow_ext$X_bar2[1001, 5]), diag(c(estSp$narrow_ext$X_bar2[1001, 2], estSp$narrow_ext$X_bar2[1001,6])), alpha=0.32, npoints=250), draw=F, newplot=F)

#Calculate ellipses for real species
sp <- NULL

X_temp <-
  X_bar[order(-X_bar$sd_PC1),]
for(i in 1:62) {
  plot(0)
  sp[[i]] <- data.frame(mixtools::ellipse( as.numeric(X_temp[i, c(2,6)]), diag(X_temp[i, c(3,7)]), alpha=0.32, npoints=250), draw=F, newplot=F)
  sp[[i]]$area <- X_temp[i, 3] + X_temp[i, 7]
}

g  <-
  ggplot() +
  geom_point(data=south, aes(x=PC1, y=PC3), shape=4)

colsList <- brewer.pal(6, "Dark2")
#Sp55 is only found at one site
for ( i in 1:62) {
  g <- g +  geom_polygon(data=sp[[i]], aes(x=X1, y=X2),fill=NA,  col=colsList[5], alpha=0.25, size=0.25)
}

g1 <- g +
  geom_path(data=ellipse_broad_avg, aes(x=X1, y=X2, col="Central Generalist"), size=2)  +
  geom_polygon(data=sampled_broad_avg, aes(x=X1, y=X2, col="Central Generalist", fill="Central Generalist"), size=1, alpha=0.6) +
  geom_path(data=ellipse_broad_ext, aes(x=X1, y=X2, col="Marginal Generalist"), size=2)  +
  geom_polygon(data=sampled_broad_ext, aes(x=X1, y=X2, col="Marginal Generalist", fill="Marginal Generalist"), size=1, alpha=0.6) +
  geom_path(data=ellipse_narrow_avg, aes(x=X1, y=X2, col="Central Specialist"), size=2)  +
  geom_polygon(data=sampled_narrow_avg, aes(x=X1, y=X2, col="Central Specialist", fill="Central Specialist"), size=1, alpha=0.6) +
  geom_path(data=ellipse_narrow_ext, aes(x=X1, y=X2, col="Marginal Specialist"), size=2)  +
  geom_polygon(data=sampled_narrow_ext, aes(x=X1, y=X2, col="Marginal Specialist", fill="Marginal Specialist"), size=1, alpha=0.6) +
  scale_color_manual(values=colsList, name="Species:") +
  scale_fill_manual(values=colsList, name="Species:") +
  geom_point(data=south, aes(x=PC1, y=PC3), shape=4)



g1 <- formatPlot(g1, "left")
pdf("./Figures/01_real_species.pdf", width=6.5,
    height=6.5, pointsize=8,
    encoding="TeXtext.enc", onefile=F)
g1
dev.off()

png("./Figures/01_Species_Niches.png", width=6.5,
    height=6.5, units="in", res=600)
g1
dev.off()

# 2. Convergence #####

makeConvergencePlot <- function(dat, species_name, position) {
  p1 <- ggplot() +
    geom_bar(data=dat, aes(x=pos, y=num, fill=model), stat="identity") +
    scale_fill_manual(values=my_bright, name="Model:") +
    ggtitle(species_name) +
    xlab("Presences") +
    scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
    ylab("Converged Models")
  p1 <- formatPlot(p1, position)
  return(p1)
}

formatPlot <- function(p1, position) {
  if(position == "left") {
    p1 <- p1 + theme(axis.text.x = element_text(size=10),
                     plot.title = element_text(size=10),
                     panel.spacing=unit(0.05, "cm"),
                     strip.text = element_text(size=5),
                     legend.text =element_text(size=8),
                     legend.title=element_text(size=8))
  }
  if(position == "right") {
    p1 <- p1 +
      theme(axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size=10),
            plot.title = element_text(size=10),
            panel.spacing=unit(0.05, "cm"),
            strip.text = element_text(size=5),
            legend.text = element_text(size=8),
            legend.title =element_text(size=8))
  }
  p1 <- p1 + theme(legend.position="top")
  return(p1)
}

p1 <- makeConvergencePlot(sample_size_broad_avg, "   Central Generalist", "left")
p2 <- makeConvergencePlot(sample_size_broad_ext, "   Marginal Generalist", "right")
p3 <- makeConvergencePlot(sample_size_narrow_avg,"   Central Specialist", "right")
p4 <- makeConvergencePlot(sample_size_narrow_ext,"   Marginal Specialist", "right")
#p1 <- p1 + theme(plot.title =element_text(hjust=0.5, size=12))

# p2 <- p2 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
# p3 <- p3 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
# p4 <- p4 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
design1 <- "
5555
1234
1234"
design1b <- "
5555
5555
1234
1234
1234"


design2 <- "
1234"



pa <- (p1 | p2 | p3| p4) + guide_area() +
  plot_layout(design=design1, guides="collect", heights=c(0.5, 2))


png("./Figures/02_Convergence.png", width=6.5,
    height=2.5, units="in", res=600)
pa
dev.off()


#3. AUC and Tjurs R2 #####
makeDiscriminationPlot <- function(dat, metric, species_name, position) {
  if (metric == "AUC") {
    p1 <- ggplot() +
    stat_summary(data=dat, aes(x=pos, y=auc, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
    stat_summary(data=dat, aes(x=pos, y=auc, color=model), fun="median", geom="line") +
    geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright, name="Model") +
    scale_shape_manual(values=c(15:18,8), name="Model") +
    ggtitle(species_name) +
    xlab("Presences") +
    scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
    ylab( "AUC") +
    coord_cartesian(ylim=c(0.4, 0.9)) }
  if( metric == "Tjurs") {
    p1 <- ggplot() +
      stat_summary(data=dat, aes(x=pos, y=TjursR2, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
      stat_summary(data=dat, aes(x=pos, y=TjursR2, color=model), fun="median", geom="line") +
      scale_color_manual(values=my_bright, name="Model") +
      scale_shape_manual(values=c(15:18,8),  name="Model") +
      ggtitle(species_name) +
      xlab("Presences") +
      scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
      coord_cartesian(ylim=c(-0.01, 0.4)) +
      ylab( expression(paste("Tjurs R"^2)))
  }
  if(metric == "RMSE") {
    p1 <- ggplot() +
      stat_summary(data=dat, aes(x=pos, y=RMSEWeighted, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
      stat_summary(data=dat, aes(x=pos, y=RMSEWeighted, color=model), fun="median", geom="line") +
      geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright, name="Model") +
      scale_shape_manual(values=c(15:18,8), name="Model") +
      ggtitle(species_name) +
      xlab("Presences") +
      scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
      coord_cartesian(ylim=c(0.6, 1)) +
      ylab( "RMSE")
  }
  p1 <- formatPlot(p1, position)
  return(p1)
}

p1 <- makeDiscriminationPlot(dat_broad_avg_converged, "AUC", "(a) Central Generalist", "left")
p2 <- makeDiscriminationPlot(dat_broad_ext_converged, "AUC", "   Marginal Generalist", "right")
p3 <- makeDiscriminationPlot(dat_narrow_avg_converged, "AUC", "   Central Specialist", "right")
p4 <- makeDiscriminationPlot(dat_narrow_ext_converged, "AUC", "   Marginal Specialist", "right")
p5 <- makeDiscriminationPlot(dat_broad_avg_converged, "Tjurs", "(b) Central Generalist", "left")
p6 <- makeDiscriminationPlot(dat_broad_ext_converged, "Tjurs", "   Marginal Generalist", "right")
p7 <- makeDiscriminationPlot(dat_narrow_avg_converged, "Tjurs", "   Central Specialist", "right")
p8 <- makeDiscriminationPlot(dat_narrow_ext_converged, "Tjurs", "   Marginal Specialist", "right")


# p1 <- p1 + theme(plot.title =element_text(hjust=0.5, size=12))
#
# p2 <- p2 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
# p3 <- p3 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
# p4 <- p4 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
# p5 <- p5 +
#   theme(
#         plot.title =element_text(hjust=0.5, size=12))
# p6 <- p6 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
# p7 <- p7 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
# p8 <- p8 +
#   theme(axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.title =element_text(hjust=0.5, size=12))
#
#
# p1 <- p1 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#
#
# #top, right, bottom, left
#
# p2 <- p2 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
# #top, right, bottom, left
# p3 <- p3 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
# #top, right, bottom, left
# p4 <- p4 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
# #top, right, bottom, left
# p5 <- p5 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
# #top, right, bottom, left
# p6 <- p6 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
# #top, right, bottom, left
# p7 <- p7 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
# #top, right, bottom, left
# p8 <- p8 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
# #top, right, bottom, left
#
# pdf("./Figures/03_Discrimination_ALL_MODELS.pdf", width=6.5,
#     height=10.5, pointsize=8,
#     encoding="TeXtext.enc", onefile=F)
# ggarrange(p1, p2, p3, p4,
#           p5, p6, p7, p8,
#           widths=c(3.9, 3, 3, 3,
#                    3.9, 3, 3, 3),
#           ncol=4, nrow=2, common.legend=T) + theme(plot.margin = margin(0.1,0.15,0.1,0.15, "cm"))
# #top, right, bottom, left
# dev.off()

design1 <- "
5555
1234
1234"
design1b <- "
5555
5555
1234
1234
1234"


design2 <- "
1234"



design1 <- "
5555
1234
1234"
design1b <- "
5555
5555
1234
1234
1234"


design2 <- "
1234"



pa <- (p1 | p2 | p3| p4) + guide_area() +
  plot_layout(design=design1, guides="collect", heights=c(0.5, 2))


pb <- (p5 + theme(legend.position="none")| p6 + theme(legend.position="none") | p7 + theme(legend.position="none") | p8 + theme(legend.position="none")) +
  plot_layout(design=design2)

pdf("./Figures/03_Discrimination.pdf", width=6.5,
    height=2.5, pointsize=8,
    encoding="TeXtext.enc", onefile=F)
pa/pb
dev.off()

png("./Figures/02_Discrimination.png", width=6.5,
    height=5, units="in", res=600)
pa/pb  + plot_layout(heights=c(1.5,1))
dev.off()

#S3 RMSE ######
p1 <- makeDiscriminationPlot(dat_broad_avg_converged, "RMSE", "   Central Generalist", "left")
p2 <- makeDiscriminationPlot(dat_broad_ext_converged, "RMSE", "   Marginal Generalist", "right")
p3 <- makeDiscriminationPlot(dat_narrow_avg_converged, "RMSE", "   Central Specialist", "right")
p4 <- makeDiscriminationPlot(dat_narrow_ext_converged, "RMSE", "   Marginal Specialist", "right")

# pdf("./Figures/S1_RMSE.pdf", width=6.5,
#     height=3.5, pointsize=8,
#     encoding="TeXtext.enc", onefile=F)
# ggarrange(p1, p2, p3, p4, widths=c(3.9, 3, 3, 3),
#           ncol=4, nrow=1, common.legend=T) + theme(plot.margin = margin(0.1,0.15,0.1,0.1, "cm"))
# #top, right, bottom, left
# dev.off()

pa <- (p1 | p2 | p3| p4) + guide_area() +
  plot_layout(design=design1, guides="collect", heights=c(0.5, 2))


png("./Figures/S3_RMSE.png", width=6.5,
    height=2.5, units="in", res=600)
pa
dev.off()


# PC Calibration #####
p1 <- ggplot() +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=PC1_rankCor, color=model, shape=model), fun="mean", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=PC1_rankCor, color=model), fun="mean", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC1") +
  ggtitle("Central \n Generalist")

p2 <- ggplot() +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=PC1_rankCor, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=PC1_rankCor, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC1") +
  ggtitle("Marginal \n Generalist")

p3 <- ggplot() +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=PC1_rankCor, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=PC1_rankCor, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC1") +
  ggtitle("Central \n Specialist")

p4 <- ggplot() +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=PC1_rankCor, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=PC1_rankCor, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC1") +
  ggtitle("Marginal \n Specialist")

p5 <- ggplot() +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=varPC2Response, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=varPC2Response, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(0, 0.4)) +
  ylab( "SD PC2 Response") +
  ggtitle("Central \n Generalist")

p6 <- ggplot() +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=varPC2Response, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=varPC2Response, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(0, 0.4)) +
  ylab( "SD PC2 Response") +
  ggtitle("Marginal \n Generalist")

p7 <- ggplot() +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=varPC2Response, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=varPC2Response, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(0, 0.4)) +
  ylab( "SD PC2 Response") +
  ggtitle("Central \n Specialist")

p8 <- ggplot() +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=varPC2Response, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=varPC2Response, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(0, 0.4)) +
  ylab( "SD PC2 Response") +
  ggtitle("Marginal \n Specialist")


p9 <- ggplot() +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=PC3_rankCor, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=PC3_rankCor, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC3") +
  ggtitle("Central \n Generalist")

p10 <- ggplot() +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=PC3_rankCor, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=PC3_rankCor, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC3") +
  ggtitle("Marginal \n Generalist")

p11 <- ggplot() +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=PC3_rankCor, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=PC3_rankCor, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC3") +
  ggtitle("Central \n Specialist")

p12 <- ggplot() +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=PC3_rankCor, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=PC3_rankCor, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Quadratic", "Hmsc", "SAM"), name="Model") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  coord_cartesian(ylim=c(-1, 1)) +
  ylab( "Rank Correlation PC3") +
  ggtitle("Marginal \n Specialist")



p2 <- p2 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

p3 <- p3 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

p4 <- p4 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

p6 <- p6 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p7 <- p7 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p8 <- p8 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p10 <- p10 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p11 <- p11 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p12 <- p12 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())



pdf("./Figures/S5_PC_alignment.pdf", width=6.5,
    height=8, pointsize=8,
    encoding="TeXtext.enc", onefile=F)
ggarrange(p1,p2, p3, p4,p5, p6, p7, p8, p9, p10, p11, p12, widths=c(
  2,1.5, 1.5, 1.5,
  2, 1.5, 1.5, 1.5,
  2, 1.5, 1.5, 1.5),
  heights= c(
    2, 2, 2, 2,
    2, 1.2, 1.2, 1.2,
    2, 1.2, 1.2, 1.2
  ), ncol=4, nrow=3, common.legend=T)

dev.off()

# Calibration Barchart #####

cor_broad_avg_long <- read.csv("../data/cor_broad_avg_long.csv")
cor_broad_ext_long <- read.csv("../data/cor_broad_ext_long.csv")
cor_narrow_avg_long <- read.csv("../data/cor_narrow_avg_long.csv")
cor_narrow_ext_long <- read.csv("../data/cor_narrow_ext_long.csv")

makeCorPlot <- function(dat, title, PC, position) {

    if(PC == "PC1") {
      p1 <- ggplot(dat, aes(x=size, y=PC1Frac, fill=PC1)) +
        geom_col() +
        facet_grid(. ~ model) +
        scale_fill_manual(values=rwb2, name = "Correlation with true response:") +
        ylab("PC1: Proportion \n of Models") +
        ggtitle(title) +
        facet_grid(~ model,
                   scales="free_x",
                   space="free_x",
                   switch="x") +
        scale_x_discrete(name = "Presences", labels=c("2","", "8","", "32", ""))}
    if(PC == "PC3") {
      p1 <- ggplot(dat, aes(x=size, y=PC3Frac, fill=PC3)) +
        geom_col() +
        facet_grid(. ~ model) +
        scale_fill_manual(values=rwb2, name = "Correlation with true response:") +
        ylab("PC3: Proportion \n of Models") +
        ggtitle(title) +
        facet_grid(~ model,
                   scales="free_x",
                   space="free_x",
                   switch="x") +
        scale_x_discrete(name = "Presences", labels=c("2","", "8","", "32", ""))
    }
  p1 <- formatPlot(p1, position)
  return(p1)
}

p1 <- makeCorPlot(cor_broad_avg_percent, "(a) Central Generalist", "PC1",  "left")
p2 <- makeCorPlot(cor_broad_ext_percent, "   Marginal Generalist", "PC1", "right")
p3 <- makeCorPlot(cor_narrow_avg_percent, "   Central Specialist", "PC1", "right")
p4 <- makeCorPlot(cor_narrow_ext_percent, "   Marginal Specialist", "PC1", "right")
p5 <- makeCorPlot(cor_broad_avg_percent, "(b) Central Generalist", "PC3", "left")
p6 <- makeCorPlot(cor_broad_ext_percent, "   Marginal Generalist", "PC3", "right")
p7 <- makeCorPlot(cor_narrow_avg_percent, "   Central Specialist", "PC3", "right")
p8 <- makeCorPlot(cor_narrow_ext_percent, "   Marginal Specialist", "PC3", "right")


makePC2Plot <- function(dat, title, position) {

  p1 <- ggplot() +
    stat_summary(data=dat, aes(x=pos, y=varPC2Response, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
    stat_summary(data=dat, aes(x=pos, y=varPC2Response, color=model), fun="median", geom="line") +
    geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") +
    scale_color_manual(values=my_bright, name="Model:") +
    scale_shape_manual(values=c(15:18,8),  name="Model:") +
    xlab("Presences") +
    scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
    coord_cartesian(ylim=c(0, 0.25)) +
    ylab( "SD PC2 Response") +
    ggtitle(title)
  p1 <- formatPlot(p1, position)
  return(p1)
}

p9 <- makePC2Plot(dat_broad_avg_converged, "(c) Central Generalist", "left")
p10 <- makePC2Plot(dat_broad_ext_converged, "Marginal Generalist", "right")
p11 <- makePC2Plot(dat_narrow_avg_converged, "Central Specialist", "right")
p12 <- makePC2Plot(dat_narrow_ext_converged, "Marginal Specialist", "right")

# pa <- (p1 | p2 | p3 | p4) + plot_layout(guides = "collect") &
#   theme(legend.position = "top")
# pb <- (p5 + theme(legend.position="none")| p6 + theme(legend.position="none") | p7 + theme(legend.position="none") | p8 + theme(legend.position="none"))

design1 <- "
5555
1234
1234"
design1b <- "
5555
5555
1234
1234
1234"


design2 <- "
1234"



pa <- (p1 | p2 | p3| p4) + guide_area() +
  plot_layout(design=design1, guides="collect", heights=c(0.5, 2))


pb <- (p5 + theme(legend.position="none")| p6 + theme(legend.position="none") | p7 + theme(legend.position="none") | p8 + theme(legend.position="none")) +
  plot_layout(design=design2)

pc <- (p9 | p10 | p11| p12) + guide_area() +
  plot_layout(design=design1, guides="collect", heights=c(0.5, 2))

#pa/pb/pc


pdf("./Figures/04_Calibration.pdf", width=6.5,
    height=5, pointsize=8,
    encoding="TeXtext.enc", onefile=F)
pa/pb/pc + plot_layout(heights=c(1.5,1, 1.5))

dev.off()

png("./Figures/04_Calibration.png", width=6.5,
    height=7.5, units="in", res=600)
pa/pb/pc  + plot_layout(heights=c(1.5,1, 1.5))
dev.off()

# Plot real species #####

temp <- data.frame(count=colSums(south[,6:67]), species=colnames(south[6:67]))

png("./Figures/S1_real_species.png", width=6.5,
    height=6.5, units="in", res=600)
ggplot(temp, aes(x=reorder(species, count), y=count)) +
  geom_bar( stat="identity") + coord_flip() + xlab("Species") + ylab("Number of Presences") +
  theme(axis.text.y = element_text(size=6))
dev.off()

temp[] %>%
  arrange(count) %>%
  kable()



