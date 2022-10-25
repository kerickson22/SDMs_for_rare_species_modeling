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
  geom_path(data=ellipse_broad_avg, aes(x=X1, y=X2, col="Generalist, central"), size=2)  +
  geom_polygon(data=sampled_broad_avg, aes(x=X1, y=X2, col="Generalist, central", fill="Generalist, central"), size=1, alpha=0.6) +
  geom_path(data=ellipse_broad_ext, aes(x=X1, y=X2, col="Generalist, marginal"), size=2)  +
  geom_polygon(data=sampled_broad_ext, aes(x=X1, y=X2, col="Generalist, marginal", fill="Generalist, marginal"), size=1, alpha=0.6) +
  geom_path(data=ellipse_narrow_avg, aes(x=X1, y=X2, col="Specialist, central"), size=2)  +
  geom_polygon(data=sampled_narrow_avg, aes(x=X1, y=X2, col="Specialist, central", fill="Specialist, central"), size=1, alpha=0.6) +
  geom_path(data=ellipse_narrow_ext, aes(x=X1, y=X2, col="Specialist, marginal"), size=2)  +
  geom_polygon(data=sampled_narrow_ext, aes(x=X1, y=X2, col="Specialist, marginal", fill="Specialist, marginal"), size=1, alpha=0.6) +
  scale_color_manual(values=colsList, name="Species") +
  scale_fill_manual(values=colsList, name="Species") +
  geom_point(data=south, aes(x=PC1, y=PC3), shape=4)


g1 <- g1 + theme(legend.position="top")

pdf("./Figures/01_real_species.pdf", width=6.5,
    height=6.5, pointsize=8,
    encoding="TeXtext.enc", onefile=F)
g1
dev.off()

# 2. Convergence #####

p1 <- ggplot() +
  geom_bar(data=sample_size_broad_avg, aes(x=pos, y=num, fill=model), stat="identity") +
  scale_fill_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Central \n Generalist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylab("Converged Models") +
  annotate("text", x=2, y=105, label="NOT FINAL", col="red")

p2 <- ggplot() +
  geom_bar(data=sample_size_broad_ext, aes(x=pos, y=num, fill=model), stat="identity") +
  scale_fill_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Marginal \n Generalist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylab("Converged Models")

p3 <- ggplot() +
  geom_bar(data=sample_size_narrow_avg, aes(x=pos, y=num, fill=model), stat="identity") +
  scale_fill_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Central \n Specialist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylab("Converged Models")

p4 <- ggplot() +
  geom_bar(data=sample_size_narrow_ext, aes(x=pos, y=num, fill=model), stat="identity") +
  scale_fill_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Marginal \n Specialist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylab("Converged Models")

p1 <- p1 + theme(plot.title =element_text(hjust=0.5, size=12))

p2 <- p2 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))

p3 <- p3 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))

p4 <- p4 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))



pdf("./Figures/02_Convergence.pdf", width=6.5,
    height=2.5, pointsize=8,
    encoding="TeXtext.enc", onefile=F)
ggarrange(p1,p2, p3, p4, widths=c(2.55, 2.17, 2.17, 2.17, 2.17), ncol=4, nrow=1, common.legend=T)
dev.off()

#3. AUC and Tjurs R2 #####

p1 <- ggplot() +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=auc, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=auc, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Central \n Generalist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(0.135, 1.06) + ylab( "AUC") +
  annotate("text", x=2, y=1.05, label="NOT FINAL: Need to rerun ESMs and Hmsc2:4", col="red")

p2 <- ggplot() +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=auc, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=auc, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Marginal \n Generalist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(0.135, 1.06) + ylab( "AUC")

p3 <-  ggplot() +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=auc, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=auc, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Central \n Specialist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(0.135, 1.06) + ylab( "AUC")


p4 <- ggplot() +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=auc, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=auc, color=model), fun="median", geom="line") +
  geom_hline(yintercept=0.5, linetype="dashed", size=.5, color="black") + scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Marginal \n Specialist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(0.135, 1.06) + ylab( "AUC")

p5 <- ggplot() +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=TjursR2, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_avg_converged, aes(x=pos, y=TjursR2, color=model), fun="median", geom="line") +
  scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Central \n Generalist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(-0.3, 0.6) + ylab( expression(paste("Tjurs R"^2)))


p6 <- ggplot() +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=TjursR2, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_broad_ext_converged, aes(x=pos, y=TjursR2, color=model), fun="median", geom="line") +
  scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Marginal \n Generalist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(-0.3, 0.6) + ylab( expression(paste("Tjurs R"^2)))


p7 <-  ggplot() +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=TjursR2, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_avg_converged, aes(x=pos, y=TjursR2, color=model), fun="median", geom="line") +
  scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Central \n Specialist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(-0.3, 0.6) + ylab( expression(paste("Tjurs R"^2)))


p8 <- ggplot() +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=TjursR2, color=model, shape=model), fun="median", fun.min=function(z) { quantile(z,0.1) }, fun.max=function(z) { quantile(z,0.9) }, size=0.5) +
  stat_summary(data=dat_narrow_ext_converged, aes(x=pos, y=TjursR2, color=model), fun="median", geom="line") +
  scale_color_manual(values=my_bright,labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  scale_shape_manual(values=c(15:18,8), labels=c("glm", "ESM:Linear", "ESM:Polynomial", "Hmsc", "SAM"), name="Model") +
  ggtitle("Marginal \n Specialist") +
  xlab("Presences") +
  scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6), labels=c("2", "4", "8", "16", "32", "64")) +
  ylim(-0.3, 0.6) + ylab( expression(paste("Tjurs R"^2)))

p1 <- p1 + theme(plot.title =element_text(hjust=0.5, size=12))

p2 <- p2 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))

p3 <- p3 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))

p4 <- p4 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))

p5 <- p5 +
  theme(
        plot.title =element_text(hjust=0.5, size=12))
p6 <- p6 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))

p7 <- p7 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))
p8 <- p8 +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title =element_text(hjust=0.5, size=12))


p1 <- p1 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))


#top, right, bottom, left

p2 <- p2 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#top, right, bottom, left
p3 <- p3 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#top, right, bottom, left
p4 <- p4 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#top, right, bottom, left
p5 <- p5 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#top, right, bottom, left
p6 <- p6 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#top, right, bottom, left
p7 <- p7 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#top, right, bottom, left
p8 <- p8 + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
#top, right, bottom, left

pdf("./Figures/03_Discrimination_ALL_MODELS.pdf", width=6.5,
    height=10.5, pointsize=8,
    encoding="TeXtext.enc", onefile=F)
ggarrange(p1, p2, p3, p4,
          p5, p6, p7, p8,
          widths=c(3.9, 3, 3, 3,
                   3.8, 3, 3, 3),
          ncol=4, nrow=2, common.legend=T) + theme(plot.margin = margin(0.1,0.15,0.1,0.1, "cm"))
#top, right, bottom, left
dev.off()