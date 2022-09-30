#04_Figures_and_Tables

#Create figures for manuscript, as well as
# tables of information to include in manuscript

path <- "/Users/curculion/Documents/GitHub"
source(paste0(path, "/SDMs_for_rare_species_modeling/code/00_Constants.R"))

# 1. Species Niches

sp <- NULL

X_bar <-
  X_bar[order(-X_bar$sd_PC1),]
for(i in 1:dim(X_bar)[1]) {
  plot(0)
  sp[[i]] <- data.frame(mixtools::ellipse( as.numeric(X_bar[i, c("mu_PC1","mu_PC3")]), diag(X_bar[i, c("sd_PC1","sd_PC3")]), alpha=0.32, npoints=250), draw=F, newplot=F)
  sp[[i]]$area <- X_bar[i, "sd_PC1"] + X_bar[i, "sd_PC3"]
}

ellipse_broad_avg   <- data.frame(mixtools::ellipse(mu_avg, sd_broad, alpha=0.32, npoints=250), draw=F, newplot=F)
ellipse_broad_ext   <- data.frame(mixtools::ellipse(mu_ext, sd_broad, alpha=0.32, npoints=250), draw=F, newplot=F)
ellipse_narrow_avg <- data.frame(mixtools::ellipse(mu_avg, sd_narrow, alpha=0.32, npoints=250), draw=F, newplot=F)
ellipse_narrow_ext <- data.frame(mixtools::ellipse(mu_ext, sd_narrow, alpha=0.32, npoints=250), draw=F, newplot=F)

g  <-
  ggplot() +
  geom_point(data=south, aes(x=PC1, y=PC3), shape=4)

#Sp55 is only found at one site
for ( i in 1:dim(X_bar)[1]) {
  g <- g +  geom_polygon(data=sp[[i]], aes(x=X1, y=X2),fill=NA,  col=nicheCols[5], alpha=0.25, size=1)
}

g1 <- g +
  geom_path(data=ellipse_broad_avg, aes(x=X1, y=X2, col="Broad Average"), size=2) +
  geom_path(data=ellipse_narrow_avg, aes(x=X1, y=X2, col = "Narrow Average"), size=2) +
  geom_path(data=ellipse_broad_ext, aes(x=X1, y=X2, col = "Broad Extreme"), size=2) +
  geom_path(data=ellipse_narrow_ext, aes(x=X1, y=X2, col = "Narrow Extreme"), size=2) +
  scale_color_manual(values=nicheCols, name="Species")
g1 <- g1 +
  geom_point(data=south, aes(x=PC1, y=PC3), shape=4, color="black")

g1 <- g1 + theme(legend.position="top")
g2 <- g1 + theme_dark() + theme(legend.position="top")

