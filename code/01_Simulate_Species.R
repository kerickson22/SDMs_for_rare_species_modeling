#01a_Visualize Simulated Species

# Script for visualizing simulated species to
# check that they are well-defined

# Outputs of this script:
# *south.csv

path <- "/Users/curculion/Documents/GitHub"
source(paste0(path, "/SDMs_for_rare_species_modeling/code/00_Constants.R"))

## Load data from Norberg et al.
## Norberg data #####
DD <- "~/Documents/GitHub/rare_species/Norberg_Code/DATA"
#Spatial information (coordinates)
S_t <- read.csv(file.path(DD,"St_1_trees.csv"), header = F)
S_v <- read.csv(file.path(DD, "Sv_1_trees.csv"), header = F)
S_f <- rbind(S_t, S_v)
names(S_f) <- c("long", "lat")

#Environmental covariates (PC1, PC2 and PC3)
X_t <- read.csv(file.path(DD,"Xt_1_trees.csv"), header = F)
X_v <- read.csv(file.path(DD,"Xv_1_trees.csv"), header = F)
X_f <- rbind(X_t, X_v)
names(X_f) <- c("PC1", "PC2", "PC3")

#Presence-Absence Data
Y_t <- read.csv(file.path(DD,"Yt_1_trees.csv"), header = F)
Y_v <- read.csv(file.path(DD,"Yv_1_trees.csv"), header = F)
Y_f <- rbind(Y_t, Y_v)
names(Y_f) <- paste0("sp", 1:63)

dat <- data.frame(S_f, X_f, Y_f)

## Spatial Extant of Study #####
# * Clip data from Norberg et al. to SE US
states <- readRDS(paste0(path, "/rare_species/data/GADM_3.6/gadm36_USA_1_sp.rds"))

SE_states <- c("Oklahoma", "Arkansas", "Tennessee", "North Carolina", "South Carolina", "Texas", "Louisiana", "Mississippi", "Alabama", "Georgia", "Florida")

SE <- states[states@data$NAME_1 %in% SE_states,]

temp <-dat
coordinates(temp) <- c(1,2)
proj4string(temp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(SE) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

temp <- temp[SE,]

#Check that no coastlines were clipped:
index <- row.names(temp)
index_full <- row.names(dat)
not <- which(!index_full %in% index)
temp2 <- dat[not,]
plot(temp)
points(temp2, col="red")

south <- data.frame(temp)
north <- data.frame(temp2)
south <- south[,1:68]


## Real Species Niches #####
X_bar <- data.frame(mu_PC1 = rep(0, 63),
                    sd_PC1 = rep(0, 63),
                    mu_PC2 = rep(0, 63),
                    sd_PC2 = rep(0, 63),
                    mu_PC3 = rep(0, 63),
                    sd_PC3 = rep(0, 63),
                    N = rep(NA, 63))


rownames(X_bar) <- paste0("sp", 1:63)

for (j in 1:63) {
  temp <-subset(south, south[,j+5]==1)
  X_bar$mu_PC1[j] <- mean(temp$PC1) #mu_PC1
  X_bar$sd_PC1[j] <- sd(temp$PC1) #sd_PC1
  X_bar$min_PC1[j] <- min(temp$PC1)
  X_bar$max_PC1[j] <- max(temp$PC1)

  X_bar$mu_PC2[j] <- mean(temp$PC2) #mu_PC2
  X_bar$sd_PC2[j] <- sd(temp$PC2) #sd_PC2
  X_bar$min_PC2[j] <- min(temp$PC2)
  X_bar$max_PC2[j] <- max(temp$PC2)

  X_bar$mu_PC3[j] <- mean(temp$PC3) #mu_PC3
  X_bar$sd_PC3[j] <- sd(temp$PC3) #sd_PC3
  X_bar$min_PC3[j] <- min(temp$PC3)
  X_bar$max_PC3[j] <- max(temp$PC3)
  X_bar$N[j] <- length(temp$PC1)

}

mu_PC1 <- mean(X_bar$mu_PC1)
sd_PC1 <-mean(X_bar$sd_PC1, na.rm=T)
mu_PC2 <- mean(X_bar$mu_PC2)
sd_PC2 <- mean(X_bar$sd_PC2, na.rm=T)
mu_PC3 <- mean(X_bar$mu_PC3)
sd_PC3 <- mean(X_bar$sd_PC3[1:63], na.rm=T)


# Remove real species from our analysis with
# fewer than 4 occurrences

south <- south[, -(which(X_bar$N == 1)+5)]
write.csv(south, file="../data/south.csv", row.names=FALSE)

X_bar <- subset(X_bar, X_bar$N >=4)
mu_PC1 <- mean(X_bar$mu_PC1)
sd_PC1 <-mean(X_bar$sd_PC1, na.rm=T)
mu_PC2 <- mean(X_bar$mu_PC2)
sd_PC2 <- mean(X_bar$sd_PC2, na.rm=T)
mu_PC3 <- mean(X_bar$mu_PC3)
sd_PC3 <- mean(X_bar$sd_PC3[1:63], na.rm=T)

## Virtual Species Niches #####


mu_avg <- c(mu_PC1, mu_PC3)
sd_broad <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.90, na.rm=T),
                        quantile(X_bar$sd_PC3, probs=0.90, na.rm=T)))

#Now, to decide on the parameters associated
# with narrow niche and extreme niche position.
# To do this, we will explore a range of possible
# choices, and pick the one such that the probability
# of generating
sd_narrow <- sd_broad*niche_reduction_factor

#Pick

tries_PC1 <- c(seq(from=-9, to=-2.5, by=0.5), seq(from=-2, to=1, by=0.5))
tries_PC3 <- c(seq(from=-10, to=-1, by=0.5), seq(from=1, to=3, by=0.5))
tries_niche_reduction_factor <- seq(from=0.1, to = 0.9, by=0.1)
tries <- expand.grid(tries_PC1,
                     tries_PC3,
                     tries_niche_reduction_factor)
names(tries) <- c(
  "mu_PC1", "mu_PC3",
  "niche_reduction_factor"
)

df_tries <- data.frame(
  mu_PC1 = tries$mu_PC1,
  mu_PC3 = tries$mu_PC3,
  niche_reduction_factor = tries$niche_reduction_factor,
  p_narrow_avg = rep(NA, length(tries$mu_PC1)),
  p_broad_ext = rep(NA, length(tries$mu_PC1)),
  p_narrow_ext = rep(NA, length(tries$mu_PC1)))

for (i in 1:length(df_tries$mu_PC1)) {
  mu_ext <- c(df_tries$mu_PC1[i], df_tries$mu_PC3[i])
  sd_narrow <- df_tries$niche_reduction_factor[i] * sd_broad

  L_narrow_avg <- dnorm(south$PC1, mean=mu_avg[1], sd=sd_narrow[1,1])*
      dnorm(south$PC3, mean=mu_avg[2], sd=sd_narrow[2,2])
  L_narrow_avg <- L_narrow_avg/sum(L_narrow_avg)
  L_narrow_avg <- sort(L_narrow_avg, decreasing=T)

  L_narrow_ext <- dnorm(south$PC1, mean=mu_ext[1], sd=sd_narrow[1,1])*
    dnorm(south$PC3, mean=mu_ext[2], sd=sd_narrow[2,2])
  L_narrow_ext <- L_narrow_ext/sum(L_narrow_ext)
  L_narrow_ext <- sort(L_narrow_ext, decreasing=T)

  L_wide_ext <- dmvnorm(south[,c(3,5)],
                        mean=mu_ext,
                        sigma=Sigma_wide )
  L_wide_ext <- L_wide_ext/sum(L_wide_ext)
  L_wide_ext <- sort(L_wide_ext, decreasing=T)

  df_tries$p128_narrow_avg[ii] <- L_narrow_avg[128]/max(L_narrow_avg)
  df_tries$p128_narrow_ext[ii] <- L_narrow_ext[128]/max(L_narrow_ext)
  df_tries$p128_wide_ext[ii] <-L_wide_ext[128]/max(L_wide_ext)

  if (ii %% 10000 == 0) { cat(paste0(ii, " of ", length(df_tries$mu_PC1), "\n"))
  }
}
