#00a_Run_Once

# This script should be run at the start of the project, and only run once.
if(Sys.info()['sysname'] == "Darwin") {
  path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
  path <- "C:/Users/kerickson/Documents/GitHub"
}

source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))


# Assemble Norberg et al. data #####

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
write.csv(X_bar, file="../data/X_bar.csv", row.names=F)



# Joint Hmsc #####
# 1) How many latent factors?
# 2) What are the posterior widths of the latent factors?
nChains <- 2
thin <- 1
samples <- 20000
transient <-2000
verbose <-1000

m <- Hmsc(Y=south[,6:67], XData=south[,3:5],
          XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
            I(PC2^2)+I(PC3^2),
          distr="probit")

m <- sampleMcmc(m, thin=thin, samples=samples,
                transient=transient,
                nChains=nChains, verbose=verbose)

save(m, file="../data/real_species_joint_Hmsc.RData")
load("../data/real_species_joint_Hmsc.RData")
mpost <- convertToCodaObject(m)

#species niches
betas <- MCMCsummary(mpost$Beta)
max(betas$Rhat) #1.09
MCMCplot(mpost$Beta, ref_ovl=TRUE)

MCMCplot(mpost$Gamma, ref_ovl = TRUE)
MCMCplot(mpost$Sigma, ref_ovl = TRUE)
#influence of traits on niches
#pg 110
# 7 of these
# Intercept, tr1
# PC1, tr1
# PC2, tr1
#...
gammas <- MCMCsummary(mpost$Gamma)
max(gammas$Rhat) #1.13
MCMCplot(mpost$Gamma, ref_ovl=TRUE)

#V: residual covariance of species niches
#pg 110 of book
#49 of these
#Intercept, Intercept
#Intercept, PC1
# Intercept, PC2
# Intercept, PC3
# Intercept, PC1^2
# Intercept, PC2^2
vs <- MCMCsummary(mpost$V)
MCMCplot(mpost$V, ref_ovl = TRUE)
max(vs$Rhat) #1.13


# Residual variance
# There are 61 of these, one for each species

sigmas <- MCMCsummary(mpost$Sigma)
max(sigmas$Rhat) #NaN, because there is no variance
MCMCplot(sigmas)

#Priors:
m$V0 #7 x 7 matrix, with 1s on diagonal, 0s elsewhere
m$f0 #8
m$mGamma  #vector of 7 zeroes
m$UGamma #same as m$V0
m$aSigma #62 ones
m$bSigma #62 "0.01"s
m$nu #NULL
m$a1 #NULL
m$b1 #NULL
m$a2 #NUll
m$b2 #NULL
m$rhopw #101 x 2
m$nuRRR #3
m$a1RRR #1
m$b1RRR #1
m$a2RRR #50
m$b2RRR #1

# With Latent Factors #####
studyDesign <- data.frame(sample=as.factor(1:dim(south)[1]))
rL = HmscRandomLevel(units = studyDesign$sample)
#Constrain latent factors to 1
rL$nfMin =1
rL$nfMax =1

m <- Hmsc(Y=south[,6:67], XData=south[,3:5],
          XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
            I(PC2^2)+I(PC3^2), studyDesign=studyDesign,
          ranLevels=list(sample=rL),
          distr="probit")
mpost2 <- convertToCodaObject(m)
samples <- 1000
transient <- 500


m <- sampleMcmc(m, thin=thin, samples=samples,
                transient=transient,
                nChains=nChains, verbose=verbose)

# SAM #####

#Picking Simulate Species Values
#Now, to decide on the parameters associated
# with narrow niche and extreme niche position.
# To do this, we will explore a range of possible
# choices, and pick the one such that the 128th pick
# of presences still has relatively high probability


tries_PC1 <- c(seq(from=-9, to=-2.5, by=0.5), seq(from=-2, to=1, by=0.5))
tries_PC3 <- c(seq(from=-10, to=-1, by=0.5), seq(from=1, to=3, by=0.5))
tries_sd_quantile_PC1 <- seq(from=0.1, to = 0.8, by=0.1)
tries_sd_quantile_PC3 <- seq(from=0.1, to = 0.8, by=0.1)

tries <- expand.grid(tries_PC1,
                     tries_PC3,
                     tries_sd_quantile_PC1,
                     tries_sd_quantile_PC3)
names(tries) <- c(
  "mu_PC1", "mu_PC3",
  "sd_quantile_PC1",
  "sd_quantile_PC3"
)

tries$dist_mu_PC1 <- abs(mu_PC1-tries$mu_PC1)
tries$dist_mu_PC3 <- abs(mu_PC3-tries$mu_PC3)

df_tries <- data.frame(
  mu_PC1 = tries$mu_PC1,
  dist_mu_PC1 = tries$dist_mu_PC1,
  mu_PC3 = tries$mu_PC3,
  dist_mu_PC3 = tries$dist_mu_PC3,
  sd_quantile_PC1 = tries$sd_quantile_PC1,
  sd_quantile_PC3 = tries$sd_quantile_PC3,
  p_128_narrow_avg = rep(NA, length(tries$mu_PC1)),
  p_128_broad_ext = rep(NA, length(tries$mu_PC1)),
  p_128_narrow_ext = rep(NA, length(tries$mu_PC1)))

for (i in 1:length(df_tries$mu_PC1)) {
  mu_ext <- c(df_tries$mu_PC1[i], df_tries$mu_PC3[i])
  sd_narrow <- diag(x=c(quantile(X_bar$sd_PC1, probs=df_tries$sd_quantile_PC1[i], na.rm=T),
                        quantile(X_bar$sd_PC3, probs=df_tries$sd_quantile_PC3[i], na.rm=T)))

  L_narrow_avg <- dnorm(south$PC1, mean=mu_avg[1], sd=sd_narrow[1,1])*
    dnorm(south$PC3, mean=mu_avg[2], sd=sd_narrow[2,2])
  L_narrow_avg <- L_narrow_avg/sum(L_narrow_avg)
  L_narrow_avg <- sort(L_narrow_avg, decreasing=T)

  L_narrow_ext <- dnorm(south$PC1, mean=mu_ext[1], sd=sd_narrow[1,1])*
    dnorm(south$PC3, mean=mu_ext[2], sd=sd_narrow[2,2])
  L_narrow_ext <- L_narrow_ext/sum(L_narrow_ext)
  L_narrow_ext <- sort(L_narrow_ext, decreasing=T)

  L_broad_ext <- dnorm(south$PC1, mean=mu_ext[1], sd=sd_broad[1,1])*
    dnorm(south$PC3, mean=mu_ext[2], sd=sd_broad[2,2])
  L_broad_ext <- L_broad_ext/sum(L_broad_ext)
  L_broad_ext <- sort(L_broad_ext, decreasing=T)

  df_tries$p_128_narrow_avg[i] <- L_narrow_avg[128]/max(L_narrow_avg)
  df_tries$p_128_narrow_ext[i] <- L_narrow_ext[128]/max(L_narrow_ext)
  df_tries$p_128_broad_ext[i] <-L_broad_ext[128]/max(L_broad_ext)

  if (i %% 3000 == 0) { cat(paste0(i, " of ", length(df_tries$mu_PC1), "\n"))
  }
}

df_tries2 <- subset(df_tries, df_tries$p_128_narrow_avg > 0.6 &
                      df_tries$p_128_broad_ext > 0.6 & df_tries$p_128_narrow_ext>0.6)

df_tries2 <- df_tries2[order(df_tries2$sd_quantile_PC1, df_tries2$sd_quantile_PC3, -df_tries2$dist_mu_PC1, -df_tries2$dist_mu_PC3),]

head(df_tries2)
df_tries3 <- subset(df_tries2, df_tries2$mu_PC1 == -3.5 &
                      df_tries2$mu_PC3 == 2)
df_tries4 <- subset(df_tries3,
                    df_tries3$sd_quantile_PC1==0.4)
#sd_narrow <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.1, na.rm=T),
# quantile(X_bar$sd_PC3, probs=0.8, na.rm=T)))
sd_narrow <- diag(c(0.9, 2.5))
mu_ext <- c(-3.5, 2)
sd_broad <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.90, na.rm=T),
                     quantile(X_bar$sd_PC3, probs=0.90, na.rm=T)))
#Discussion point: choose different values of sd?

# ESM
verboten <- c('PC1:PC2', 'PC2:PC3', "PC1:PC3")
formulae <- makeFormulae(Y ~ PC1 + PC2 +
                           PC3+I(PC1^2) +
                           I(PC2^2)+I(PC3^2),
                         verboten=verboten,
                         interceptOnly=F,
                         intercept=T)
## Fill out Formula Matrix with coefficient
# positions. (Hand-scored)
# formulaMatrix <- data.frame(
#   formula = as.character(formulae),
#   Int = rep(NA, 26),
#   PC1 = rep(NA, 26),
#   PC2 = rep(NA, 26),
#   PC3 = rep(NA, 26),
#   PC1_2 = rep(NA, 26),
#   PC2_2 = rep(NA, 26),
#   PC3_2 = rep(NA, 26)
# )
#
# write.csv(formulaMatrix, file="../data/formulaMatrix.csv")

# Estimate realized species niches #####

replicates2 <- paste0("rep", 1:1000)
NPresences <- c(64)



estSp <- vector("list", length(species))
names(estSp) <- species
for(s in 1:length(species)){
  estSp[[s]] <- vector("list", 1)
  names(estSp[[s]]) <- "size64"
  for(n in 1:length(sizes)) {
    estSp[[s]] <- vector("list", length(replicates2))
    names(estSp[[s]]) <- replicates2

  }}




for (s in 1:length(species)) {
  switch(species[s],
         "broad_avg"={
           estSp[[s]]$L <-dnorm(south[,3], mu_avg[1],
                               sd=sd_broad[1,1]) * dnorm(south[,5],
                                                         mu_avg[2], sd=sd_broad[2,2])
         },
         "broad_ext" = {
           estSp[[s]]$L <-dnorm(south[,3], mu_ext[1],
                               sd=sd_broad[1,1]) * dnorm(south[,5],
                                                         mu_ext[2], sd=sd_broad[2,2])

         },
         "narrow_avg" = {
           estSp[[s]]$L <-dnorm(south[,3],
                               mu_avg[1], sd=sd_narrow[1,1]) *
             dnorm(south[,5], mu_avg[2], sd=sd_narrow[2,2])

         },
         "narrow_ext" = {
           estSp[[s]]$L <-dnorm(south[,3],
                               mu_ext[1],sd=sd_narrow[1,1]) *
             dnorm(south[,5], mu_ext[2], sd=sd_narrow[2,2])

         })#end of switch



    for(r in 1:length(replicates2)){
      estSp[[s]][[r]]$index <- matrix(ncol=64)
      estSp[[s]][[r]]$YSim <- matrix(0, nrow=NSites, ncol=1)
      estSp[[s]][[r]]$YSim  <- data.frame(estSp[[s]][[r]]$YSim )
      names(estSp[[s]][[r]]$YSim ) <- "SimSp"

      estSp[[s]][[r]]$index <- sample(1:NSites, 64, replace=F, prob=estSp[[s]]$L/sum(estSp[[s]]$L))

      for(site in 1:NSites) {
        if(site %in% estSp[[s]][[r]]$index) {
          estSp[[s]][[r]]$YSim[site,] <- 1
        }
      }









    }
}


for (s in 1:length(species)) {


    X_bar2 <- data.frame(mu_PC1 = rep(0,length(replicates2)+1),
                         sd_PC1 = rep(0, length(replicates2)+1),
                         mu_PC2 = rep(0, length(replicates2)+1),
                         sd_PC2 = rep(0, length(replicates2)+1),
                         mu_PC3 = rep(0, length(replicates2)+1),
                         sd_PC3 = rep(0, length(replicates2)+1),
                         N = rep(NA, length(replicates2)+1))


    rownames(X_bar2) <- c(paste0("sp", 1:length(replicates2)), "Mean")

    for (j in 1:length(replicates2)) {
      temp <- cbind(south, estSp[[s]][[j]]$YSim)
      temp <-subset(temp, temp$SimSp==1)

      X_bar2$mu_PC1[j] <- mean(temp$PC1) #mu_PC1
      X_bar2$sd_PC1[j] <- sd(temp$PC1) #sd_PC1
      X_bar2$min_PC1[j] <- min(temp$PC1)
      X_bar2$max_PC1[j] <- max(temp$PC1)

      X_bar2$mu_PC2[j] <- mean(temp$PC2) #mu_PC2
      X_bar2$sd_PC2[j] <- sd(temp$PC2) #sd_PC2
      X_bar2$min_PC2[j] <- min(temp$PC2)
      X_bar2$max_PC2[j] <- max(temp$PC2)

      X_bar2$mu_PC3[j] <- mean(temp$PC3) #mu_PC3
      X_bar2$sd_PC3[j] <- sd(temp$PC3) #sd_PC3
      X_bar2$min_PC3[j] <- min(temp$PC3)
      X_bar2$max_PC3[j] <- max(temp$PC3)
      X_bar2$N[j] <- length(temp$PC1)

    }

    X_bar2$mu_PC1[length(replicates2)+1] <- mean(X_bar2$mu_PC1[1:length(replicates2)])
    X_bar2$sd_PC1[length(replicates2)+1] <- mean(X_bar2$sd_PC1[1:length(replicates2)], na.rm=T)
    X_bar2$mu_PC2 <- mean(X_bar2$mu_PC2[1:length(replicates2)])
    X_bar2$sd_PC2[length(replicates2)+1] <- mean(X_bar2$sd_PC2[1:length(replicates2)], na.rm=T)
    X_bar2$mu_PC3[length(replicates2)+1] <- mean(X_bar2$mu_PC3[1:length(replicates2)])
    X_bar2$sd_PC3[length(replicates2)+1] <- mean(X_bar2$sd_PC3[1:length(replicates2)], na.rm=T)
    estSp[[s]]$X_bar2 <- X_bar2
  }


save(estSp, file="../data/estSp.RData")
rm(X_bar2)


