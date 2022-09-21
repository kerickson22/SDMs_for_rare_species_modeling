#Evaluate calibration for PC2 for all models
library(arm)
modelType <- c("HMSC_single_simple", "HMSC_joint_simple", "ESM", "ESM_bivariate")


memory.limit(2^32)
load("../data/south.RData")
NSites <- length(south$long)
row.names(south) <- 1:NSites
replicates <- paste0("rep", 1:30)
NPresences <- c(64, 32, 16, 8, 4, 2)
sizes <- paste0("size", NPresences)# Desired sample size
species <- c("wide_avg", "wide_ext",
             "narrow_avg", "narrow_ext")
splits <- paste0("split", 1:3)

sites <- 1:NSites



load("../data/X_bar.RData")
mu_ext <- c(-3.5, 2)
Sigma_narrow <- diag(c(0.9, 2.5))

mu_avg <- c(X_bar$mu_PC1[64], X_bar$mu_PC3[64])
Sigma_wide <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.90, na.rm=T)^2,  quantile(X_bar$sd_PC3, probs=0.90, na.rm=T)^2))

#Load simulated data
load(paste0("../data/sims_data.RData"))

x1 <- mu_avg[1]
x2 <- south$PC2
x3 <-mu_avg[2]

response <- expand.grid(x1, x2, x3)
names(response) <- c("PC1", "PC2", "PC3")








#Joint Simple
load("../data/models/HMSC_joint_simple/results2.RData")
for (i in 1:length(results$species)) {
  for ( j in 1:length(response$PC1)) {
    response$L[j] <- results$b1_1[i]*response$PC1[j] + results$b1_2[i]*response$PC1[j]*response$PC1[j] +
      results$b2_1[i]*response$PC2[j] + results$b2_2[i]*response$PC2[j]*response$PC2[j] +
      results$b3_1[i]*response$PC3[j] + results$b3_2[i]*response$PC3[j]*response$PC3[j]
     }
 results$varPC2Response[i]<- sd(invlogit(response$L))
}

save(results, file="../data/models/HMSC_joint_simple/results3.RData")
rm(results)

#Single Simple
load("../data/models/HMSC_single_simple/results3.RData")
for (i in 1:length(results$species)) {
  for ( j in 1:length(response$PC1)) {
    response$L[j] <- results$b1_1[i]*response$PC1[j] + results$b1_2[i]*response$PC1[j]*response$PC1[j] +
      results$b2_1[i]*response$PC2[j] + results$b2_2[i]*response$PC2[j]*response$PC2[j] +
      results$b3_1[i]*response$PC3[j] + results$b3_2[i]*response$PC3[j]*response$PC3[j]
  }
  results$varPC2Response[i]<- sd(invlogit(response$L))
}
save(results, file="../data/models/HMSC_single_simple/results3.RData")
rm(results)

#ESM
load("../data/models/ESM/results2.RData")
for (i in 1:length(results$species)) {
  for ( j in 1:length(response$PC1)) {
    response$L[j] <- results$b1_1[i]*response$PC1[j] + results$b1_2[i]*response$PC1[j]*response$PC1[j] +
      results$b2_1[i]*response$PC2[j] + results$b2_2[i]*response$PC2[j]*response$PC2[j] +
      results$b3_1[i]*response$PC3[j] + results$b3_2[i]*response$PC3[j]*response$PC3[j]
  }
  results$sdPC2Response[i]<- sd(invlogit(response$L))
}
save(results, file="../data/models/ESM/results3.RData")

rm(results)

#ESM_bivariate
load("../data/models/ESM_bivariate/results2.RData")
for (i in 1:length(results$species)) {
  for ( j in 1:length(response$PC1)) {
    response$L[j] <- results$b1_1[i]*response$PC1[j] +
      results$b2_1[i]*response$PC2[j]  +
      results$b3_1[i]*response$PC3[j]
  }
  results$sdPC2Response[i]<- sd(invlogit(response$L))
}
save(results, file="../data/models/ESM_bivariate/results3.RData")
rm(results)
