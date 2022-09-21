#Ensembles of Small Models
# K. D. Erickson
# 11/24/21

library(enmSdm)
library(omnibus)
library(PRROC)
library(statisfactory)
#remotes::install_github('adamlilith/statisfactory', dependencies=TRUE)
# memory.limit(2^32) #on Windows machine
load("../data/south.RData")
NSites <- length(south$long)
row.names(south) <- 1:NSites
replicates <- paste0("rep", 1:30)
NPresences <- c(64, 32, 16, 8, 4, 2)
sizes <- paste0("size", NPresences)# Desired sample size
species <- c("wide_avg", "wide_ext",
             "narrow_avg", "narrow_ext")
splits <- paste0("split", 1:3)


modelType  <- c("ESM", "ESM_bivariate")


sites <- 1:NSites

x1 <- seq(from=-9, to = 13, length.out=100)
x2 <- seq(from=-9, to = 16, length.out=100)
x3 <- seq(from=-11, to = 4, length.out=100)

load("../data/X_bar.RData")
mu_ext <- c(-3.5, 2)
Sigma_narrow <- diag(c(0.9, 2.5))

mu_avg <- c(X_bar$mu_PC1[64], X_bar$mu_PC3[64])
Sigma_wide <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.90, na.rm=T)^2,  quantile(X_bar$sd_PC3, probs=0.90, na.rm=T)^2))

#Load simulated data
load(paste0("../data/sims_data.RData"))


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


formulaMatrix <- read.csv("../data/formulaMatrix.csv",
                          header=T)

computeRMSEWeighted = function(Y, predY) {
  RMSE <- sqrt((mean((Y[Y==1]-predY[Y==1])^2, na.rm=TRUE) +
                  mean((Y[Y==0]-predY[Y==0])^2, na.rm=TRUE)))
  return(RMSE)
}


computeTjurR2 = function(Y, predY) {

    R2 = mean(predY[which(Y == 1)]) - mean(predY[which(Y == 0)])

  return(R2)
}

weights_df <- expand.grid(species, sizes, replicates)
names(weights_df) <- c("species", "sizes", "replicates")
for(i in 1:length(formulaMatrix$X)) {
weights_df[,3+i] <- rep(0, length(weights_df$species))
}
names(weights_df)[4:29] <- paste0("m", 1:26)

#ESM (with full model)
timeStart <- Sys.time()

for(s in 1:length(species)) {
  for(n in 1:length(sizes)) {
    for(r in 1:length(replicates)) {
      model <-list()
      model$seed <- .GlobalEnv$.Random.seed
      model$NPresences <- sizes[n]
      model$species <- species[s]
      model$replicate <- replicates[r]
      model$modelType <- modelType[1]

      switch(as.character(model$species),
             "wide_avg"={
               mu <-mu_avg
               Sigma <- Sigma_wide
             },
             "wide_ext" = {
               mu <-mu_ext
               Sigma <- Sigma_wide
             },
             "narrow_avg" = {
               mu <-mu_avg
               Sigma <- Sigma_narrow
             },
             "narrow_ext" = {
               mu <-mu_ext
               Sigma <- Sigma_narrow
             })#end of switch

      real_curve_PC1 <- dnorm(x1, mean=mu[1],
                              sd=Sigma[1,1])
      real_curve_PC3 <- dnorm(x3, mean=mu[2],
                              sd=Sigma[2,2])


      Y <- cbind(south[,6:68],as.matrix(sims[[s]][[n]][[r]]$YSim))
      XData <- south[,3:5]


      partition1 <- sims[[s]][[n]][[r]]$split1
      partition2 <- sims[[s]][[n]][[r]]$split2
      partition3 <- sims[[s]][[n]][[r]]$split3

      partitions <- vector(mode = "list", length = 3)
      partitions[[1]][[1]] <- cbind(Y=Y[partition1==1,]$SimSp, XData[partition1==1,])
      partitions[[1]][[2]] <- cbind(Y=Y[partition1==2,]$SimSp, XData[partition1==2,])
      partitions[[2]][[1]] <- cbind(Y=Y[partition2==1,]$SimSp, XData[partition2==1,])
      partitions[[2]][[2]] <- cbind(Y=Y[partition2==2,]$SimSp, XData[partition2==2,])
      partitions[[3]][[1]] <- cbind(Y=Y[partition3==1,]$SimSp, XData[partition3==1,])
      partitions[[3]][[2]] <- cbind(Y=Y[partition3==2,]$SimSp, XData[partition3==2,])


      #Information required for Weighted AUC calculation
      NPresent <- sum(Y[,64]/2)
      NAbsent <- (dim(Y)[1]/2)-NPresent
      presWeights <- rep(1/NPresent, NPresent)
      contrastWeights <- rep(1/NAbsent, NAbsent)

      weight_vec1 <- vector("list", 3)
      weight_vec2 <- vector("list", 3)
      pred_sum1   <- vector("list", 3)
      pred_sum2   <- vector("list", 3)
      pred_esm1   <- vector("list", 3)
      pred_esm2   <- vector("list", 3)
      B0_sum1     <- vector("list", 3)
      B1_sum1     <- vector("list", 3)
      B2_sum1     <- vector("list", 3)
      B3_sum1     <- vector("list", 3)
      B4_sum1     <- vector("list", 3)
      B5_sum1     <- vector("list", 3)
      B6_sum1     <- vector("list", 3)

      B0_sum2     <- vector("list", 3)
      B1_sum2     <- vector("list", 3)
      B2_sum2     <- vector("list", 3)
      B3_sum2     <- vector("list", 3)
      B4_sum2     <- vector("list", 3)
      B5_sum2     <- vector("list", 3)
      B6_sum2     <- vector("list", 3)

      B0_esm1     <- vector("list", 3)
      B1_esm1     <- vector("list", 3)
      B2_esm1     <- vector("list", 3)
      B3_esm1     <- vector("list", 3)
      B4_esm1     <- vector("list", 3)
      B5_esm1     <- vector("list", 3)
      B6_esm1     <- vector("list", 3)

      B0_esm2     <- vector("list", 3)
      B1_esm2     <- vector("list", 3)
      B2_esm2     <- vector("list", 3)
      B3_esm2     <- vector("list", 3)
      B4_esm2     <- vector("list", 3)
      B5_esm2     <- vector("list", 3)
      B6_esm2     <- vector("list", 3)

      aucWeighted_ESM1 <- vector("list", 3)
      aucWeighted_ESM2 <- vector("list", 3)

      RMSEWeighted_ESM1 <- vector("list", 3)
      RMSEWeighted_ESM2 <- vector("list", 3)

      PR_ESM1 <- vector("list", 3)
      PR_ESM2 <- vector("list", 3)

      TjurR2_ESM1 <- vector("list", 3)
      TjurR2_ESM2 <- vector("list", 3)

      model_curve_PC1_ESM1 <- vector("list", 3)
      model_curve_PC1_ESM2 <- vector("list", 3)

      model_curve_PC2_ESM1 <- vector("list", 3)
      model_curve_PC2_ESM2 <- vector("list", 3)

      model_curve_PC3_ESM1 <- vector("list", 3)
      model_curve_PC3_ESM2 <- vector("list", 3)

      PC1_rankCor_ESM1 <- vector("list", 3)
      PC1_rankCor_ESM2 <- vector("list", 3)

      PC2_sd_ESM1 <- vector("list", 3)
      PC2_sd_ESM2 <- vector("list", 3)

      PC3_rankCor_ESM1 <- vector("list", 3)
      PC3_rankCor_ESM2 <- vector("list", 3)

      for(k in 1:3){
        train1 <- partitions[[k]][[1]]
        test1  <- partitions[[k]][[2]]
        train2 <- test1
        test2  <- train1

        weight_vec1[[k]] <- vector()
        weight_vec2[[k]] <- vector() # vector with weights of single bivariate models

        pred_sum1[[k]] <- rep(0, 0.5*dim(Y)[1])
        pred_sum2[[k]] <- rep(0, 0.5*dim(Y)[1]) # for ensemble prediction of bivariate models

        if(model$NPresences == "size2") {
          for(f in 1:3) {

            #A: Train model on Fold 1, test model on Fold 2
            mod1 <- glm(as.formula(formulaMatrix$formula[f]),
                        family=binomial, data=train1)
            X1 <- test1[,2:4]
            newX <- as.data.frame(X1[,which(!is.na(formulaMatrix[f, 4:6]))])
            names(newX) <- names(X1[which(!is.na(formulaMatrix[f, 4:6]))])
            pred1 <- predict(mod1, newX,
                             type="response")

            #Evaluate bivariate model and calculate its weight
            aucWeighted1 <- aucWeighted(pres=pred1[test1$Y==1],
                                        contrast=pred1[test1$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)
            weight1 <- 2*aucWeighted1 -1 #Schoners' D
            if(weight1 < 0) {weight1 <- 0}
            weight_vec1[[k]] <- c(weight_vec1[[k]], weight1)
            pred_sum1[[k]] <- pred_sum1[[k]] + (pred1 * weight1)
            if(!is.na(formulaMatrix$Int[f]))   {B0_sum1[[k]] <- sum(B0_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$Int[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1[f]))   {B1_sum1[[k]] <- sum(B1_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2[f]))   {B2_sum1[[k]] <- sum(B2_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3[f]))   {B3_sum1[[k]] <- sum(B3_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum1[[k]] <- sum(B4_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum1[[k]] <- sum(B5_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum1[[k]] <- sum(B6_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3_2[f]], na.rm=T)}

            #B: Train model on Fold 2, test model on Fold 1
            mod2 <- glm(as.formula(formulaMatrix$formula[f]),
                        family=binomial, data=train2)
            X2 <- test2[,2:4]
            newX <- as.data.frame(X2[,which(!is.na(formulaMatrix[f, 4:6]))])
            names(newX) <- names(X2[which(!is.na(formulaMatrix[f, 4:6]))])
            pred2 <- predict(mod2, newX,
                             type="response")

            #Evaluate bivariate model and calculate its weight
            aucWeighted2 <- aucWeighted(pres=pred2[test2$Y==1],
                                        contrast=pred2[test2$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)
            weight2 <- 2*aucWeighted2 -1 #Schoners' D
            if(weight2 < 0) {weight2 <- 0}
            weight_vec2[[k]] <- c(weight_vec2[[k]], weight2)
            pred_sum2[[k]] <- pred_sum2[[k]] + (pred2 * weight2)
            if(!is.na(formulaMatrix$Int[f]))   {B0_sum2[[k]] <- sum(B0_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$Int[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1[f]))   {B1_sum2[[k]] <- sum(B1_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2[f]))   {B2_sum2[[k]] <- sum(B2_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3[f]))   {B3_sum2[[k]] <- sum(B3_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum2[[k]] <- sum(B4_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum2[[k]] <- sum(B5_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum2[[k]] <- sum(B6_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3_2[f]], na.rm=T)}
          }
        }
        if(model$NPresences == "size4") {
          for(f in 1:9) {

            #A: Train model on Fold 1, test model on Fold 2
            mod1 <- glm(as.formula(formulaMatrix$formula[f]),
                        family=binomial, data=train1)
            X1 <- test1[,2:4]
            newX <- as.data.frame(X1[,which(!is.na(formulaMatrix[f, 4:6]))])
            names(newX) <- names(X1[which(!is.na(formulaMatrix[f, 4:6]))])
            pred1 <- predict(mod1, newX,
                             type="response")

            #Evaluate bivariate model and calculate its weight
            aucWeighted1 <- aucWeighted(pres=pred1[test1$Y==1],
                                        contrast=pred1[test1$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)
            weight1 <- 2*aucWeighted1 -1 #Schoners' D
            if(weight1 < 0) {weight1 <- 0}
            weight_vec1[[k]] <- c(weight_vec1[[k]], weight1)
            pred_sum1[[k]] <- pred_sum1[[k]] + (pred1 * weight1)
            if(!is.na(formulaMatrix$Int[f]))   {B0_sum1[[k]] <- sum(B0_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$Int[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1[f]))   {B1_sum1[[k]] <- sum(B1_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2[f]))   {B2_sum1[[k]] <- sum(B2_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3[f]))   {B3_sum1[[k]] <- sum(B3_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum1[[k]] <- sum(B4_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum1[[k]] <- sum(B5_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum1[[k]] <- sum(B6_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3_2[f]], na.rm=T)}

            #B: Train model on Fold 2, test model on Fold 1
            mod2 <- glm(as.formula(formulaMatrix$formula[f]),
                        family=binomial, data=train2)
            X2 <- test2[,2:4]
            newX <- as.data.frame(X2[,which(!is.na(formulaMatrix[f, 4:6]))])
            names(newX) <- names(X2[which(!is.na(formulaMatrix[f, 4:6]))])
            pred2 <- predict(mod2, newX,
                             type="response")

            #Evaluate bivariate model and calculate its weight
            aucWeighted2 <- aucWeighted(pres=pred2[test2$Y==1],
                                        contrast=pred2[test2$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)
            weight2 <- 2*aucWeighted2 -1 #Schoners' D
            if(weight2 < 0) {weight2 <- 0}
            weight_vec2[[k]] <- c(weight_vec2[[k]], weight2)
            pred_sum2[[k]] <- pred_sum2[[k]] + (pred2 * weight2)
            if(!is.na(formulaMatrix$Int[f]))   {B0_sum2[[k]] <- sum(B0_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$Int[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1[f]))   {B1_sum2[[k]] <- sum(B1_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2[f]))   {B2_sum2[[k]] <- sum(B2_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3[f]))   {B3_sum2[[k]] <- sum(B3_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum2[[k]] <- sum(B4_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum2[[k]] <- sum(B5_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum2[[k]] <- sum(B6_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3_2[f]], na.rm=T)}
          }
        }
        else{
        for(f in 1:length(formulaMatrix$formula)) {

          #A: Train model on Fold 1, test model on Fold 2
          mod1 <- glm(as.formula(formulaMatrix$formula[f]),
                      family=binomial, data=train1)
          X1 <- test1[,2:4]
          newX <- as.data.frame(X1[,which(!is.na(formulaMatrix[f, 4:6]))])
          names(newX) <- names(X1[which(!is.na(formulaMatrix[f, 4:6]))])
          pred1 <- predict(mod1, newX,
                           type="response")

          #Evaluate bivariate model and calculate its weight
          aucWeighted1 <- aucWeighted(pres=pred1[test1$Y==1],
                                      contrast=pred1[test1$Y==0],
                                      presWeight = presWeights,
                                      contrastWeight = contrastWeights)
          weight1 <- 2*aucWeighted1 -1 #Schoners' D
          if(weight1 < 0) {weight1 <- 0}
          weight_vec1[[k]] <- c(weight_vec1[[k]], weight1)
          pred_sum1[[k]] <- pred_sum1[[k]] + (pred1 * weight1)
          if(!is.na(formulaMatrix$Int[f]))   {B0_sum1[[k]] <- sum(B0_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$Int[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1[f]))   {B1_sum1[[k]] <- sum(B1_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2[f]))   {B2_sum1[[k]] <- sum(B2_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3[f]))   {B3_sum1[[k]] <- sum(B3_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum1[[k]] <- sum(B4_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum1[[k]] <- sum(B5_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum1[[k]] <- sum(B6_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3_2[f]], na.rm=T)}

          #B: Train model on Fold 2, test model on Fold 1
          mod2 <- glm(as.formula(formulaMatrix$formula[f]),
                      family=binomial, data=train2)
          X2 <- test2[,2:4]
          newX <- as.data.frame(X2[,which(!is.na(formulaMatrix[f, 4:6]))])
          names(newX) <- names(X2[which(!is.na(formulaMatrix[f, 4:6]))])
          pred2 <- predict(mod2, newX,
                           type="response")

          #Evaluate bivariate model and calculate its weight
          aucWeighted2 <- aucWeighted(pres=pred2[test2$Y==1],
                                      contrast=pred2[test2$Y==0],
                                      presWeight = presWeights,
                                      contrastWeight = contrastWeights)
          weight2 <- 2*aucWeighted2 -1 #Schoners' D
          if(weight2 < 0) {weight2 <- 0}
          weight_vec2[[k]] <- c(weight_vec2[[k]], weight2)
          pred_sum2[[k]] <- pred_sum2[[k]] + (pred2 * weight2)
          if(!is.na(formulaMatrix$Int[f]))   {B0_sum2[[k]] <- sum(B0_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$Int[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1[f]))   {B1_sum2[[k]] <- sum(B1_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2[f]))   {B2_sum2[[k]] <- sum(B2_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3[f]))   {B3_sum2[[k]] <- sum(B3_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum2[[k]] <- sum(B4_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum2[[k]] <- sum(B5_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum2[[k]] <- sum(B6_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3_2[f]], na.rm=T)}
        }
        }
        if(sum(weight_vec1[[k]])==0 |
           sum(weight_vec2[[k]]) ==0){
          pred_esm1[[k]] <- NA
          pred_esm2[[k]] <- NA

          B0_esm1[[k]] <- NA
          B1_esm1[[k]] <- NA
          B2_esm1[[k]] <- NA
          B3_esm1[[k]] <- NA
          B4_esm1[[k]] <- NA
          B5_esm1[[k]] <- NA
          B6_esm1[[k]] <- NA

          B0_esm2[[k]] <- NA
          B1_esm2[[k]] <- NA
          B2_esm2[[k]] <- NA
          B3_esm2[[k]] <- NA
          B4_esm2[[k]] <- NA
          B5_esm2[[k]] <- NA
          B6_esm2[[k]] <- NA

          aucWeighted_ESM1[[k]] <- NA
          aucWeighted_ESM2[[k]] <- NA
          RMSEWeighted_ESM1[[k]] <- NA
          RMSEWeighted_ESM2[[k]] <- NA
          PR_ESM1[[k]] <- NA
          PR_ESM2[[k]] <- NA
          TjurR2_ESM1[[k]] <- NA
          TjurR2_ESM2[[k]] <- NA
          model_curve_PC1_ESM1[[k]] <- NA
          model_curve_PC1_ESM2[[k]] <- NA
          model_curve_PC2_ESM1[[k]] <- NA
          model_curve_PC2_ESM2[[k]] <- NA
          model_curve_PC3_ESM1[[k]] <- NA
          model_curve_PC3_ESM2[[k]] <- NA

          PC1_rankCor_ESM1[[k]] <- NA
          PC1_rankCor_ESM2[[k]] <- NA

          PC2_sd_ESM1[[k]] <- NA
          PC2_sd_ESM2[[k]] <- NA

          PC3_rankCor_ESM1[[k]] <- NA
          PC3_rankCor_ESM2[[k]] <- NA
          }
        else{
        # calculate ESM prediction
        pred_esm1[[k]] <- pred_sum1[[k]]/sum(weight_vec1[[k]])
        pred_esm2[[k]] <- pred_sum2[[k]]/sum(weight_vec2[[k]])

        #Calculate betas
        B0_esm1[[k]] <- B0_sum1[[k]]/sum(weight_vec1[[k]])
        B1_esm1[[k]] <- B1_sum1[[k]]/sum(weight_vec1[[k]])
        B2_esm1[[k]] <- B2_sum1[[k]]/sum(weight_vec1[[k]])
        B3_esm1[[k]] <- B3_sum1[[k]]/sum(weight_vec1[[k]])
        B4_esm1[[k]] <- B4_sum1[[k]]/sum(weight_vec1[[k]])
        B5_esm1[[k]] <- B5_sum1[[k]]/sum(weight_vec1[[k]])
        B6_esm1[[k]] <- B6_sum1[[k]]/sum(weight_vec1[[k]])

        B0_esm2[[k]] <- B0_sum2[[k]]/sum(weight_vec2[[k]])
        B1_esm2[[k]] <- B1_sum2[[k]]/sum(weight_vec2[[k]])
        B2_esm2[[k]] <- B2_sum2[[k]]/sum(weight_vec2[[k]])
        B3_esm2[[k]] <- B3_sum2[[k]]/sum(weight_vec2[[k]])
        B4_esm2[[k]] <- B4_sum2[[k]]/sum(weight_vec2[[k]])
        B5_esm2[[k]] <- B5_sum2[[k]]/sum(weight_vec2[[k]])
        B6_esm2[[k]] <- B6_sum2[[k]]/sum(weight_vec2[[k]])

        #evaluate ESM
        aucWeighted_ESM1[[k]] <- aucWeighted(pres=pred_esm1[[k]][test1$Y==1],
                                        contrast=pred_esm1[[k]][test1$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)
        aucWeighted_ESM2[[k]] <- aucWeighted(pres=pred_esm2[[k]][test2$Y==1],
                                        contrast=pred_esm2[[k]][test2$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)

        RMSEWeighted_ESM1[[k]] <- computeRMSEWeighted(Y=test1$Y, predY=pred_esm1[[k]])
        RMSEWeighted_ESM2[[k]] <- computeRMSEWeighted(Y=test2$Y, predY=pred_esm2[[k]])

        PR_ESM1[[k]] <- pr.curve(scores.class0=pred_esm1[[k]],
                           weights.class0 = test1$Y,
                           curve=TRUE)$auc.integral
        PR_ESM2[[k]] <- pr.curve(scores.class0=pred_esm2[[k]],
                                 weights.class0 = test2$Y,
                                 curve=TRUE)$auc.integral

        TjurR2_ESM1[[k]] <- computeTjurR2(Y=test1$Y,
                                          predY=pred_esm1[[k]])
        TjurR2_ESM2[[k]] <- computeTjurR2(Y=test2,
                                          predY=pred_esm2[[k]])
       #Species Response Curves
        model_curve_PC1_ESM1[[k]] <- pnorm(B1_esm1[[k]]*x1 + B4_esm1[[k]]*x1*x1)
        model_curve_PC1_ESM2[[k]] <- pnorm(B1_esm2[[k]]*x1 + B4_esm2[[k]]*x1*x2)

        model_curve_PC2_ESM1[[k]] <- pnorm(B2_esm1[[k]]*x2 + B5_esm1[[k]]*x2*x2)
        model_curve_PC2_ESM2[[k]] <- pnorm(B2_esm2[[k]]*x2 + B5_esm2[[k]]*x2*x2)

        model_curve_PC3_ESM1[[k]] <- pnorm(B3_esm1[[k]]*x3 + B6_esm1[[k]]*x3*x3)
        model_curve_PC3_ESM2[[k]] <- pnorm(B3_esm2[[k]]*x3 + B6_esm2[[k]]*x3*x3)

        if(sd(model_curve_PC1_ESM1[[k]])>0) {PC1_rankCor_ESM1[[k]] <- compareResponse(model_curve_PC1_ESM1[[k]],
          real_curve_PC1,
          data =data.frame(x1),
          graph=F)$rankCor}
        if(sd(model_curve_PC1_ESM2[[k]])>0) {PC1_rankCor_ESM2[[k]] <- compareResponse(model_curve_PC1_ESM2[[k]],
          real_curve_PC1,
          data =data.frame(x1),
          graph=F)$rankCor}

        PC2_sd_ESM1[[k]] <- sd(model_curve_PC2_ESM1[[k]])
        PC2_sd_ESM2[[k]] <- sd(model_curve_PC2_ESM2[[k]])

        if(sd(model_curve_PC3_ESM1[[k]])>0) {PC3_rankCor_ESM1[[k]] <- compareResponse(model_curve_PC3_ESM1[[k]],
          real_curve_PC3,
          data =data.frame(x3),
          graph=F)$rankCor}
        if(sd(model_curve_PC3_ESM2[[k]])>0) {PC3_rankCor_ESM2[[k]] <- compareResponse(model_curve_PC3_ESM2[[k]],
          real_curve_PC3,
          data =data.frame(x3),
          graph=F)$rankCor}
        }
      }

      weights_df[i,4:29] <- (weight_vec1[[1]] + weight_vec1[[2]]+weight_vec1[[3]] +
                          weight_vec2[[1]] + weight_vec2[[2]]+weight_vec2[[3]])/6

      weight_vec2avg <- mean(unlist(lapply(weight_vec2, sapply, mean)))


      aucWeighted_ESM1avg <- mean(unlist(lapply(aucWeighted_ESM1,sapply,mean)))
      aucWeighted_ESM2avg <- mean(unlist(lapply(aucWeighted_ESM2, sapply, mean)))
      if(is.na(aucWeighted_ESM1avg) | is.na(aucWeighted_ESM2avg)) {model$aucWeighted <- NA}
      else{model$aucWeighted <- mean(aucWeighted_ESM1avg, aucWeighted_ESM2avg)}
      RMSEWeighted_ESM1avg <- mean(unlist(lapply(RMSEWeighted_ESM1,sapply,mean)))
      RMSEWeighted_ESM2avg <-   mean(unlist(lapply(RMSEWeighted_ESM2,sapply,mean)))
      if(is.na(RMSEWeighted_ESM1avg) | is.na(RMSEWeighted_ESM2avg)) {model$RMSEWeighted <- NA}
      else{model$RMSEWeighted <- mean(RMSEWeighted_ESM1avg,
                                      RMSEWeighted_ESM2avg)}
      PR_ESM1avg <- mean(unlist(lapply(PR_ESM1, sapply, mean)))
      PR_ESM2avg <- mean(unlist(lapply(PR_ESM2, sapply, mean)))
      if(is.na(PR_ESM1avg) | is.na(PR_ESM2avg)) {model$PR <- NA}
      else{model$PR <- mean(PR_ESM1avg,
                       PR_ESM2avg)}
      TjurR2_ESM1avg <- mean(unlist(lapply(TjurR2_ESM1, sapply, mean)))
      TjurR2_ESM2avg <- mean(unlist(lapply(TjurR2_ESM2, sapply, mean)))
      if(is.na(TjurR2_ESM1avg) | is.na(TjurR2_ESM2avg)) {model$TjurR2 <- NA}
      else{model$TjurR2 <- mean(TjurR2_ESM1avg,
                          TjurR2_ESM2avg)}

      PC1_rankCor_ESM1avg <- mean(unlist(lapply(PC1_rankCor_ESM1, sapply, mean)))
      PC1_rankCor_ESM2avg <- mean(unlist(lapply(PC1_rankCor_ESM2, sapply, mean)))
      if(is.na(PC1_rankCor_ESM1avg) | is.na(PC1_rankCor_ESM2avg)) {model$PC1_rankCor <- NA}
      else{model$PC1_rankCor <- mean(PC1_rankCor_ESM1avg, PC1_rankCor_ESM2avg)}

      PC2_sd_ESM1avg <- mean(unlist(lapply(PC2_sd_ESM1, sapply, mean)))
      PC2_sd_ESM2avg <- mean(unlist(lapply(PC2_sd_ESM2, sapply, mean)))
      if(is.na(PC2_sd_ESM1avg) | is.na(PC2_sd_ESM2avg)) {model$PC2_sd <- NA}
      else{model$PC2_sd <- mean(PC2_sd_ESM1avg,
                                PC2_sd_ESM2avg)}

      PC3_rankCor_ESM1avg <- mean(unlist(lapply(PC3_rankCor_ESM1, sapply, mean)))
      PC3_rankCor_ESM2avg <- mean(unlist(lapply(PC3_rankCor_ESM2, sapply, mean)))
      if(is.na(PC3_rankCor_ESM1avg) | is.na(PC3_rankCor_ESM2avg)) {model$PC3_rankCor <- NA}
      else{model$PC3_rankCor <- mean(PC3_rankCor_ESM1avg, PC3_rankCor_ESM2avg)}

      #Linear term for PC1
      b1_1_ESM1_avg <- mean(unlist(lapply(B1_esm1, sapply, mean)))
      b1_1_ESM2_avg <- mean(unlist(lapply(B1_esm2, sapply, mean)))
      if(is.na(b1_1_ESM1_avg) | is.na(b1_1_ESM2_avg)) {model$b1_1 <- NA}
      else{model$b1_1 <- mean(b1_1_ESM1_avg, b1_1_ESM2_avg)}

      #Quadratic term for PC1
      b1_2_ESM1_avg <- mean(unlist(lapply(B4_esm1, sapply, mean)))
      b1_2_ESM2_avg <- mean(unlist(lapply(B4_esm2, sapply, mean)))
      if(is.na(b1_2_ESM1_avg) | is.na(b1_2_ESM2_avg)) {model$b1_2 <- NA}
      else{model$b1_2 <- mean(b1_2_ESM1_avg, b1_2_ESM2_avg)}

      #Linear term for PC2
      b2_1_ESM1_avg <- mean(unlist(lapply(B2_esm1, sapply, mean)))
      b2_1_ESM2_avg <- mean(unlist(lapply(B2_esm2, sapply, mean)))
      if(is.na(b2_1_ESM1_avg) | is.na(b2_1_ESM2_avg)) {model$b2_1 <- NA}
      else{model$b2_1 <- mean(b2_1_ESM1_avg, b2_1_ESM2_avg)}

      #Quadratic term for PC2
      b2_2_ESM1_avg <- mean(unlist(lapply(B5_esm1, sapply, mean)))
      b2_2_ESM2_avg <- mean(unlist(lapply(B5_esm2, sapply, mean)))
      if(is.na(b2_2_ESM1_avg) | is.na(b2_2_ESM2_avg)) {model$b2_2 <- NA}
      else{model$b2_2 <- mean(b2_2_ESM1_avg, b2_2_ESM2_avg)}

      #Linear term for PC3
      b3_1_ESM1_avg <- mean(unlist(lapply(B3_esm1, sapply, mean)))
      b3_1_ESM2_avg <- mean(unlist(lapply(B3_esm2, sapply, mean)))
      if(is.na(b3_1_ESM1_avg) | is.na(b3_1_ESM2_avg)) {model$b3_1 <- NA}
      else{model$b3_1 <- mean(b3_1_ESM1_avg, b3_1_ESM2_avg)}

      #Quadratic term for PC3
      b3_2_ESM1_avg <- mean(unlist(lapply(B6_esm1, sapply, mean)))
      b3_2_ESM2_avg <- mean(unlist(lapply(B6_esm2, sapply, mean)))
      if(is.na(b3_2_ESM1_avg) | is.na(b3_2_ESM2_avg)) {model$b3_2 <- NA}
      else{model$b3_2 <- mean(b3_2_ESM1_avg, b3_2_ESM2_avg)}


      save(model, file=paste0("../data/models/",
                              modelType[1], "/",species[s], "/",
                              sizes[n], "/", "model_",replicates[r],
                              ".RData"))
      cat(species[s], sizes[n], replicates[r], "\n")

     i <- i + 1}
  }
}
runTime <- Sys.time() - timeStart #1.09 mins

#ESM_bivariate
timeStart <- Sys.time()
for(s in 1:length(species)) {
  for(n in 1:length(sizes)) {
    for(r in 1:length(replicates)) {
      model <-list()
      model$seed <- .GlobalEnv$.Random.seed
      model$NPresences <- sizes[n]
      model$species <- species[s]
      model$replicate <- replicates[r]
      model$modelType <- modelType[2]

      switch(as.character(model$species),
             "wide_avg"={
               mu <-mu_avg
               Sigma <- Sigma_wide
             },
             "wide_ext" = {
               mu <-mu_ext
               Sigma <- Sigma_wide
             },
             "narrow_avg" = {
               mu <-mu_avg
               Sigma <- Sigma_narrow
             },
             "narrow_ext" = {
               mu <-mu_ext
               Sigma <- Sigma_narrow
             })#end of switch

      real_curve_PC1 <- dnorm(x1, mean=mu[1],
                              sd=Sigma[1,1])
      real_curve_PC3 <- dnorm(x3, mean=mu[2],
                              sd=Sigma[2,2])


      Y <- cbind(south[,6:68],as.matrix(sims[[s]][[n]][[r]]$YSim))
      XData <- south[,3:5]


      partition1 <- sims[[s]][[n]][[r]]$split1
      partition2 <- sims[[s]][[n]][[r]]$split2
      partition3 <- sims[[s]][[n]][[r]]$split3

      partitions <- vector(mode = "list", length = 3)
      partitions[[1]][[1]] <- cbind(Y=Y[partition1==1,]$SimSp, XData[partition1==1,])
      partitions[[1]][[2]] <- cbind(Y=Y[partition1==2,]$SimSp, XData[partition1==2,])
      partitions[[2]][[1]] <- cbind(Y=Y[partition2==1,]$SimSp, XData[partition2==1,])
      partitions[[2]][[2]] <- cbind(Y=Y[partition2==2,]$SimSp, XData[partition2==2,])
      partitions[[3]][[1]] <- cbind(Y=Y[partition3==1,]$SimSp, XData[partition3==1,])
      partitions[[3]][[2]] <- cbind(Y=Y[partition3==2,]$SimSp, XData[partition3==2,])


      #Information required for Weighted AUC calculation
      NPresent <- sum(Y[,64]/2)
      NAbsent <- (dim(Y)[1]/2)-NPresent
      presWeights <- rep(1/NPresent, NPresent)
      contrastWeights <- rep(1/NAbsent, NAbsent)

      weight_vec1 <- vector("list", 3)
      weight_vec2 <- vector("list", 3)
      pred_sum1   <- vector("list", 3)
      pred_sum2   <- vector("list", 3)
      pred_esm1   <- vector("list", 3)
      pred_esm2   <- vector("list", 3)
      B0_sum1     <- vector("list", 3)
      B1_sum1     <- vector("list", 3)
      B2_sum1     <- vector("list", 3)
      B3_sum1     <- vector("list", 3)
      B4_sum1     <- vector("list", 3)
      B5_sum1     <- vector("list", 3)
      B6_sum1     <- vector("list", 3)

      B0_sum2     <- vector("list", 3)
      B1_sum2     <- vector("list", 3)
      B2_sum2     <- vector("list", 3)
      B3_sum2     <- vector("list", 3)
      B4_sum2     <- vector("list", 3)
      B5_sum2     <- vector("list", 3)
      B6_sum2     <- vector("list", 3)

      B0_esm1     <- vector("list", 3)
      B1_esm1     <- vector("list", 3)
      B2_esm1     <- vector("list", 3)
      B3_esm1     <- vector("list", 3)
      B4_esm1     <- vector("list", 3)
      B5_esm1     <- vector("list", 3)
      B6_esm1     <- vector("list", 3)

      B0_esm2     <- vector("list", 3)
      B1_esm2     <- vector("list", 3)
      B2_esm2     <- vector("list", 3)
      B3_esm2     <- vector("list", 3)
      B4_esm2     <- vector("list", 3)
      B5_esm2     <- vector("list", 3)
      B6_esm2     <- vector("list", 3)

      aucWeighted_ESM1 <- vector("list", 3)
      aucWeighted_ESM2 <- vector("list", 3)

      RMSEWeighted_ESM1 <- vector("list", 3)
      RMSEWeighted_ESM2 <- vector("list", 3)

      PR_ESM1 <- vector("list", 3)
      PR_ESM2 <- vector("list", 3)

      TjurR2_ESM1 <- vector("list", 3)
      TjurR2_ESM2 <- vector("list", 3)

      model_curve_PC1_ESM1 <- vector("list", 3)
      model_curve_PC1_ESM2 <- vector("list", 3)

      model_curve_PC2_ESM1 <- vector("list", 3)
      model_curve_PC2_ESM2 <- vector("list", 3)

      model_curve_PC3_ESM1 <- vector("list", 3)
      model_curve_PC3_ESM2 <- vector("list", 3)

      PC1_rankCor_ESM1 <- vector("list", 3)
      PC1_rankCor_ESM2 <- vector("list", 3)

      PC2_sd_ESM1 <- vector("list", 3)
      PC2_sd_ESM2 <- vector("list", 3)

      PC3_rankCor_ESM1 <- vector("list", 3)
      PC3_rankCor_ESM2 <- vector("list", 3)

      for(k in 1:3){
        train1 <- partitions[[k]][[1]]
        test1  <- partitions[[k]][[2]]
        train2 <- test1
        test2  <- train1

        weight_vec1[[k]] <- vector()
        weight_vec2[[k]] <- vector() # vector with weights of single bivariate models

        pred_sum1[[k]] <- rep(0, 0.5*dim(Y)[1])
        pred_sum2[[k]] <- rep(0, 0.5*dim(Y)[1]) # for ensemble prediction of bivariate models

        if(model$NPresences == "size2") {
          for(f in 1:3) {

            #A: Train model on Fold 1, test model on Fold 2
            mod1 <- glm(as.formula(formulaMatrix$formula[f]),
                        family=binomial, data=train1)
            X1 <- test1[,2:4]
            newX <- as.data.frame(X1[,which(!is.na(formulaMatrix[f, 4:6]))])
            names(newX) <- names(X1[which(!is.na(formulaMatrix[f, 4:6]))])
            pred1 <- predict(mod1, newX,
                             type="response")

            #Evaluate bivariate model and calculate its weight
            aucWeighted1 <- aucWeighted(pres=pred1[test1$Y==1],
                                        contrast=pred1[test1$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)
            weight1 <- 2*aucWeighted1 -1 #Schoners' D
            if(weight1 < 0) {weight1 <- 0}
            weight_vec1[[k]] <- c(weight_vec1[[k]], weight1)
            pred_sum1[[k]] <- pred_sum1[[k]] + (pred1 * weight1)
            if(!is.na(formulaMatrix$Int[f]))   {B0_sum1[[k]] <- sum(B0_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$Int[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1[f]))   {B1_sum1[[k]] <- sum(B1_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2[f]))   {B2_sum1[[k]] <- sum(B2_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3[f]))   {B3_sum1[[k]] <- sum(B3_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum1[[k]] <- sum(B4_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum1[[k]] <- sum(B5_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum1[[k]] <- sum(B6_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3_2[f]], na.rm=T)}

            #B: Train model on Fold 2, test model on Fold 1
            mod2 <- glm(as.formula(formulaMatrix$formula[f]),
                        family=binomial, data=train2)
            X2 <- test2[,2:4]
            newX <- as.data.frame(X2[,which(!is.na(formulaMatrix[f, 4:6]))])
            names(newX) <- names(X2[which(!is.na(formulaMatrix[f, 4:6]))])
            pred2 <- predict(mod2, newX,
                             type="response")

            #Evaluate bivariate model and calculate its weight
            aucWeighted2 <- aucWeighted(pres=pred2[test2$Y==1],
                                        contrast=pred2[test2$Y==0],
                                        presWeight = presWeights,
                                        contrastWeight = contrastWeights)
            weight2 <- 2*aucWeighted2 -1 #Schoners' D
            if(weight2 < 0) {weight2 <- 0}
            weight_vec2[[k]] <- c(weight_vec2[[k]], weight2)
            pred_sum2[[k]] <- pred_sum2[[k]] + (pred2 * weight2)
            if(!is.na(formulaMatrix$Int[f]))   {B0_sum2[[k]] <- sum(B0_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$Int[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1[f]))   {B1_sum2[[k]] <- sum(B1_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2[f]))   {B2_sum2[[k]] <- sum(B2_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3[f]))   {B3_sum2[[k]] <- sum(B3_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum2[[k]] <- sum(B4_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum2[[k]] <- sum(B5_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2_2[f]], na.rm=T)}
            if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum2[[k]] <- sum(B6_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3_2[f]], na.rm=T)}
          }
        }
        else{
        for(f in 1:6) {

          #A: Train model on Fold 1, test model on Fold 2
          mod1 <- glm(as.formula(formulaMatrix$formula[f]),
                      family=binomial, data=train1)
          X1 <- test1[,2:4]
          newX <- as.data.frame(X1[,which(!is.na(formulaMatrix[f, 4:6]))])
          names(newX) <- names(X1[which(!is.na(formulaMatrix[f, 4:6]))])
          pred1 <- predict(mod1, newX,
                           type="response")

          #Evaluate bivariate model and calculate its weight
          aucWeighted1 <- aucWeighted(pres=pred1[test1$Y==1],
                                      contrast=pred1[test1$Y==0],
                                      presWeight = presWeights,
                                      contrastWeight = contrastWeights)
          weight1 <- 2*aucWeighted1 -1 #Schoners' D
          if(weight1 < 0) {weight1 <- 0}
          weight_vec1[[k]] <- c(weight_vec1[[k]], weight1)
          pred_sum1[[k]] <- pred_sum1[[k]] + (pred1 * weight1)
          if(!is.na(formulaMatrix$Int[f]))   {B0_sum1[[k]] <- sum(B0_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$Int[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1[f]))   {B1_sum1[[k]] <- sum(B1_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2[f]))   {B2_sum1[[k]] <- sum(B2_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3[f]))   {B3_sum1[[k]] <- sum(B3_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum1[[k]] <- sum(B4_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC1_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum1[[k]] <- sum(B5_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC2_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum1[[k]] <- sum(B6_sum1[[k]], weight1*coefficients(mod1)[formulaMatrix$PC3_2[f]], na.rm=T)}

          #B: Train model on Fold 2, test model on Fold 1
          mod2 <- glm(as.formula(formulaMatrix$formula[f]),
                      family=binomial, data=train2)
          X2 <- test2[,2:4]
          newX <- as.data.frame(X2[,which(!is.na(formulaMatrix[f, 4:6]))])
          names(newX) <- names(X2[which(!is.na(formulaMatrix[f, 4:6]))])
          pred2 <- predict(mod2, newX,
                           type="response")

          #Evaluate bivariate model and calculate its weight
          aucWeighted2 <- aucWeighted(pres=pred2[test2$Y==1],
                                      contrast=pred2[test2$Y==0],
                                      presWeight = presWeights,
                                      contrastWeight = contrastWeights)
          weight2 <- 2*aucWeighted2 -1 #Schoners' D
          if(weight2 < 0) {weight2 <- 0}
          weight_vec2[[k]] <- c(weight_vec2[[k]], weight2)
          pred_sum2[[k]] <- pred_sum2[[k]] + (pred2 * weight2)
          if(!is.na(formulaMatrix$Int[f]))   {B0_sum2[[k]] <- sum(B0_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$Int[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1[f]))   {B1_sum2[[k]] <- sum(B1_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2[f]))   {B2_sum2[[k]] <- sum(B2_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3[f]))   {B3_sum2[[k]] <- sum(B3_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum2[[k]] <- sum(B4_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC1_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum2[[k]] <- sum(B5_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC2_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum2[[k]] <- sum(B6_sum2[[k]], weight2*coefficients(mod2)[formulaMatrix$PC3_2[f]], na.rm=T)}
        }}
        if(sum(weight_vec1[[k]])==0 |
           sum(weight_vec2[[k]]) ==0){
          pred_esm1[[k]] <- NA
          pred_esm2[[k]] <- NA

          B0_esm1[[k]] <- NA
          B1_esm1[[k]] <- NA
          B2_esm1[[k]] <- NA
          B3_esm1[[k]] <- NA
          B4_esm1[[k]] <- NA
          B5_esm1[[k]] <- NA
          B6_esm1[[k]] <- NA

          B0_esm2[[k]] <- NA
          B1_esm2[[k]] <- NA
          B2_esm2[[k]] <- NA
          B3_esm2[[k]] <- NA
          B4_esm2[[k]] <- NA
          B5_esm2[[k]] <- NA
          B6_esm2[[k]] <- NA

          aucWeighted_ESM1[[k]] <- NA
          aucWeighted_ESM2[[k]] <- NA
          RMSEWeighted_ESM1[[k]] <- NA
          RMSEWeighted_ESM2[[k]] <- NA
          PR_ESM1[[k]] <- NA
          PR_ESM2[[k]] <- NA
          TjurR2_ESM1[[k]] <- NA
          TjurR2_ESM2[[k]] <- NA
          model_curve_PC1_ESM1[[k]] <- NA
          model_curve_PC1_ESM2[[k]] <- NA
          model_curve_PC2_ESM1[[k]] <- NA
          model_curve_PC2_ESM2[[k]] <- NA
          model_curve_PC3_ESM1[[k]] <- NA
          model_curve_PC3_ESM2[[k]] <- NA

          PC1_rankCor_ESM1[[k]] <- NA
          PC1_rankCor_ESM2[[k]] <- NA

          PC2_sd_ESM1[[k]] <- NA
          PC2_sd_ESM2[[k]] <- NA

          PC3_rankCor_ESM1[[k]] <- NA
          PC3_rankCor_ESM2[[k]] <- NA
        }
        else{
          # calculate ESM prediction
          pred_esm1[[k]] <- pred_sum1[[k]]/sum(weight_vec1[[k]])
          pred_esm2[[k]] <- pred_sum2[[k]]/sum(weight_vec2[[k]])

          #Calculate betas
          B0_esm1[[k]] <- B0_sum1[[k]]/sum(weight_vec1[[k]])
          B1_esm1[[k]] <- B1_sum1[[k]]/sum(weight_vec1[[k]])
          B2_esm1[[k]] <- B2_sum1[[k]]/sum(weight_vec1[[k]])
          B3_esm1[[k]] <- B3_sum1[[k]]/sum(weight_vec1[[k]])
          B4_esm1[[k]] <- B4_sum1[[k]]/sum(weight_vec1[[k]])
          B5_esm1[[k]] <- B5_sum1[[k]]/sum(weight_vec1[[k]])
          B6_esm1[[k]] <- B6_sum1[[k]]/sum(weight_vec1[[k]])

          B0_esm2[[k]] <- B0_sum2[[k]]/sum(weight_vec2[[k]])
          B1_esm2[[k]] <- B1_sum2[[k]]/sum(weight_vec2[[k]])
          B2_esm2[[k]] <- B2_sum2[[k]]/sum(weight_vec2[[k]])
          B3_esm2[[k]] <- B3_sum2[[k]]/sum(weight_vec2[[k]])
          B4_esm2[[k]] <- B4_sum2[[k]]/sum(weight_vec2[[k]])
          B5_esm2[[k]] <- B5_sum2[[k]]/sum(weight_vec2[[k]])
          B6_esm2[[k]] <- B6_sum2[[k]]/sum(weight_vec2[[k]])

          #evaluate ESM
          aucWeighted_ESM1[[k]] <- aucWeighted(pres=pred_esm1[[k]][test1$Y==1],
                                               contrast=pred_esm1[[k]][test1$Y==0],
                                               presWeight = presWeights,
                                               contrastWeight = contrastWeights)
          aucWeighted_ESM2[[k]] <- aucWeighted(pres=pred_esm2[[k]][test2$Y==1],
                                               contrast=pred_esm2[[k]][test2$Y==0],
                                               presWeight = presWeights,
                                               contrastWeight = contrastWeights)

          RMSEWeighted_ESM1[[k]] <- computeRMSEWeighted(Y=test1$Y, predY=pred_esm1[[k]])
          RMSEWeighted_ESM2[[k]] <- computeRMSEWeighted(Y=test2$Y, predY=pred_esm2[[k]])

          PR_ESM1[[k]] <- pr.curve(scores.class0=pred_esm1[[k]],
                                   weights.class0 = test1$Y,
                                   curve=TRUE)$auc.integral
          PR_ESM2[[k]] <- pr.curve(scores.class0=pred_esm2[[k]],
                                   weights.class0 = test2$Y,
                                   curve=TRUE)$auc.integral

          TjurR2_ESM1[[k]] <- computeTjurR2(Y=test1$Y,
                                            predY=pred_esm1[[k]])
          TjurR2_ESM2[[k]] <- computeTjurR2(Y=test2,
                                            predY=pred_esm2[[k]])
          #Species Response Curves
          model_curve_PC1_ESM1[[k]] <- pnorm(B1_esm1[[k]]*x1)
          model_curve_PC1_ESM2[[k]] <- pnorm(B1_esm2[[k]]*x1)

          model_curve_PC2_ESM1[[k]] <- pnorm(B2_esm1[[k]]*x2)
          model_curve_PC2_ESM2[[k]] <- pnorm(B2_esm2[[k]]*x2)

          model_curve_PC3_ESM1[[k]] <- pnorm(B3_esm1[[k]]*x3)
          model_curve_PC3_ESM2[[k]] <- pnorm(B3_esm2[[k]]*x3)

          if(sd(model_curve_PC1_ESM1[[k]])>0) {PC1_rankCor_ESM1[[k]] <- compareResponse(model_curve_PC1_ESM1[[k]],
                                                                                        real_curve_PC1,
                                                                                        data =data.frame(x1),
                                                                                        graph=F)$rankCor}
          if(sd(model_curve_PC1_ESM2[[k]])>0) {PC1_rankCor_ESM2[[k]] <- compareResponse(model_curve_PC1_ESM2[[k]],
                                                                                        real_curve_PC1,
                                                                                        data =data.frame(x1),
                                                                                        graph=F)$rankCor}

          PC2_sd_ESM1[[k]] <- sd(model_curve_PC2_ESM1[[k]])
          PC2_sd_ESM2[[k]] <- sd(model_curve_PC2_ESM2[[k]])

          if(sd(model_curve_PC3_ESM1[[k]])>0) {PC3_rankCor_ESM1[[k]] <- compareResponse(model_curve_PC3_ESM1[[k]],
                                                                                        real_curve_PC3,
                                                                                        data =data.frame(x3),
                                                                                        graph=F)$rankCor}
          if(sd(model_curve_PC3_ESM2[[k]])>0) {PC3_rankCor_ESM2[[k]] <- compareResponse(model_curve_PC3_ESM2[[k]],
                                                                                        real_curve_PC3,
                                                                                        data =data.frame(x3),
                                                                                        graph=F)$rankCor}
        }
      }

      aucWeighted_ESM1avg <- mean(unlist(lapply(aucWeighted_ESM1,sapply,mean)))
      aucWeighted_ESM2avg <- mean(unlist(lapply(aucWeighted_ESM2, sapply, mean)))
      if(is.na(aucWeighted_ESM1avg) | is.na(aucWeighted_ESM2avg)) {model$aucWeighted <- NA}
      else{model$aucWeighted <- mean(aucWeighted_ESM1avg, aucWeighted_ESM2avg)}
      RMSEWeighted_ESM1avg <- mean(unlist(lapply(RMSEWeighted_ESM1,sapply,mean)))
      RMSEWeighted_ESM2avg <-   mean(unlist(lapply(RMSEWeighted_ESM2,sapply,mean)))
      if(is.na(RMSEWeighted_ESM1avg) | is.na(RMSEWeighted_ESM2avg)) {model$RMSEWeighted <- NA}
      else{model$RMSEWeighted <- mean(RMSEWeighted_ESM1avg,
                                      RMSEWeighted_ESM2avg)}
      PR_ESM1avg <- mean(unlist(lapply(PR_ESM1, sapply, mean)))
      PR_ESM2avg <- mean(unlist(lapply(PR_ESM2, sapply, mean)))
      if(is.na(PR_ESM1avg) | is.na(PR_ESM2avg)) {model$PR <- NA}
      else{model$PR <- mean(PR_ESM1avg,
                            PR_ESM2avg)}
      TjurR2_ESM1avg <- mean(unlist(lapply(TjurR2_ESM1, sapply, mean)))
      TjurR2_ESM2avg <- mean(unlist(lapply(TjurR2_ESM2, sapply, mean)))
      if(is.na(TjurR2_ESM1avg) | is.na(TjurR2_ESM2avg)) {model$TjurR2 <- NA}
      else{model$TjurR2 <- mean(TjurR2_ESM1avg,
                                TjurR2_ESM2avg)}

      PC1_rankCor_ESM1avg <- mean(unlist(lapply(PC1_rankCor_ESM1, sapply, mean)))
      PC1_rankCor_ESM2avg <- mean(unlist(lapply(PC1_rankCor_ESM2, sapply, mean)))
      if(is.na(PC1_rankCor_ESM1avg) | is.na(PC1_rankCor_ESM2avg)) {model$PC1_rankCor <- NA}
      else{model$PC1_rankCor <- mean(PC1_rankCor_ESM1avg, PC1_rankCor_ESM2avg)}

      PC2_sd_ESM1avg <- mean(unlist(lapply(PC2_sd_ESM1, sapply, mean)))
      PC2_sd_ESM2avg <- mean(unlist(lapply(PC2_sd_ESM2, sapply, mean)))
      if(is.na(PC2_sd_ESM1avg) | is.na(PC2_sd_ESM2avg)) {model$PC2_sd <- NA}
      else{model$PC2_sd <- mean(PC2_sd_ESM1avg,
                                PC2_sd_ESM2avg)}

      PC3_rankCor_ESM1avg <- mean(unlist(lapply(PC3_rankCor_ESM1, sapply, mean)))
      PC3_rankCor_ESM2avg <- mean(unlist(lapply(PC3_rankCor_ESM2, sapply, mean)))
      if(is.na(PC3_rankCor_ESM1avg) | is.na(PC3_rankCor_ESM2avg)) {model$PC3_rankCor <- NA}
      else{model$PC3_rankCor <- mean(PC3_rankCor_ESM1avg, PC3_rankCor_ESM2avg)}

      #Linear term for PC1
      b1_1_ESM1_avg <- mean(unlist(lapply(B1_esm1, sapply, mean)))
      b1_1_ESM2_avg <- mean(unlist(lapply(B1_esm2, sapply, mean)))
      if(is.na(b1_1_ESM1_avg) | is.na(b1_1_ESM2_avg)) {model$b1_1 <- NA}
      else{model$b1_1 <- mean(b1_1_ESM1_avg, b1_1_ESM2_avg)}

      #Quadratic term for PC1
      b1_2_ESM1_avg <- mean(unlist(lapply(B4_esm1, sapply, mean)))
      b1_2_ESM2_avg <- mean(unlist(lapply(B4_esm2, sapply, mean)))
      if(is.na(b1_2_ESM1_avg) | is.na(b1_2_ESM2_avg)) {model$b1_2 <- NA}
      else{model$b1_2 <- mean(b1_2_ESM1_avg, b1_2_ESM2_avg)}

      #Linear term for PC2
      b2_1_ESM1_avg <- mean(unlist(lapply(B2_esm1, sapply, mean)))
      b2_1_ESM2_avg <- mean(unlist(lapply(B2_esm2, sapply, mean)))
      if(is.na(b2_1_ESM1_avg) | is.na(b2_1_ESM2_avg)) {model$b2_1 <- NA}
      else{model$b2_1 <- mean(b2_1_ESM1_avg, b2_1_ESM2_avg)}

      #Quadratic term for PC2
      b2_2_ESM1_avg <- mean(unlist(lapply(B5_esm1, sapply, mean)))
      b2_2_ESM2_avg <- mean(unlist(lapply(B5_esm2, sapply, mean)))
      if(is.na(b2_2_ESM1_avg) | is.na(b2_2_ESM2_avg)) {model$b2_2 <- NA}
      else{model$b2_2 <- mean(b2_2_ESM1_avg, b2_2_ESM2_avg)}

      #Linear term for PC3
      b3_1_ESM1_avg <- mean(unlist(lapply(B3_esm1, sapply, mean)))
      b3_1_ESM2_avg <- mean(unlist(lapply(B3_esm2, sapply, mean)))
      if(is.na(b3_1_ESM1_avg) | is.na(b3_1_ESM2_avg)) {model$b3_1 <- NA}
      else{model$b3_1 <- mean(b3_1_ESM1_avg, b3_1_ESM2_avg)}

      #Quadratic term for PC3
      b3_2_ESM1_avg <- mean(unlist(lapply(B6_esm1, sapply, mean)))
      b3_2_ESM2_avg <- mean(unlist(lapply(B6_esm2, sapply, mean)))
      if(is.na(b3_2_ESM1_avg) | is.na(b3_2_ESM2_avg)) {model$b3_2 <- NA}
      else{model$b3_2 <- mean(b3_2_ESM1_avg, b3_2_ESM2_avg)}


      save(model, file=paste0("../data/models/",
                              modelType[2], "/",species[s], "/",
                              sizes[n], "/", "model_",replicates[r],
                              ".RData"))
      cat(species[s], sizes[n], replicates[r], "\n")

    }
  }
}
runTime <- Sys.time() - timeStart #1.63 mins

#Process ESM models
results <- expand.grid(species, sizes, replicates)
results <- results[order(results$Var1, results$Var2, results$Var3),]
names(results) <- c("species", "size", "rep")
results$AUC_weighted <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat <- rep(NA, length(results$species))
results$samples <- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$pr_davis_goadrich <- rep(NA, length(results$species))
results$model <- rep("ESM", length(results$species))
results$PC2_sd <- rep(NA, length(results$species))
results$b1_1 <- rep(NA, length(results$species))
results$b1_2 <- rep(NA, length(results$species))
results$b2_1 <- rep(NA, length(results$species))
results$b2_2 <- rep(NA, length(results$species))
results$b3_1 <- rep(NA, length(results$species))
results$b3_2 <- rep(NA,length(results$species))

timeStart <- Sys.time()
for ( i in 1:length(results$species)){
  if(file.exists(paste0("../data/models/",
                        modelType[1], "/",results$species[i], "/",
                        results$size[i], "/", "model_",results$rep[i],
                        ".RData"))){
    load(paste0("../data/models/",
                modelType[1], "/",results$species[i], "/",
                results$size[i], "/", "model_",results$rep[i],
                ".RData"))

    results$PC1_rankCor[i] <- model$PC1_rankCor
    results$PC2_sd[i] <- model$PC2_sd
    results$PC3_rankCor[i] <- model$PC3_rankCor
    results$AUC_weighted[i] <- model$aucWeighted
    results$RMSEWeighted[i] <- model$RMSEWeighted
    results$pr_integral[i] <- model$PR
    results$TjursR2[i] <- model$TjurR2
    results$b1_1[i] <- model$b1_1
    results$b1_2[i] <- model$b1_2
    results$b2_1[i] <- model$b2_1
    results$b2_2[i] <- model$b2_2
    results$b3_1[i] <- model$b3_1
    results$b3_2[i] <- model$b3_2
}


  if(i %% 100== 0) { cat(paste0(i, " of ", length(results$species), "\n"))}
}

timeEnd <- Sys.time() - timeStart #1.17s

save(results, file="../data/models/ESM/results2.RData")


#Process ESM_bivariate results
results <- expand.grid(species, sizes, replicates)
results <- results[order(results$Var1, results$Var2, results$Var3),]
names(results) <- c("species", "size", "rep")
results$AUC_weighted <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat <- rep(NA, length(results$species))
results$samples <- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$pr_davis_goadrich <- rep(NA, length(results$species))
results$model <- rep("ESM", length(results$species))
results$PC2_sd <- rep(NA, length(results$species))
results$b1_1 <- rep(NA, length(results$species))
results$b1_2 <- rep(NA, length(results$species))
results$b2_1 <- rep(NA, length(results$species))
results$b2_2 <- rep(NA, length(results$species))
results$b3_1 <- rep(NA, length(results$species))
results$b3_2 <- rep(NA,length(results$species))

timeStart <- Sys.time()
for ( i in 1:length(results$species)){
  if(file.exists(paste0("../data/models/",
                        modelType[2], "/",results$species[i], "/",
                        results$size[i], "/", "model_",results$rep[i],
                        ".RData"))){
    load(paste0("../data/models/",
                modelType[2], "/",results$species[i], "/",
                results$size[i], "/", "model_",results$rep[i],
                ".RData"))

    results$PC1_rankCor[i] <- model$PC1_rankCor
    results$PC2_sd[i] <- model$PC2_sd
    results$PC3_rankCor[i] <- model$PC3_rankCor
    results$AUC_weighted[i] <- model$aucWeighted
    results$RMSEWeighted[i] <- model$RMSEWeighted
    results$pr_integral[i] <- model$PR
    results$TjursR2[i] <- model$TjurR2
    results$b1_1[i] <- model$b1_1
    results$b1_2[i] <- model$b1_2
    results$b2_1[i] <- model$b2_1
    results$b2_2[i] <- model$b2_2
    results$b3_1[i] <- model$b3_1
    results$b3_2[i] <- model$b3_2
  }


  if(i %% 100== 0) { cat(paste0(i, " of ", length(results$species), "\n"))}
}

timeEnd <- Sys.time() - timeStart #1.14s

save(results, file="../data/models/ESM_bivariate/results2.RData")


###Explore potential overfitting in full ESM
save(weights_df, file="../data/models/ESM/weights_df.RData")

#Sample Size 2 Overfitting

#Ensure weights add up to 1 for each row
weights_df[,4:29] <- weights_df[,4:29]/sum(weights_df[,4:29])

weights_df$cumOverfit <- rep(0, length(weights_df$species))

for (i in 1:length(weights_df$species)) {
  if(weights_df$sizes[i] == "size2") {
    weights_df$cumOverfit[i] <- sum(weights_df[i,3+ 4:26])
  }
  if(weights_df$sizes[i] == "size4") {
    weights_df$cumOverfit[i] <- sum(weights_df[i, 3+17:26])
  }

}

weights_df_overfit <- subset(weights_df, weights_df$cumOverfit >0)

ggplot(weights_df_overfit, aes(x=cumOverfit,
      color=sizes, fill=sizes)) +
  geom_histogram(alpha=0.5,
                 position="identity") +
  xlab("Weight of Overfit Models")