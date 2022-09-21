### Evaluate models
library(MCMCvis)
library(enmSdm)
library(PRROC)

# A script for pulling out model evaluation criteria
# from the models that were run in `HMSC_single_simple` and
#

load("../data/south.RData")
load(paste0("../data/sims_data.RData"))
NSites <- length(south$long)
row.names(south) <- 1:NSites
replicates <- paste0("rep", 1:30)
NPresences <- c(64, 32, 16, 8, 4, 2)
sizes <- paste0("size", NPresences)# Desired sample size
species <- c("wide_avg", "wide_ext",
             "narrow_avg", "narrow_ext")
splits <- paste0("split", 1:3)
modelType  <- "HMSC_single_simple"


sites <- 1:NSites
x1 <- seq(from=-9, to = 13, length.out=100)
x2 <- seq(from=-9, to = 16, length.out=100)
x3 <- seq(from=-11, to = 4, length.out=100)


results <- expand.grid(species, sizes, replicates)
results <- results[order(results$Var1, results$Var2, results$Var3),]
names(results) <- c("species", "size", "rep")
results$AUC_weighted <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat<- rep(NA, length(results$species))
results$PC2_sd <- rep(NA, length(results$species))
results$b1_1 <- rep(NA, length(results$species))
results$b1_2 <- rep(NA, length(results$species))
results$b2_1 <- rep(NA, length(results$species))
results$b2_2 <- rep(NA, length(results$species))
results$b3_1 <- rep(NA, length(results$species))
results$b3_2 <- rep(NA, length(results$species))

load("../data/X_bar.RData")
mu_ext <- c(-3.5, 2)
Sigma_narrow <- diag(c(0.9, 2.5))

mu_avg <- c(X_bar$mu_PC1[64], X_bar$mu_PC3[64])
Sigma_wide <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.90, na.rm=T)^2,  quantile(X_bar$sd_PC3, probs=0.90, na.rm=T)^2))

computeRMSEWeighted = function(Y, predY) {
  RMSE <- sqrt((mean((Y[Y==1]-predY[Y==1])^2, na.rm=TRUE) +
                  mean((Y[Y==0]-predY[Y==0])^2, na.rm=TRUE)))
  return(RMSE)
}

timeStart <- Sys.time()
for ( i in 1:length(results$species)){

      load(paste0("../data/models/",
                 modelType[1], "/",results$species[i], "/",
                 results$size[i], "/", "model_",results$rep[i],
                 ".RData"))
      results$samples[i] <- model$samples
      results$Rhat[i] <- model$maxRHat
      results$PC1_rankCor[i] <- model$PC1_rankCor
      results$PC2_sd[i] <- model$PC2_sd
      results$PC3_rankCor[i] <- model$PC3_rankCor
      #results$AUC[i] <- mean(model$split1_metrics$AUC,
                             #model$split2_metrics$AUC,
                             #model$split3_metrics$AUC)
      results$AUC_weighted[i] <- mean(model$aucWeighted_1_1,
                                       model$aucWeighted_1_2,
                                       model$aucWeighted_2_1,
                                       model$aucWeighted_2_2,
                                       model$aucWeighted_3_1,
                                       model$aucWeighted_3_2)
      results$TjursR2[i] <- mean(model$split1_metrics$TjurR2,
                                 model$split2_metrics$TjurR2,
                                 model$split3_metrics$TjurR2)
      #results$RMSE[i] <- mean(model$split1_metrics$RMSE,
                             # model$split2_metrics$RMSE,
                              #model$split3_metrics$RMSE)

      Y <- as.matrix(sims[[results$species[i]]][[results$size[i]]][[results$rep[i]]]$YSim)
      partition1 <- sims[[results$species[i]]][[results$size[i]]][[results$rep[i]]]$split1
      partition2 <- sims[[results$species[i]]][[results$size[i]]][[results$rep[i]]]$split2
      partition3 <- sims[[results$species[i]]][[results$size[i]]][[results$rep[i]]]$split3

      Y1_1 <- Y[partition1==1]
      Y1_2 <- Y[partition1==2]
      Y2_1 <- Y[partition2==1]
      Y2_2 <- Y[partition2==2]
      Y3_1 <- Y[partition3==1]
      Y3_2 <- Y[partition3==2]

      preds1_mean1 <- model$preds1_mean[partition1==1]
      preds1_mean2 <- model$preds1_mean[partition1==2]
      preds2_mean1 <- model$preds2_mean[partition2==1]
      preds2_mean2 <- model$preds2_mean[partition2==2]
      preds3_mean1 <- model$preds3_mean[partition3==1]
      preds3_mean2 <- model$preds3_mean[partition3==2]



      RMSEWeighted_1_1 <- computeRMSEWeighted(Y1_1, preds1_mean1)
      RMSEWeighted_1_2 <- computeRMSEWeighted(Y1_2, preds1_mean2)
      RMSEWeighted_2_1 <- computeRMSEWeighted(Y2_1, preds2_mean1)
      RMSEWeighted_2_2 <- computeRMSEWeighted(Y2_2, preds2_mean2)
      RMSEWeighted_3_1 <- computeRMSEWeighted(Y3_1, preds3_mean1)
      RMSEWeighted_3_2 <- computeRMSEWeighted(Y3_2, preds3_mean2)

      results$RMSEWeighted[i] <- mean(RMSEWeighted_1_1,
                                      RMSEWeighted_1_2,
                                      RMSEWeighted_2_1,
                                      RMSEWeighted_2_2,
                                      RMSEWeighted_3_1,
                                      RMSEWeighted_3_2)
      PR_1_1 <- pr.curve(scores.class0=preds1_mean1,
                         weights.class0 = Y1_1,
                         curve=TRUE) #is this correct?
      PR_1_2 <- pr.curve(scores.class0=preds1_mean2,
                         weights.class0 = Y1_2,
                         curve=TRUE)
      PR_2_1 <- pr.curve(scores.class0=preds2_mean1,
                         weights.class0 = Y2_1,
                         curve=TRUE) #is this correct?
      PR_2_2 <- pr.curve(scores.class0=preds2_mean2,
                         weights.class0 = Y2_2,
                         curve=TRUE)
      PR_3_1 <- pr.curve(scores.class0=preds3_mean1,
                         weights.class0 = Y3_1,
                         curve=TRUE) #is this correct?
      PR_3_2 <- pr.curve(scores.class0=preds3_mean2,
                         weights.class0 = Y3_2,
                         curve=TRUE)
      results$pr_integral[i] <- mean(PR_1_1$auc.integral,
                                     PR_1_2$auc.integral,
                                     PR_2_1$auc.integral,
                                     PR_2_2$auc.integral,
                                     PR_3_1$auc.integral,
                                     PR_3_2$auc.integral)
      results$pr_davis_goadrich[i] <- mean(PR_1_1$auc.davis.goadrich,
                                           PR_1_2$auc.davis.goadrich,
                                           PR_2_1$auc.davis.goadrich,
                                           PR_2_2$auc.davis.goadrich,
                                           PR_3_1$auc.davis.goadrich,
                                           PR_3_2$auc.davis.goadrich)
    results$b1_1[i] <- model$b1_1
    results$b1_2[i] <- model$b1_2
    results$b2_1[i] <- model$b2_1
    results$b2_2[i] <- model$b2_2
    results$b3_1[i] <- model$b3_1
    results$b3_2[i] <- model$b3_2

if(i %% 20== 0) { cat(paste0(i, " of ", length(results$species), "\n"))}
}
timeEnd <- Sys.time() - timeStart

save(results,
     file="../data/models/HMSC_single_simple/results2.RData")
