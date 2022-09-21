### Evaluate models
library(MCMCvis)
library(enmSdm)
library(PRROC)
memory.limit(2^32) 
# A script for pulling out model evaluation criteria
# from the models that were run in `HMSC_single_simple` and
# `HMSC_joint_simple`

load("../data/south.RData")
#Load simulated data
load(paste0("../data/sims_data.RData"))
NSites <- length(south$long)
row.names(south) <- 1:NSites
replicates <- paste0("rep", 1:30)
NPresences <- c(64, 32, 16, 8, 4, 2)
sizes <- paste0("size", NPresences)# Desired sample size
species <- c("wide_avg", "wide_ext",
             "narrow_avg", "narrow_ext")
splits <- paste0("split", 1:3)
modelType  <- "HMSC_joint_simple"


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
results$samples <- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$pr_davis_goadrich <- rep(NA, length(results$species))
results$PC2_sd <- rep(NA, length(results$species))


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
for ( i in thing){
  if(file.exists(paste0("../data/models/",
                        modelType[1], "/",results$species[i], "/",
                        results$size[i], "/", "model_",results$rep[i],
                        ".RData"))){
    load(paste0("../data/models/",
                modelType[1], "/",results$species[i], "/",
                results$size[i], "/", "model_",results$rep[i],
                ".RData"))
    results$samples[i] <- model$samples
    results$Rhat[i] <- model$maxRHat
    results$PC1_rankCor[i] <- model$PC1_rankCor
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
    
    results$TjursR2[i] <- mean(model$split1_metrics$TjurR2[64],
                               model$split2_metrics$TjurR2[64],
                               model$split3_metrics$TjurR2[64])
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
    
    if(names(model)[22] == "m.post1_1"){
    #Split 1
    betas1_1 <- MCMCsummary(model$m.post1_1$Beta)
    betas1_2 <- MCMCsummary(model$m.post1_2$Beta)

    
    
    #Split2
    betas2_1 <- MCMCsummary(model$m.post2_1$Beta)
    betas2_2 <- MCMCsummary(model$m.post2_2$Beta)

    #Split3
    betas3_1 <- MCMCsummary(model$m.post3_1$Beta)
    betas3_2 <- MCMCsummary(model$m.post3_2$Beta)
    }
    
    if(names(model)[28]== "beta_chain1_1") {
      betas1_1 <- MCMCsummary(model$beta_chain1_1)
      betas1_2 <- MCMCsummary(model$beta_chain1_2)
      
      #Split2
      betas2_1 <- MCMCsummary(model$beta_chain2_1)
      betas2_2 <- MCMCsummary(model$beta_chain2_2)
      
      #Split3
      betas3_1 <- MCMCsummary(model$beta_chain3_1)
      betas3_2 <- MCMCsummary(model$beta_chain3_2)
    }

    model_curve_PC2_m1_1 <- pnorm(betas1_1$mean[444]*x2 + betas1_1$mean[447]*x2*x2)
    model_curve_PC2_m1_2 <- pnorm(betas1_2$mean[444]*x2 + betas1_2$mean[447]*x2*x2)
    
    model_curve_PC2_m2_1 <- pnorm(betas2_1$mean[444]*x2 + betas2_1$mean[447]*x2*x2)
    model_curve_PC2_m2_2 <- pnorm(betas2_2$mean[444]*x2 + betas2_2$mean[447]*x2*x2)
    
    model_curve_PC2_m3_1 <- pnorm(betas3_1$mean[444]*x2 + betas3_1$mean[447]*x2*x2)
    model_curve_PC2_m3_2 <- pnorm(betas3_2$mean[444]*x2 + betas3_2$mean[447]*x2*x2)
    
    results$PC2_sd[i] <- mean(c(sd(model_curve_PC2_m1_1),
                           sd(model_curve_PC2_m1_2),
                           sd(model_curve_PC2_m2_1),
                           sd(model_curve_PC2_m2_2),
                           sd(model_curve_PC2_m3_1),
                           sd(model_curve_PC2_m3_2)))
    
    
  }else{print("already done")}
      

if(i %% 2== 0) { cat(paste0(i, " of ", length(results$species), "\n"))}
}
timeEnd <- Sys.time() - timeStart

save(results, file="../data/models/HMSC_joint_simple/results.RData")

load("../data/models/HMSC_joint_simple/results.RData")
for ( i in 343:length(results$species)){
  if(file.exists(paste0("../data/models/",
                        modelType[1], "/",results$species[i], "/",
                        results$size[i], "/", "model_",results$rep[i],
                        ".RData"))){
    load(paste0("../data/models/",
                modelType[1], "/",results$species[i], "/",
                results$size[i], "/", "model_",results$rep[i],
                ".RData"))
    if(names(model)[22] == "m.post1_1"){
      #Split 1
      betas1_1 <- MCMCsummary(model$m.post1_1$Beta)
      betas1_2 <- MCMCsummary(model$m.post1_2$Beta)
      
      
      
      #Split2
      betas2_1 <- MCMCsummary(model$m.post2_1$Beta)
      betas2_2 <- MCMCsummary(model$m.post2_2$Beta)
      
      #Split3
      betas3_1 <- MCMCsummary(model$m.post3_1$Beta)
      betas3_2 <- MCMCsummary(model$m.post3_2$Beta)
    }
    
    if(names(model)[28]== "beta_chain1_1") {
      betas1_1 <- MCMCsummary(model$beta_chain1_1)
      betas1_2 <- MCMCsummary(model$beta_chain1_2)
      
      #Split2
      betas2_1 <- MCMCsummary(model$beta_chain2_1)
      betas2_2 <- MCMCsummary(model$beta_chain2_2)
      
      #Split3
      betas3_1 <- MCMCsummary(model$beta_chain3_1)
      betas3_2 <- MCMCsummary(model$beta_chain3_2)
    }
    
    #linear term for PC1
    results$b1_1[i] <- mean(betas1_1$mean[443], betas1_2$mean[443],
                             betas2_1$mean[443], betas2_2$mean[443],
                             betas3_1$mean[443], betas3_2$mean[443], na.rm=T)
    #quadratic term for PC1
    results$b1_2[i] <- mean(betas1_1$mean[446], betas1_2$mean[446],
                            betas2_1$mean[446], betas2_2$mean[446],
                            betas3_1$mean[446], betas3_2$mean[446], na.rm=T)
    #linear term for PC2
    results$b2_1[i] <- mean(betas1_1$mean[444], betas1_2$mean[444],
                            betas2_1$mean[444], betas2_2$mean[444],
                            betas3_1$mean[444], betas3_2$mean[444], na.rm=T)
    #quadratic term for PC2
    results$b2_2[i] <- mean(betas1_1$mean[447], betas1_2$mean[447],
                            betas2_1$mean[447], betas2_2$mean[447],
                            betas3_1$mean[447], betas3_2$mean[447], na.rm=T)
    #linear term for PC3
    results$b3_1[i] <- mean(betas1_1$mean[445], betas1_2$mean[445],
                            betas2_1$mean[445], betas2_2$mean[445],
                            betas3_1$mean[445], betas3_2$mean[445], na.rm=T)
    
    results$b3_2[i] <- mean(betas1_1$mean[448], betas1_2$mean[448],
                            betas2_1$mean[448], betas2_2$mean[448],
                            betas3_1$mean[448], betas3_2$mean[448], na.rm=T)
    if(i %% 2== 0) { cat(paste0(i, " of ", length(results$species), "\n"))}
  }}

save(results, file="../data/models/HMSC_joint_simple/results2.RData")
