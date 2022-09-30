#02_Model_6_SAM

#In this script we run SAM model, using the
# number of aggregates identified in script
# `02_Model_6a_Real_Species_SAM`

path <- "/Users/curculion/Documents/GitHub"
source(paste0(path, "/SDMs_for_rare_species_modeling/code/00_Constants.R"))

modelType  <- models[6]



timeStart <- Sys.time()
for (s in 1:length(species)){
  for (n in 1:length(sizes)){
    for( r in 1:length(replicates)){
      timeStart1 <- Sys.time()
      model <-list()
      model$NPresesences <- sizes[n]
      model$species <- species[s]
      model$replicate <- replicates[r]

      model$modelType <- modelType[1]
      Y <- cbind(south[,6:67],as.matrix(sims[[s]][[n]][[r]]$YSim))
      XData <- south[,3:5]

      Y_train <- Y[sims[[s]][[n]][[r]]$train,]
      Y_test <- Y[sims[[s]][[n]][[r]]$test,]
      X_train <- XData[sims[[s]][[n]][[r]]$train,]
      X_test <- XData[sims[[s]][[n]][[r]]$test,]
      archetype_formula <- as.formula(paste0(paste0('cbind(',paste(colnames(Y),collapse = ", "),") ~ PC1 + PC2 + PC3 + I(PC1^2) + I(PC2^2) + I(PC3^2)")))

      ## Species formula
      species_formula <- ~ 1

      m <- species_mix(archetype_formula = archetype_formula, # Archetype formula
                              species_formula = species_formula,    # Species formula
                              data = cbind(X_train,Y_train),            # Data
                              nArchetypes = 6,              # Number of groups (mixtures) to fit
                              family = 'bernoulli',         # Which family to use
                              control = list(quiet = TRUE))



      model$preds <- predict(m, newdata=X_test,
                                     type="response", prediction.type="species")
      preds.mean <- model$preds[,"SimSp"]


      presWeights <- rep(1/64, 64)
      contrastWeights <- rep(1/64, 64)


      model$auc <- aucWeighted(pres=preds.mean[1:64],
                               contrast=preds.mean[65:128],
                               presWeight = presWeights,
                               contrastWeight = contrastWeights)
      model$RMSEWeighted <- computeRMSEWeighted(Y_test["SimSp"], preds.mean)
      model$TjursR2 <- computeR2(Y_test["SimSp"
      ], preds.mean)
      #preds2 <- computePredictedValues(m) #fits to training data
      #preds.mean <- apply(model$preds, FUN='mean', margin='1')

      betas <- m$beta[grep(names(which(round(m$tau["SimSp",])==1)), names(m$beta))]
      names(betas) <- c("PC1", "PC2", "PC3",
                        "I(PC1^2)",
                        "I(PC2^2)",
                        "I(PC3^2)")

      model_curve_PC1 <- pnorm(betas[1]*x1 + betas[4]*x1*x1)
      model_curve_PC2 <- pnorm(betas[2]*x2 + betas[5]*x2*x2)
      model_curve_PC3 <- pnorm(betas[3]*x3 + betas[6]*x3*x3)

      real_curve_PC1 <- dnorm(x1, mean=sims[[s]]$mu[1],
                              sd=sims[[s]]$sd[1,1])
      real_curve_PC3 <- dnorm(x3, mean=sims[[s]]$mu[2],
                              sd=sims[[s]]$sd[2,2])
      model$PC1_cor<- compareResponse(model_curve_PC1,
                                      real_curve_PC1,
                                      data =data.frame(x1),
                                      graph=F)$cor
      model$PC3_cor<- compareResponse(model_curve_PC3,
                                      real_curve_PC3,
                                      data =data.frame(x3),
                                      graph=F)$cor
      model$PC1_rankCor <- compareResponse(model_curve_PC1,
                                           real_curve_PC1,
                                           data =data.frame(x1),
                                           graph=F)$rankCor
      model$PC3_rankCor <- compareResponse(model_curve_PC3,
                                           real_curve_PC3,
                                           data =data.frame(x3),
                                           graph=F)$rankCor
      #Calculate PC2 response
      x1_PC2_response <- mu_avg[1]
      x2_PC2_response <- south$PC2
      x3_PC2_response <-mu_avg[2]

      response <- expand.grid(x1_PC2_response,
                              x2_PC2_response,
                              x3_PC2_response)
      names(response) <- c("PC1", "PC2", "PC3")

      #grep("SimSp", rownames(betas))
      for ( j in 1:length(response$PC1)) {
        response$L[j] <- m$alpha["SimSp"] + betas[1]*response$PC1[j] +
          betas[2]*response$PC2[j] + betas[3]*response$PC3[j] +
          betas[4]*response$PC1[j]*response$PC1[j] +
          betas[5]*response$PC2[j]*response$PC2[j] +
          betas[6]*response$PC3[j]*response$PC3[j]
      }
      model$varPC2Response<- sd(invlogit(response$L))
      model$timeEnd <- Sys.time() - timeStart1
      save(model, file=paste0("../data/models/",
                              modelType[1], "/",species[s], "/",
                              sizes[n], "/", "model_",replicates[r],
                              ".RData"))
    }
  }
}
timeEnd <- Sys.time()-timeStart


