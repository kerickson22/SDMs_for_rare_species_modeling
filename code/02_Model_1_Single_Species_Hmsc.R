#02_Model_1_Single_Species_Hmsc

#In this script we use the single-species variant
# of Hmsc (essentially a Bayesian glm).

path <- "/Users/curculion/Documents/GitHub"
source(paste0(path, "/SDMs_for_rare_species_modeling/code/00_Constants.R"))

modelType  <- models[1]
nChains <- 2
thin <- 1
samples <- 15000
transient <-2000
verbose <-1000


timeStart <- Sys.time()
for (s in 1:length(species)){
  for (n in 1:length(sizes)){
    for( r in 1:length(replicates)){
      timeStart1 <- Sys.time()
      model <-list()
      model$NPresesences <- sizes[n]
      model$species <- species[s]
      model$replicate <- replicates[r]
      model$nChains <- nChains
      model$samples <- samples
      model$transient <- transient
      model$thin <-thin
      model$modelType <- modelType[1]
      Y <- as.matrix(sims[[s]][[n]][[r]]$YSim)
      XData <- south[,3:5]

      Y_train <- Y[sims[[s]][[n]][[r]]$train]
      Y_test <- Y[sims[[s]][[n]][[r]]$test]
      X_train <- XData[sims[[s]][[n]][[r]]$train,]
      X_test <- XData[sims[[s]][[n]][[r]]$test,]
      m <- Hmsc(Y=Y_train, XData=X_train,
                      XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                        I(PC2^2)+I(PC3^2),
                      distr="probit")
      m <- sampleMcmc(m, thin=thin, samples=samples,
                 transient=transient,nParallel=2,
                 nChains=nChains, verbose=verbose)
      model$preds <- predict(m,XData=X_test, expected=F) #expected=T, probs, expected =F, 0,1
      Matrix_x <- matrix(unlist(model$preds), ncol = nChains*samples, byrow = TRUE)
      preds.mean <- apply(Matrix_x, FUN="mean", MARGIN=1)
     # preds.mean <- data.frame(preds.mean)
     # rownames(preds.mean) <- rownames(X_test)
      presWeights <- rep(1/64, 64)
      contrastWeights <- rep(1/64, 64)


      model$auc <- aucWeighted(pres=preds.mean[1:64],
                                           contrast=preds.mean[65:128],
                                           presWeight = presWeights,
                                           contrastWeight = contrastWeights)
      model$RMSEWeighted <- computeRMSEWeighted(Y_test, preds.mean)
      model$TjursR2 <- computeR2(Y_test, preds.mean)
      #preds2 <- computePredictedValues(m) #fits to training data
      #preds.mean <- apply(model$preds, FUN='mean', margin='1')
      mpost <- convertToCodaObject(m)
      betas <- MCMCsummary(mpost$Beta)
      model$maxRhat <- max(betas$Rhat)

      model_curve_PC1 <- pnorm(betas$mean[2]*x1 + betas$mean[5]*x1*x1)
      model_curve_PC2 <- pnorm(betas$mean[3]*x2 + betas$mean[6]*x2*x2)
      model_curve_PC3 <- pnorm(betas$mean[4]*x3 + betas$mean[7]*x3*x3)

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


        for ( j in 1:length(response$PC1)) {
          response$L[j] <- betas$mean[1] + betas$mean[2]*response$PC1[j] +
            betas$mean[3]*response$PC2[j] + betas$mean[4]*response$PC3[j] +
            betas$mean[5]*response$PC1[j]*response$PC1[j] +
            betas$mean[6]*response$PC2[j]*response$PC2[j] +
            betas$mean[7]*response$PC3[j]*response$PC3[j]
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
