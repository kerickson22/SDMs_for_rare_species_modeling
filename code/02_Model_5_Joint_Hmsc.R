#02_Model_5_Joint_Hmsc

#In this script we run Hmsc for joint species with
# latent factors

if(Sys.info()['sysname'] == "Darwin") {
path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
path <- "C:/Users/kerickson/Documents/GitHub"
path2 <- "H:/Global Change Program/Research/ENMs - Modeling Methods for Rare Species (Kelley Erickson)/rare_species/data"
}


source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

modelType  <- models[5]
nChains <- 2
thin <- 1
samples <- 20000
transient <-2000
verbose <-1000


timeStart <- Sys.time()

for( r in repStart:repEnd){
  for (s in 1:length(species)){
    for (n in 1:length(sizes)){

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
      Y <- cbind(south[,6:67],as.matrix(sims[[s]][[n]][[r]]$YSim))
      XData <- south[,3:5]

      Y_train <- Y[sims[[s]][[n]][[r]]$train,]
      Y_test <- Y[sims[[s]][[n]][[r]]$test,]
      X_train <- XData[sims[[s]][[n]][[r]]$train,]
      X_test <- XData[sims[[s]][[n]][[r]]$test,]

	if (n ==1) {
  		 m <- Hmsc(Y=Y_train, XData=X_train,
                XFormula = ~PC1 + PC2 + PC3,
                distr="probit")
	}
	if(n >1) {

      m <- Hmsc(Y=Y_train, XData=X_train,
                XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                  I(PC2^2)+I(PC3^2),
                distr="probit")
	}

      m <- sampleMcmc(m, thin=thin, samples=samples,
                      transient=transient,
                      nChains=nChains, verbose=verbose)
      model$preds <- predict(m,XData=X_test, expected=F) #expected=T, probs, expected =F, 0,1

      thing <- lapply(
        X = model$preds,
        FUN = function(x){x[,63]}
      )
       Matrix_x <- matrix(unlist(thing), ncol = nChains*samples, byrow = TRUE)
      rm(thing)

       preds.mean <- apply(Matrix_x, FUN="mean", MARGIN=1)

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
      mpost <- convertToCodaObject(m)
      betas <- MCMCsummary(mpost$Beta)
      model$maxRhat <- max(betas$Rhat)
      #grep("SimSp", rownames(betas))
	if(n == 1) {
		model_curve_PC1 <- pnorm(betas$mean[250]*x1)
     		model_curve_PC2 <- pnorm(betas$mean[251]*x2)
      	model_curve_PC3 <- pnorm(betas$mean[252]*x3)
	}

	if(n > 1) {
      model_curve_PC1 <- pnorm(betas$mean[436]*x1 + betas$mean[439]*x1*x1)
      model_curve_PC2 <- pnorm(betas$mean[437]*x2 + betas$mean[440]*x2*x2)
      model_curve_PC3 <- pnorm(betas$mean[438]*x3 + betas$mean[441]*x3*x3)
	}

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
	if(n == 1) {
		for ( j in 1:length(response$PC1)) {
        		response$L[j] <- betas$mean[249] + betas$mean[250]*response$PC1[j] +
          		betas$mean[251]*response$PC2[j] + betas$mean[252]*response$PC3[j]      }

	}
	if(n > 1) {
      	for ( j in 1:length(response$PC1)) {
        		response$L[j] <- betas$mean[435] + betas$mean[436]*response$PC1[j] +
          		betas$mean[437]*response$PC2[j] + betas$mean[438]*response$PC3[j] +
          		betas$mean[439]*response$PC1[j]*response$PC1[j] +
          		betas$mean[440]*response$PC2[j]*response$PC2[j] +
          		betas$mean[441]*response$PC3[j]*response$PC3[j]
      	}
	}

      model$varPC2Response<- sd(invlogit(response$L))

      #Precision Recall
      PR <- pr.curve(scores.class0=preds.mean,
                         weights.class0 = Y_test[["SimSp"]],
                         curve=TRUE)
      model$pr_integral <- PR$auc.integral
      model$timeEnd <- Sys.time() - timeStart1
      save(model, file=paste0(path2, "/models/",
                              modelType[1], "/",species[s], "/",
                              sizes[n], "/", "model_",replicates[r],
                              ".RData"))
     
        status <- model$timeEnd
        save(status, file=paste0(path, "/SDMs_for_rare_species_modeling/data/models/",
                              modelType[1], "/status/finished_", species[s], "_", sizes[n],"_", replicates[r], ".RData" ))
        myfile <- paste0(path, "/SDMs_for_rare_species_modeling/data/models/",
                              modelType[1], "/status/finished_", species[s], "_", sizes[n],"_", replicates[r], ".RData" ) %>%
          drive_upload(paste0("status_updates_for_Hmsc_joint/finished_",species[s], "_", sizes[n],"_", replicates[r], ".RData"  ))

      
    }
  }
}

