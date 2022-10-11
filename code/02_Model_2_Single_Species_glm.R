#02_Model_2_Single_Species_glm

#In this script we use a basic glm



if(Sys.info()['sysname'] == "Darwin") {
  path <- "/Users/curculion/Documents/GitHub"
  path2 <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
  path <- "C:/Users/kerickson/Documents/GitHub"
  path2 <- "H:/Global Change Program/Research/ENMs - Modeling Methods for Rare Species (Kelley Erickson)/rare_species/data"
}


source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

modelType  <- models[2]

timeStart <- Sys.time()

for( r in 1:30){
  for (s in 1:length(species)){
    for (n in 1:length(sizes)){

      timeStart1 <- Sys.time()
      model <-list()
      model$NPresesences <- sizes[n]
      model$species <- species[s]
      model$replicate <- replicates[r]

      model$modelType <- modelType[1]
      Y <- as.matrix(sims[[s]][[n]][[r]]$YSim)
      XData <- south[,3:5]

      Y_train <- Y[sims[[s]][[n]][[r]]$train]
      Y_test <- Y[sims[[s]][[n]][[r]]$test]
      X_train <- XData[sims[[s]][[n]][[r]]$train,]
      X_test <- XData[sims[[s]][[n]][[r]]$test,]
      train <- cbind(Y_train, X_train)
      if (n==1) {
        m <- glm(Y_train ~ PC1 + PC2 + PC3,
                 data = train, family=binomial())

      }
      if(n >1) {
        m <- glm(Y_train ~ PC1 + PC2 + PC3 +
                   I(PC1^2) + I(PC2^2) + I(PC3^2),
                 data=train, family=binomial())
      }


      model$preds <- predict(m,newdata=X_test, type="response")




      presWeights <- rep(1/64, 64)
      contrastWeights <- rep(1/64, 64)


      model$auc <- aucWeighted(pres=model$preds[1:64],
                               contrast=model$preds[65:128],
                               presWeight = presWeights,
                               contrastWeight = contrastWeights)

       model$RMSEWeighted <- computeRMSEWeighted(Y_test, model$preds)
      model$TjursR2 <- computeR2(Y_test, model$preds)
      #preds2 <- computePredictedValues(m) #fits to training data
      #preds.mean <- apply(model$preds, FUN='mean', margin='1')



      if(n == 1) {
        model_curve_PC1 <- pnorm(m$coefficients[2]*x1)
        model_curve_PC2 <- pnorm(m$coefficients[3]*x2)
        model_curve_PC3 <- pnorm(m$coefficients[4]*x3)
      }

      if(n > 1) {
        model_curve_PC1 <- pnorm(m$coefficients[2]*x1 + m$coefficients[5]*x1*x1)
        model_curve_PC2 <- pnorm(m$coefficients[3]*x2 + m$coefficients[6]*x2*x2)
        model_curve_PC3 <- pnorm(m$coefficients[4]*x3 + m$coefficients[7]*x3*x3)
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
          response$L[j] <- m$coefficients[1] + m$coefficients[2]*response$PC1[j] +
            m$coefficients[3]*response$PC2[j] + m$coefficients[4]*response$PC3[j]      }

      }
      if(n > 1) {
        for ( j in 1:length(response$PC1)) {
          response$L[j] <- m$coefficients[1] + m$coefficients[2]*response$PC1[j] +
           m$coefficients[3]*response$PC2[j] + m$coefficients[4]*response$PC3[j] +
            m$coefficients[5]*response$PC1[j]*response$PC1[j] +
            m$coefficients[6]*response$PC2[j]*response$PC2[j] +
            m$coefficients[7]*response$PC3[j]*response$PC3[j]
        }
      }

      model$varPC2Response<- sd(invlogit(response$L))

      #Precision Recall
      PR <- pr.curve(scores.class0=model$preds,
                     weights.class0 = Y_test,
                     curve=TRUE)
      model$pr_integral <- PR$auc.integral
      model$timeEnd <- Sys.time() - timeStart1
      save(model, file=paste0(path2, "/models/",
                              modelType[1], "/",species[s], "/",
                              sizes[n], "/", "model_",replicates[r],
                              ".RData"))

    }
  }
}
timeEnd <- Sys.time()-timeStart

