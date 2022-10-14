#02_Model_4_ESM_Complex

#In this script we allow the submodels to
# have more complicated terms than simply bivariate

if(Sys.info()['sysname'] == "Darwin") {
  path <- "/Users/curculion/Documents/GitHub"
  path2 <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
  path <- "C:/Users/kerickson/Documents/GitHub"
  path2 <- "H:/Global Change Program/Research/ENMs - Modeling Methods for Rare Species (Kelley Erickson)/rare_species/data"
}


source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

modelType  <- models[4]

weights_df <- expand.grid(species, sizes, replicates)
names(weights_df) <- c("species", "sizes", "replicates")
for(i in 1:length(formulaMatrix$X)) {
  weights_df[,3+i] <- rep(0, length(weights_df$species))
}
names(weights_df)[4:29] <- paste0("m", 1:26)


timeStart <- Sys.time()

for( r in 1:100){
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
      names(train) <- c("Y", "PC1", "PC2", "PC3")
      presWeights <- rep(1/64, 64)
      contrastWeights <- rep(1/64, 64)

      if (n==1) {
        weight_vec <- NULL
        pred_sum <-rep(0, 0.5*length(Y_test))
        B0_sum <- rep(0, 0.5*length(Y_test))
        B1_sum <- rep(0, 0.5*length(Y_test))
        B2_sum <- rep(0, 0.5*length(Y_test))
        B3_sum <- rep(0, 0.5*length(Y_test))
        for (f in 1:6) {
          mod1 <- glm(as.formula(formulaMatrix$formula[f]),
                      family=binomial, data=train)
          pred <- predict(mod1, X_test,
                          type="response")


          auc_submod <- aucWeighted(pres=pred[1:64],
                                    contrast=pred[65:128],
                                    presWeight = presWeights,
                                    contrastWeight = contrastWeights)
          weight <- 2*auc_submod -1 #Schoners' D
          if(weight < 0) {weight <- 0}
          weight_vec <- c(weight_vec, weight)
          pred_sum <- pred_sum + (pred * weight)
          if(!is.na(formulaMatrix$Int[f]))   {B0_sum <- sum(B0_sum, weight*coefficients(mod1)[formulaMatrix$Int[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1[f]))   {B1_sum <- sum(B1_sum, weight*coefficients(mod1)[formulaMatrix$PC1[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2[f]))   {B2_sum <- sum(B2_sum, weight*coefficients(mod1)[formulaMatrix$PC2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3[f]))   {B3_sum <- sum(B3_sum, weight*coefficients(mod1)[formulaMatrix$PC3[f]], na.rm=T)}


        }
      }
      if (n>1) {
        weight_vec <- NULL
        pred_sum <-rep(0, 0.5*length(Y_test))
        B0_sum <- rep(0, 0.5*length(Y_test))
        B1_sum <- rep(0, 0.5*length(Y_test))
        B2_sum <- rep(0, 0.5*length(Y_test))
        B3_sum <- rep(0, 0.5*length(Y_test))
        B4_sum <- rep(0, 0.5*length(Y_test))
        B5_sum <- rep(0, 0.5*length(Y_test))
        B6_sum <- rep(0, 0.5*length(Y_test))
        for (f in 1:length(formulaMatrix$formula)) {
          mod1 <- glm(as.formula(formulaMatrix$formula[f]),
                      family=binomial, data=train)
          pred <- predict(mod1, X_test,
                          type="response")


          auc_submod <- aucWeighted(pres=pred[1:64],
                                    contrast=pred[65:128],
                                    presWeight = presWeights,
                                    contrastWeight = contrastWeights)
          weight <- 2*auc_submod -1 #Schoners' D
          if(weight < 0) {weight <- 0}
          weight_vec <- c(weight_vec, weight)
          pred_sum <- pred_sum + (pred * weight)
          if(!is.na(formulaMatrix$Int[f]))   {B0_sum <- sum(B0_sum, weight*coefficients(mod1)[formulaMatrix$Int[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1[f]))   {B1_sum <- sum(B1_sum, weight*coefficients(mod1)[formulaMatrix$PC1[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2[f]))   {B2_sum <- sum(B2_sum, weight*coefficients(mod1)[formulaMatrix$PC2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3[f]))   {B3_sum <- sum(B3_sum, weight*coefficients(mod1)[formulaMatrix$PC3[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC1_2[f])) {B4_sum <- sum(B4_sum, weight*coefficients(mod1)[formulaMatrix$PC1_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC2_2[f])) {B5_sum <- sum(B5_sum, weight*coefficients(mod1)[formulaMatrix$PC2_2[f]], na.rm=T)}
          if(!is.na(formulaMatrix$PC3_2[f])) {B6_sum <- sum(B6_sum, weight*coefficients(mod1)[formulaMatrix$PC3_2[f]], na.rm=T)}

        }
      }

      if(sum(weight_vec) == 0){
        pred_esm <- NA
        B0_esm <- NA
        B1_esm <- NA
        B2_esm <- NA
        B3_esm <- NA
        B4_esm <- NA
        B5_esm <- NA
        B6_esm <- NA
        aucWeighted_ESM <- NA
        RMSEWeighted_ESM <- NA
        PR_ESM <- NA
        TjurR2_ESM <- NA
        model_curve_PC1_ESM <- NA
        model_curve_PC2_ESM <- NA
        model_curve_PC3_ESM <- NA
        PC1_rankCor_ESM <- NA
        PC2_sd_ESM <- NA
        PC3_rankCor_ESM <- NA
        model$auc <- NA
        model$RMSEWeighted <- NA
        model$TjursR2 <- NA
        model$PC1_rankCor <- NA
        model$PC3_rankCor <- NA
        model$pr_integral <- NA
        model$varPC2Response <- NA

      }  else{

        # calculate ESM prediction
        pred_esm <- pred_sum/sum(weight_vec)


        #Calculate betas
        B0_esm <- B0_sum/sum(weight_vec)
        B1_esm <- B1_sum/sum(weight_vec)
        B2_esm <- B2_sum/sum(weight_vec)
        B3_esm <- B3_sum/sum(weight_vec)
        if(n >1) {
          B4_esm <- B4_sum/sum(weight_vec)
          B5_esm <- B5_sum/sum(weight_vec)
          B6_esm <- B6_sum/sum(weight_vec)}


        #evaluate ESM
        model$auc <- aucWeighted(pres=pred_esm[Y_test==1],
                                 contrast=pred_esm[Y_test==0],
                                 presWeight = presWeights,
                                 contrastWeight = contrastWeights)


        model$RMSEWeighted <- computeRMSEWeighted(Y_test, predY=pred_esm)

        model$pr_integral<- pr.curve(scores.class0=pred_esm,
                                     weights.class0 = Y_test,
                                     curve=TRUE)$auc.integral

        model$TjursR2 <- computeR2(Y=Y_test,
                                   predY=pred_esm)

        #Species Response Curves
        real_curve_PC1 <- dnorm(x1, mean=sims[[s]]$mu[1],
                                sd=sims[[s]]$sd[1,1])
        real_curve_PC3 <- dnorm(x3, mean=sims[[s]]$mu[2],
                                sd=sims[[s]]$sd[2,2])
        if(n ==1) {
          model_curve_PC1_ESM <- pnorm(B1_esm*x1)
          model_curve_PC2_ESM <- pnorm(B2_esm*x2)
          model_curve_PC3_ESM <- pnorm(B3_esm*x3)
        }
        if(n > 1){
          model_curve_PC1_ESM <- pnorm(B1_esm*x1 + B4_esm*x1*x1)
          model_curve_PC2_ESM <- pnorm(B2_esm*x2 + B5_esm*x2*x2)
          model_curve_PC3_ESM <- pnorm(B3_esm*x3 + B6_esm*x3*x3)
        }



        if(sd(model_curve_PC1_ESM)>0) {
          model$PC1_rankCor <- compareResponse(model_curve_PC1_ESM,
                                               real_curve_PC1,data =data.frame(x1),graph=F)$rankCor
        } else{model$PC1_rankCor <- NA}

        if(sd(model_curve_PC3_ESM)>0) {
          model$PC3_rankCor <- compareResponse(model_curve_PC3_ESM,
                                               real_curve_PC3, data =data.frame(x3), graph=F)$rankCor
        } else{model$PC3_rankCor <- NA}
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
            response$L[j] <- B0_esm + B1_esm*response$PC1[j] +
              B2_esm*response$PC2[j] + B3_esm*response$PC3[j]      }

        }
        if(n > 1) {
          for ( j in 1:length(response$PC1)) {
            response$L[j] <- B0_esm + B1_esm*response$PC1[j] +
              B2_esm*response$PC2[j] + B3_esm*response$PC3[j] +
              B4_esm*response$PC1[j]*response$PC1[j] +
              B5_esm*response$PC2[j]*response$PC2[j] +
              B6_esm*response$PC3[j]*response$PC3[j]
          }
        }

        model$varPC2Response<- sd(invlogit(response$L))


      }








      model$timeEnd <- Sys.time() - timeStart1
      save(model, file=paste0(path2, "/models/",
                              modelType[1], "/",species[s], "/",
                              sizes[n], "/", "model_",replicates[r],
                              ".RData"))

    }
  }
}
timeEnd <- Sys.time()-timeStart
