#02_Model_6_SAM

#In this script we run SAM model, using the
# number of aggregates identified in script
# `02_Model_6a_Real_Species_SAM`
if(Sys.info()['sysname'] == "Darwin") {
  path <- "/Users/curculion/Documents/GitHub"
  path2 <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
  path <- "C:/Users/kerickson/Documents/GitHub"
  path2 <- "H:/Global Change Program/Research/ENMs - Modeling Methods for Rare Species (Kelley Erickson)/rare_species/data"
}


source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

modelType  <- models[6]

timeStart <- Sys.time()

for( r in 1:length(replicates)){
  for (s in 1:length(species)){
    for (n in 1:length(sizes)){
      if(!file.exists(paste0(path2, "/models/",
                             modelType[1], "/",species[s], "/",
                             sizes[n], "/", "model_",replicates[r],
                             ".RData"))) {

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
      train <- cbind(Y_train, X_train)

      archetype_formula <- as.formula(paste0(paste0('cbind(',paste(colnames(Y_train),collapse = ", "),") ~ PC1 + PC2 + PC3 + I(PC1^2) + I(PC2^2) + I(PC3^2)")))

      ## Species formula
      species_formula <- ~ 1


      val <- tryCatch({m <- species_mix(archetype_formula = archetype_formula, # Archetype formula
                       species_formula = species_formula,    # Species formula
                       data = train,            # Data
                       nArchetypes = 6,              # Number of groups (mixtures) to fit
                       family = 'bernoulli',
                       titbits = TRUE,
                       control = list(maxits = 100000))
      model$preds <- predict(m, newdata=X_test,
                      type="response", prediction.type="species")[,"SimSp"]

      presWeights <- rep(1/64, 64)
      contrastWeights <- rep(1/64, 64)


      model$auc <- aucWeighted(pres=model$preds[1:64],
                               contrast=model$preds[65:128],
                               presWeight = presWeights,
                               contrastWeight = contrastWeights)

      model$RMSEWeighted <- computeRMSEWeighted(Y_test["SimSp"], model$preds)
      model$TjursR2 <- computeR2(Y_test["SimSp"], model$preds)


      b1 <- m$beta[1]*m$tau["SimSp",1] +
        m$beta[2]*m$tau["SimSp", 2] +
        m$beta[3]*m$tau["SimSp", 3] +
        m$beta[4]*m$tau["SimSp", 4] +
        m$beta[5]*m$tau["SimSp",5] +
        m$beta[6]*m$tau["SimSp", 6]

      b2 <- m$beta[7]*m$tau["SimSp",1] +
        m$beta[8]*m$tau["SimSp",2] +
        m$beta[9]*m$tau["SimSp",3] +
        m$beta[10]*m$tau["SimSp", 4] +
        m$beta[11]*m$tau["SimSp", 5] +
        m$beta[12]*m$tau["SimSp",6]

      b3 <- m$beta[13]*m$tau["SimSp",1] +
        m$beta[14]*m$tau["SimSp",2] +
        m$beta[15]*m$tau["SimSp",3] +
        m$beta[16]*m$tau["SimSp",4] +
        m$beta[17]*m$tau["SimSp", 5] +
        m$beta[18]*m$tau["SimSp", 6]

      b4 <- m$beta[19]*m$tau["SimSp",1] +
        m$beta[20]*m$tau["SimSp",2] +
        m$beta[21]*m$tau["SimSp",3] +
        m$beta[22]*m$tau["SimSp",4] +
        m$beta[23]*m$tau["SimSp", 5] +
        m$beta[24]*m$tau["SimSp", 6]

      b5 <- m$beta[25]*m$tau["SimSp",1] +
        m$beta[26]*m$tau["SimSp",2] +
        m$beta[27]*m$tau["SimSp",3] +
        m$beta[28]*m$tau["SimSp",4] +
        m$beta[29]*m$tau["SimSp", 5] +
        m$beta[30]*m$tau["SimSp", 6]

      b6 <- m$beta[31]*m$tau["SimSp",1] +
        m$beta[32]*m$tau["SimSp",2] +
        m$beta[33]*m$tau["SimSp",3] +
        m$beta[34]*m$tau["SimSp",4] +
        m$beta[35]*m$tau["SimSp", 5] +
        m$beta[36]*m$tau["SimSp", 6]

        model_curve_PC1 <- pnorm(b1*x1 + b4*x1*x1)
        model_curve_PC2 <- pnorm(b2*x2 + b5*x2*x2)
        model_curve_PC3 <- pnorm(b3*x3 + b6*x3*x3)

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
          response$L[j] <- m$alpha["SimSp"] + b1*response$PC1[j] +
            b2*response$PC2[j] + b3*response$PC3[j] +
            b4*response$PC1[j]*response$PC1[j] +
            b5*response$PC2[j]*response$PC2[j] +
            b6*response$PC3[j]*response$PC3[j]
        }


      model$varPC2Response<- sd(invlogit(response$L))

      #Precision Recall
      PR <- pr.curve(scores.class0=model$preds,
                     weights.class0 = Y_test[["SimSp"]],
                     curve=FALSE)
      model$pr_integral <- PR$auc.integral}, warning=function(w) w, error=function(w) w)
      if(class(val)[1] == "simpleWarning" | class(val)[1]=="simpleError"){
        model$auc <- NA
        model$RMSEWeighted <- NA
        model$TjursR2 <- NA
        model$PC1_cor <- NA
        model$PC3_cor <- NA
        model$P1_rankCor <- NA
        model$PC3_rankCor <- NA
        model$varPC2Response <- NA
        model$pr_integral <- NA
      }


      model$timeEnd <- Sys.time() - timeStart1
   save(model, file=paste0(path2, "/models/",
                              modelType[1], "/",species[s], "/",
                              sizes[n], "/", "model_",replicates[r],
                              ".RData"))}

    }
  }
}
timeEnd <- Sys.time()-timeStart




