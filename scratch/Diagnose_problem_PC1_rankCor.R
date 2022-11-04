#Check Rank Correlation Calculations

s <- 1
n <- 2

model_curve_PC1 <- list()
real_curve_PC1 <- dnorm(x1, mean=sims[[s]]$mu[1],
                        sd=sims[[s]]$sd[1,1])

orderOfX <- order(x1)
x <- x1[orderOfX]
pred2 <- real_curve_PC1[orderOfX]
pred1 <- list()


for (r in 1:length(replicates)){
  model <-list()
  model$replicate <- replicates[r]
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


  if (n>1) {
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

    #Species Response Curves
    model_curve_PC1[[r]] <- pnorm(B1_esm*x1)
    pred1[[r]] <- model_curve_PC1[[r]][orderOfX]
  }}


sim <- compareNiches(x1=pred1, x2=pred2, logit=!adjust, na.rm=T)


