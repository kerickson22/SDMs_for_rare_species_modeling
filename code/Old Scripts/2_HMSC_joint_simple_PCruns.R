# Script for running all of the HMSC-variants
library(Hmsc)
#library(snow)
library(MCMCvis)
library(enmSdm)
library(omnibus)
memory.limit(2^32)
load("../data/south.RData")
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

load("../data/X_bar.RData")
mu_ext <- c(-3.5, 2)
Sigma_narrow <- diag(c(0.9, 2.5))

mu_avg <- c(X_bar$mu_PC1[64], X_bar$mu_PC3[64])
Sigma_wide <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.90, na.rm=T)^2,  quantile(X_bar$sd_PC3, probs=0.90, na.rm=T)^2))

#Load simulated data
load(paste0("../data/sims_data.RData"))


computeRMSEWeighted = function(Y, predY) {
  RMSE <- sqrt((mean((Y[Y==1]-predY[Y==1])^2, na.rm=TRUE) +
                  mean((Y[Y==0]-predY[Y==0])^2, na.rm=TRUE)))
  return(RMSE)
}



nChains <- 2
thin <- 1
samples <- 10000
#samples <- 10
transient <-2000
#transient <- 5
verbose <-1000
#verbose <- 2
runif(1)

timeStart <- Sys.time()
s <- 4
n <- 1
for( r in 7:18){
      if(file.exists(paste0("../data/models/",
                         modelType[1], "/",species[s], "/",
                         sizes[n], "/", "model_",replicates[r],
                         ".RData"))){print("already done")
        }else{
      timeStart1 <- Sys.time()
        model <-list()
        model$seed <- .GlobalEnv$.Random.seed
        model$NPresesences <- sizes[n]
        model$species <- species[s]
        model$replicate <- replicates[r]
        model$nChains <- nChains
        model$samples <- samples
        model$transient <- transient
        model$thin <-thin
        model$modelType <- modelType[1]




        Y <- cbind(south[,6:68],as.matrix(sims[[s]][[n]][[r]]$YSim))
        XData <- south[,3:5]
        m <- Hmsc(Y=Y, XData=XData,
                        XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                          I(PC2^2)+I(PC3^2),
                        distr="probit")
        m <- sampleMcmc(m, thin=thin, samples=samples,
                              transient=transient,
                              nChains=nChains, verbose=verbose)

        partition1 <- sims[[s]][[n]][[r]]$split1
        partition2 <- sims[[s]][[n]][[r]]$split2
        partition3 <- sims[[s]][[n]][[r]]$split3

        preds1 <- computePredictedValues(m,
                                        partition=partition1)
        preds1_sp <- preds1[,64,]
        model$preds1_mean <- apply(preds1_sp, 1, "mean")
        preds2 <- computePredictedValues(m,
                                         partition=partition2)
        preds2_sp <- preds2[,64,]
        model$preds2_mean <- apply(preds2_sp, 1, "mean")

         preds3 <- computePredictedValues(m,
                                         partition=partition3)
        preds3_sp <- preds3[,64,]
        model$preds3_mean <- apply(preds3_sp, 1, "mean")


        NPresent <- sum(Y[,64]/2)
        NAbsent <- (dim(Y)[1]/2)-NPresent
        presWeights <- rep(1/NPresent, NPresent)
        contrastWeights <- rep(1/NAbsent, NAbsent)

        #Split 1
        Y1_1 <- Y[partition1==1,64]
        Y1_2 <- Y[partition1==2, 64]

        preds1_mean1 <- model$preds1_mean[partition1==1]
        preds1_mean2 <- model$preds1_mean[partition1==2]



        preds1_mean1_present <- preds1_mean1[Y1_1==1]
        preds1_mean1_absent  <- preds1_mean1[Y1_1==0]

        preds1_mean2_present <- preds1_mean2[Y1_2==1]
        preds1_mean2_absent  <- preds1_mean2[Y1_2==0]


        model$aucWeighted_1_1 <- aucWeighted(pres=preds1_mean1_present,
                              contrast=preds1_mean1_absent,
                              presWeight = presWeights,
                              contrastWeight = contrastWeights)
        model$aucWeighted_1_2 <- aucWeighted(pres=preds1_mean2_present,
                              contrast=preds1_mean2_absent,
                              presWeight = presWeights,
                              contrastWeight = contrastWeights)
        model$RMSEWeighted_1_1 <- computeRMSEWeighted(Y1_1, preds1_mean1)
        model$RMSEWeighted_1_2 <- computeRMSEWeighted(Y1_2, preds1_mean2)
        

        #Split 2
        Y2_1 <- Y[partition1==1,64]
        Y2_2 <- Y[partition1==2, 64]
        model$preds2_mean <- apply(preds2, 1, "mean")
        preds2_mean1 <- model$preds2_mean[partition2==1]
        preds2_mean2 <- model$preds2_mean[partition2==2]

        preds2_mean1_present <- preds2_mean1[Y2_1==1]
        preds2_mean1_absent  <- preds2_mean1[Y2_1==0]

        preds2_mean2_present <- preds2_mean2[Y2_2==1]
        preds2_mean2_absent  <- preds2_mean2[Y2_2==0]


        model$aucWeighted_2_1 <- aucWeighted(pres=preds2_mean1_present,
                                             contrast=preds2_mean1_absent,
                                             presWeight = presWeights,
                                             contrastWeight = contrastWeights)
        model$aucWeighted_2_2 <- aucWeighted(pres=preds2_mean2_present,
                                             contrast=preds2_mean2_absent,
                                             presWeight = presWeights,
                                             contrastWeight = contrastWeights)
        model$RMSEWeighted_2_1 <- computeRMSEWeighted(Y2_1, preds2_mean1)
        model$RMSEWeighted_2_2 <- computeRMSEWeighted(Y2_2, preds2_mean2)
        
        #Split 3
        Y3_1 <- Y[partition1==1,64]
        Y3_2 <- Y[partition1==2, 64]
        model$preds3_mean <- apply(preds3, 1, "mean")
        preds3_mean1 <- model$preds3_mean[partition3==1]
        preds3_mean2 <- model$preds3_mean[partition3==2]

        preds3_mean1_present <- preds3_mean1[Y3_1==1]
        preds3_mean1_absent  <- preds3_mean1[Y3_1==0]

        preds3_mean2_present <- preds3_mean2[Y3_2==1]
        preds3_mean2_absent  <- preds3_mean2[Y3_2==0]


        model$aucWeighted_3_1 <- aucWeighted(pres=preds3_mean1_present,
                                             contrast=preds3_mean1_absent,
                                             presWeight = presWeights,
                                             contrastWeight = contrastWeights)
        model$aucWeighted_3_2 <- aucWeighted(pres=preds3_mean2_present,
                                             contrast=preds3_mean2_absent,
                                             presWeight = presWeights,
                                             contrastWeight = contrastWeights)

        model$RMSEWeighted_3_1 <- computeRMSEWeighted(Y3_1, preds3_mean1)
        model$RMSEWeighted_3_2 <- computeRMSEWeighted(Y3_2, preds3_mean2)
        
        model$split1_metrics <- evaluateModelFit(m, preds1)
        model$split2_metrics <- evaluateModelFit(m, preds2)
        model$split3_metrics <- evaluateModelFit(m, preds3)


       m1_1 <-Hmsc(Y=Y[partition1==1,], XData=XData[partition1==1,],
                                    XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                                      I(PC2^2)+I(PC3^2),
                                    distr="probit")
       m1_2 <-Hmsc(Y=Y[partition1==2,], XData=XData[partition1==2,],
                   XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                     I(PC2^2)+I(PC3^2),
                   distr="probit")
       m2_1 <-Hmsc(Y=Y[partition2==1,], XData=XData[partition2==1,],
                   XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                     I(PC2^2)+I(PC3^2),
                   distr="probit")
       m2_2 <-Hmsc(Y=Y[partition2==2,], XData=XData[partition2==2,],
                   XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                     I(PC2^2)+I(PC3^2),
                   distr="probit")
       m3_1 <-Hmsc(Y=Y[partition3==1,], XData=XData[partition3==1,],
                   XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                     I(PC2^2)+I(PC3^2),
                   distr="probit")
       m3_2 <-Hmsc(Y=Y[partition3==2,], XData=XData[partition3==2,],
                   XFormula = ~PC1 + PC2 + PC3 + I(PC1^2)+
                     I(PC2^2)+I(PC3^2),
                   distr="probit")
       m1_1 <- sampleMcmc(m1_1, thin=thin, samples=samples,
                             transient=transient,
                             nChains=nChains, verbose=verbose)
       m1_2 <- sampleMcmc(m1_2, thin=thin, samples=samples,
                                transient=transient, 
                                nChains=nChains, verbose=verbose)
       m2_1 <- sampleMcmc(m2_1, thin=thin, samples=samples,
                                transient=transient, 
                                nChains=nChains, verbose=verbose)
       m2_2 <- sampleMcmc(m2_2, thin=thin, samples=samples,
                                transient=transient, 
                                nChains=nChains, verbose=verbose)
       m3_1 <- sampleMcmc(m3_1, thin=thin, samples=samples,
                                transient=transient, 
                                nChains=nChains, verbose=verbose)
       m3_2 <- sampleMcmc(m1_1, thin=thin, samples=samples,
                                transient=transient,
                                nChains=nChains, verbose=verbose)

       m.post1_1 <- convertToCodaObject(m1_1)
       m.post1_2 <- convertToCodaObject(m1_2)
       m.post2_1 <- convertToCodaObject(m2_1)
       m.post2_2 <- convertToCodaObject(m2_2)
       m.post3_1 <- convertToCodaObject(m3_1)
       m.post3_2 <- convertToCodaObject(m3_2)

      #Split 1
       betas1_1 <- MCMCsummary(m.post1_1$Beta)
       betas1_2 <- MCMCsummary(m.post1_2$Beta)
       model$beta_chain1_1 <- m.post1_1$Beta
       model$beta_chain1_2 <- m.post1_2$Beta
  

       #Split2
       betas2_1 <- MCMCsummary(m.post2_1$Beta)
       betas2_2 <- MCMCsummary(m.post2_2$Beta)
       model$beta_chain2_1 <- m.post2_1$Beta
       model$beta_chain2_2 <- m.post2_2$Beta


       #Split3
       betas3_1 <- MCMCsummary(m.post3_1$Beta)
       betas3_2 <- MCMCsummary(m.post3_2$Beta)
       model$beta_chain3_1 <- m.post3_1$Beta
       model$beta_chain3_2 <- m.post3_2$Beta
    
       model$maxRHat <- max(betas1_1$Rhat, betas1_2$Rhat,
                            betas2_1$Rhat, betas2_2$Rhat,
                            betas3_1$Rhat, betas3_2$Rhat,
                            na.rm=T)


       model_curve_PC1_m1_1 <- pnorm(betas1_1$mean[443]*x1 + betas1_1$mean[446]*x1*x1)
       model_curve_PC2_m1_1 <- pnorm(betas1_1$mean[444]*x2 + betas1_1$mean[447]*x2*x2)
       model_curve_PC3_m1_1 <- pnorm(betas1_1$mean[445]*x3 + betas1_1$mean[448]*x3*x3)
       model_curve_PC1_m1_2 <- pnorm(betas1_2$mean[443]*x1 + betas1_2$mean[446]*x1*x1)
       model_curve_PC2_m1_2 <- pnorm(betas1_2$mean[444]*x2 + betas1_2$mean[447]*x2*x2)
       model_curve_PC3_m1_2 <- pnorm(betas1_2$mean[445]*x3 + betas1_2$mean[448]*x3*x3)

       model_curve_PC1_m2_1 <- pnorm(betas2_1$mean[443]*x1 + betas2_1$mean[446]*x1*x1)
       model_curve_PC2_m2_1 <- pnorm(betas2_1$mean[444]*x2 + betas2_1$mean[447]*x2*x2)
       model_curve_PC3_m2_1 <- pnorm(betas2_1$mean[445]*x3 + betas2_1$mean[448]*x3*x3)
       model_curve_PC1_m2_2 <- pnorm(betas2_2$mean[443]*x1 + betas2_2$mean[446]*x1*x1)
       model_curve_PC2_m2_2 <- pnorm(betas2_2$mean[444]*x2 + betas2_2$mean[447]*x2*x2)
       model_curve_PC3_m2_2 <- pnorm(betas2_2$mean[445]*x3 + betas2_2$mean[448]*x3*x3)

       model_curve_PC1_m3_1 <- pnorm(betas3_1$mean[443]*x1 + betas3_1$mean[446]*x1*x1)
       model_curve_PC2_m3_1 <- pnorm(betas3_1$mean[444]*x2 + betas3_1$mean[447]*x2*x2)
       model_curve_PC3_m3_1 <- pnorm(betas3_1$mean[445]*x3 + betas3_1$mean[448]*x3*x3)
       model_curve_PC1_m3_2 <- pnorm(betas3_2$mean[443]*x1 + betas3_2$mean[446]*x1*x1)
       model_curve_PC2_m3_2 <- pnorm(betas3_2$mean[444]*x2 + betas3_2$mean[447]*x2*x2)
       model_curve_PC3_m3_2 <- pnorm(betas3_2$mean[445]*x3 + betas3_2$mean[448]*x3*x3)

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

       model$PC1_cor<- mean(compareResponse(model_curve_PC1_m1_1,
                                                  real_curve_PC1,
                                                  data =data.frame(x1),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC1_m1_2,
                                                  real_curve_PC1,
                                                  data =data.frame(x1),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC1_m2_1,
                                                  real_curve_PC1,
                                                  data =data.frame(x1),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC1_m2_2,
                                                  real_curve_PC1,
                                                  data =data.frame(x1),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC1_m3_1,
                                                  real_curve_PC1,
                                                  data =data.frame(x1),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC1_m3_2,
                                                  real_curve_PC1,
                                                  data =data.frame(x1),
                                                  graph=F)$cor
       )

       model$PC3_cor <- mean(compareResponse(model_curve_PC3_m1_1,
                                                  real_curve_PC3,
                                                  data =data.frame(x3),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC3_m1_2,
                                                  real_curve_PC3,
                                                  data =data.frame(x3),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC3_m2_1,
                                                  real_curve_PC3,
                                                  data =data.frame(x3),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC3_m2_2,
                                                  real_curve_PC3,
                                                  data =data.frame(x3),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC3_m3_1,
                                                  real_curve_PC3,
                                                  data =data.frame(x3),
                                                  graph=F)$cor,
                                  compareResponse(model_curve_PC3_m3_2,
                                                  real_curve_PC3,
                                                  data =data.frame(x3),
                                                  graph=F)$cor
       )

       model$PC1_rankCor <- mean(compareResponse(model_curve_PC1_m1_1,
                                                      real_curve_PC1,
                                                      data =data.frame(x1),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC1_m1_2,
                                                      real_curve_PC1,
                                                      data =data.frame(x1),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC1_m2_1,
                                                      real_curve_PC1,
                                                      data =data.frame(x1),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC1_m2_2,
                                                      real_curve_PC1,
                                                      data =data.frame(x1),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC1_m3_1,
                                                      real_curve_PC1,
                                                      data =data.frame(x1),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC1_m3_2,
                                                      real_curve_PC1,
                                                      data =data.frame(x1),
                                                      graph=F)$rankCor
       )

       model$PC3_rankCor <- mean(compareResponse(model_curve_PC3_m1_1,
                                                      real_curve_PC3,
                                                      data =data.frame(x3),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC3_m1_2,
                                                      real_curve_PC3,
                                                      data =data.frame(x3),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC3_m2_1,
                                                      real_curve_PC3,
                                                      data =data.frame(x3),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC3_m2_2,
                                                      real_curve_PC3,
                                                      data =data.frame(x3),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC3_m3_1,
                                                      real_curve_PC3,
                                                      data =data.frame(x3),
                                                      graph=F)$rankCor,
                                      compareResponse(model_curve_PC3_m3_2,
                                                      real_curve_PC3,
                                                      data =data.frame(x3),
                                                      graph=F)$rankCor
       )
        model$timeEnd <- Sys.time() - timeStart1
       save(model, file=paste0("../data/models/",
          modelType[1], "/",species[s], "/",
          sizes[n], "/", "model_",replicates[r],
          ".RData"))
      }
      }
    

  
timeEnd <- Sys.time()-timeStart
save(timeEnd, file=paste0("../data/models/",
                          modelType[1], "/time.RData"))