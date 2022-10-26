thing1 <- mvtnorm::dmvnorm(south[,c(3,5)],
            mean=mu_avg,
            sigma=var_narrow)

thing2 <- dnorm(south[,3], mean=mu_avg[1], sd=sd_narrow[1,1])*
  dnorm(south[,5], mean=mu_avg[2], sd=sd_narrow[2,2])

thing3 <- dnorm(south[,3], mean=mu_avg[1], sd=var_narrow[1,1])*
  dnorm(south[,5], mean=mu_avg[2], sd=var_narrow[2,2])

sum(thing1-thing2)
#There is a minor difference: -2.58*10^-15

#What about a difference in the predict function?



Y_train <- south[1:237,6:67]
Y_test <- south[238:474, 6:67]
Y <- south[,6:67]
XData_train <- south[1:237, 3:5]
XData_test <- south[238:474, 3:5]
XData <- south[,3:5]

partition <- c(rep(1, 237), rep(2,237))

thin <- 1
samples <- 500
transient <- 200
nChains <- 1
verbose <- 100
m <- Hmsc(Y=Y, XData=XData,
          XFormula = ~PC1 +PC2 + PC3 + I(PC1^2)+
            I(PC2^2) + I(PC3^2),
          distr="probit")

m_a <- sampleMcmc(m, thin=thin, samples=samples,
                transient=transient,
                nChains=nChains, verbose=verbose,
                initPar = "fixed effects")
m_b <- sampleMcmc(m, thin=thin, samples=samples,
                 transient=transient, nChains=nChains,
                 verbose=verbose)



preds_a <- computePredictedValues(m_a,
                                 partition=partition)
preds_b <- computePredictedValues(m_b, partition=partition)

predsa_sp <- preds_a[,62,] #dims 474 x 500
predsb_sp <- preds_b[,62,]



predsa_mean <- apply(predsa_sp, 1, "mean") #length 474
predsb_mean <- apply(predsb_sp, 1, "mean")


NPresent <- sum(Y[,62])
NAbsent <- 474 - NPresent
NPresentfold1 <- sum(Y[1:237,62])
NAbsent1 <- 237-NPresent1
NPresent2 <- sum(Y[238:474,62])
NAbsent2 <- 237-NPresent2
presWeights <- rep(1/NPresent, NPresent)
presWeights1 <- rep(1/NPresent1, NPresent1)
presWeights2 <- rep(1/NPresent2, NPresent2)
contrastWeights1 <- rep(1/NAbsent1, NAbsent1)
contrastWeights2 <- rep(1/NAbsent2, NAbsent2)
contrastWeights <- rep(1/NAbsent, NAbsent)
Y1 <- Y[partition==1,62]
Y2 <- Y[partition==2, 62]

predsa_mean_fold1 <- predsa_mean[1:237]
predsa_mean_fold2 <- predsa_mean[238:474]

predsb_mean_fold1 <- predsb_mean[1:237]
predsb_mean_fold2 <- predsb_mean[238:474]

predsa_mean_present <- predsa_mean[Y[,62]==1]
predsa_mean_absent <- predsa_mean[Y[,62]==0]
predsa_mean_fold1_present <- predsa_mean_fold1[Y1==1]
predsa_mean_fold1_absent  <- predsa_mean_fold1[Y1==0]

predsa_mean_fold2_present <- predsa_mean_fold2[Y2==1]
predsa_mean_fold2_absent  <- predsa_mean_fold2[Y2==0]

predsb_mean_present <- predsb_mean[Y[,62]==1]
predsb_mean_absent <- predsb_mean[Y[,62]==0]
predsb_mean_fold1_present <- predsb_mean_fold1[Y1==1]
predsb_mean_fold1_absent  <- predsb_mean_fold1[Y1==0]

predsb_mean_fold2_present <- predsb_mean_fold2[Y2==1]
predsb_mean_fold2_absent  <- predsb_mean_fold2[Y2==0]

auc_a <- aucWeighted(pres=predsa_mean_present,
                    contrast=predsa_mean_absent,
                    presWeight = presWeights,
                    contrastWeight = contrastWeights)
auc_b <- aucWeighted(pres=predsb_mean_present,
                     contrast =predsb_mean_absent,
                     presWeight <- presWeights,
                     contrastWeight = contrastWeights)

#Now to do it the new way
preds <- predict(m_a,XData=XData_test, expected=F) #expected=T, probs, expected =F, 0,1

thing <- lapply(
  X = preds,
  FUN = function(x){x[,62]}
)

Matrix_x <- matrix(unlist(thing), ncol = nChains*samples, byrow = TRUE)
rm(thing)

preds.mean <- apply(Matrix_x, FUN="mean", MARGIN=1)

predsa_mean[1:237]-preds.mean
auc_new_a <- aucWeighted(pres=preds.mean[Y[1:237,62] == 1],
            contrast=preds.mean[Y[1:237,62]==0],
            presWeight = presWeights1,
            contrastWeight = contrastWeights1)

auc_new_b <- aucWeighted(pres=pres.mean[Y[238:474,62] == 1],
                         contrast=preds.mean[Y[238:474,62] ==2])


#### New test #####

#computePredictedValues
preds_a <- computePredictedValues(m_a,
                                  partition=partition)
#Dimensions 474 sites x 62 species x 500 iterations

#Select just species 62
predsa_sp <- preds_a[,62,]

#Average over the iterations
predsa_mean <- apply(predsa_sp, 1, "mean") #length 474

NPresent <- sum(Y[,62])
NAbsent <- 474 - NPresent
#NPresentfold1 <- sum(Y[1:237,62])
#NAbsent1 <- 237-NPresent1
#NPresent2 <- sum(Y[238:474,62])
#NAbsent2 <- 237-NPresent2
presWeights <- rep(1/NPresent, NPresent)
#presWeights1 <- rep(1/NPresent1, NPresent1)
#presWeights2 <- rep(1/NPresent2, NPresent2)
#contrastWeights1 <- rep(1/NAbsent1, NAbsent1)
#contrastWeights2 <- rep(1/NAbsent2, NAbsent2)
contrastWeights <- rep(1/NAbsent, NAbsent)

Y1 <- Y[partition==1,62]
Y2 <- Y[partition==2, 62]

predsa_mean_present <- predsa_mean[Y[,62]==1]
predsa_mean_absent <- predsa_mean[Y[,62]==0]


predsa_mean_fold1 <- predsa_mean[1:237]
predsa_mean_fold2 <- predsa_mean[238:474]

aucWeighted(pres=predsa_mean_present,
            contrast=predsa_mean_absent,
            presWeight = presWeights,
            contrastWeight = contrastWeights)

#Now do for Hmsc.predict

preds <- predict(m_a,XDataNew=XData,
                 expected=T,
                 predictEtaMean=T,
                 mcmcStep=1000) #expected=T, probs, expected =F, 0,1

thing <- lapply(
  X = preds,
  FUN = function(x){x[,62]}
)

Matrix_x <- matrix(unlist(thing), ncol = nChains*samples, byrow = TRUE)
rm(thing)

preds.mean <- apply(Matrix_x, FUN="mean", MARGIN=1)
preds.mean_present <- preds.mean[Y[,62]==1]
preds.mean_absent <- preds.mean[Y[,62]==0]

aucWeighted(pres=preds.mean_present,
            contrast=preds.mean_absent,
            presWeight = presWeights,
            contrastWeight = contrastWeights)

gradient <- prepareGradient(m_a, XDataNew = XData)
preds <- predict(m_a, Gradient=gradient, expected=T)


#From function
postList <- poolMcmcChains(m_a$postList)
#500
preds <- predict(m_a, postList, expected=T)
pred_array <- abind(preds, along=3)