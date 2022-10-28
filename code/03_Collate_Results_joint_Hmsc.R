#03_Collate_Results

#In this script, we assemble all of the results
# from the individual models into one data frame

if(Sys.info()['sysname'] == "Darwin") {
  path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
  path <- "C:/Users/kerickson/Documents/GitHub"
  path2 <- "H:/Global Change Program/Research/ENMs - Modeling Methods for Rare Species (Kelley Erickson)/rare_species/data"
}


source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

#Joint_Hmsc #####
modelType  <- models[4]
results <- expand.grid(species, sizes, replicates)
results <- results[order(results$Var1, results$Var2, results$Var3),]
names(results) <- c("species", "size", "rep")

results$auc <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat<- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$varPC2Response <- rep(NA, length(results$species))

for (i in 1:length(results$species)) {
  if(file.exists(paste0(path2, "/models/",
                        modelType[1], "/",results$species[i], "/",
                        results$size[i], "/", "model_",results$rep[i],
                        ".RData"))) {
    load(paste0(path2, "/models/",
                modelType[1], "/",results$species[i], "/",
                results$size[i], "/", "model_",results$rep[i],
                ".RData"))

    results$auc[i] <- model$auc
    results$RMSEWeighted[i] <- model$RMSEWeighted
    results$TjursR2[i] <- model$TjursR2
    results$PC1_rankCor[i] <- model$PC1_rankCor
    results$PC3_rankCor[i] <- model$PC3_rankCor
    results$Rhat[i]<- model$maxRhat
    results$pr_integral[i] <- model$pr_integral
    results$varPC2Response[i] <- model$varPC2Response
  }
}

save(results, file=paste0(path2, "/models/",
                          modelType[1], "/results.RData" ))
myfile <- paste0(path2, "/models/",
                 modelType[1], "/results.RData" ) %>%
  drive_upload(paste0("status_updates_for_Hmsc_joint/Hmsc_joint_results.RData"))
