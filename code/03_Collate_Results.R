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

# Single-species glm #####

#glm should only be used for n>2
modelType  <- models[1]
results <- expand.grid(species, sizes[3:6], replicates)
results <- results[order(results$Var1, results$Var2, results$Var3), ]
names(results) <- c("species", "size", "rep")
results$model <- rep(modelType, length(results$species))
results$auc <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat <- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$varPC2Response <- rep(NA, length(results$species))

for (i in 1:length(results$species)) {
  if (file.exists(paste0(path2, "/models/", modelType[1], "/", results$species[i], "/", results$size[i], "/", "model_", results$rep[i], ".RData"))) {
    load(
      paste0(
        path2,
        "/models/",
        modelType[1],
        "/",
        results$species[i],
        "/",
        results$size[i],
        "/",
        "model_",
        results$rep[i],
        ".RData"
      )
    )

    results$auc[i] <- model$auc
    results$RMSEWeighted[i] <- model$RMSEWeighted
    results$TjursR2[i] <- model$TjursR2
    results$PC1_rankCor[i] <- model$PC1_rankCor
    results$PC3_rankCor[i] <- model$PC3_rankCor
    results$pr_integral[i] <- model$pr_integral
    results$varPC2Response[i] <- model$varPC2Response
  }
}

save(results, file=paste0(path2, "/models/",
                          modelType[1], "/results.RData" ))

# ESM_linear #####

modelType  <- models[2]
results <- expand.grid(species, sizes, replicates)
results <- results[order(results$Var1, results$Var2, results$Var3), ]
names(results) <- c("species", "size", "rep")

results$auc <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat <- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$varPC2Response <- rep(NA, length(results$species))

for (i in 1:length(results$species)) {
  if (file.exists(paste0(path2, "/models/", modelType[1], "/", results$species[i], "/", results$size[i], "/", "model_", results$rep[i], ".RData"))) {
    load(
      paste0(
        path2,
        "/models/",
        modelType[1],
        "/",
        results$species[i],
        "/",
        results$size[i],
        "/",
        "model_",
        results$rep[i],
        ".RData"
      )
    )

    results$auc[i] <- model$auc
    results$RMSEWeighted[i] <- model$RMSEWeighted
    results$TjursR2[i] <- model$TjursR2
    results$PC1_rankCor[i] <- model$PC1_rankCor
    results$PC3_rankCor[i] <- model$PC3_rankCor
    results$pr_integral[i] <- model$pr_integral
    results$varPC2Response[i] <- model$varPC2Response
  }
}

save(results, file=paste0(path2, "/models/",
                          modelType[1], "/results.RData" ))

# ESM_complex #####
#ESM complex should only be used for n>=8
modelType  <- models[3]
results <- expand.grid(species, sizes[3:6], replicates)
results <- results[order(results$Var1, results$Var2, results$Var3), ]
names(results) <- c("species", "size", "rep")

results$auc <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat <- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$varPC2Response <- rep(NA, length(results$species))

for (i in 1:length(results$species)) {
  if (file.exists(paste0(path2, "/models/", modelType[1], "/", results$species[i], "/", results$size[i], "/", "model_", results$rep[i], ".RData"))) {
    load(
      paste0(
        path2,
        "/models/",
        modelType[1],
        "/",
        results$species[i],
        "/",
        results$size[i],
        "/",
        "model_",
        results$rep[i],
        ".RData"
      )
    )

    results$auc[i] <- model$auc
    results$RMSEWeighted[i] <- model$RMSEWeighted
    results$TjursR2[i] <- model$TjursR2
    results$PC1_rankCor[i] <- model$PC1_rankCor
    results$PC3_rankCor[i] <- model$PC3_rankCor
    results$pr_integral[i] <- model$pr_integral
    results$varPC2Response[i] <- model$varPC2Response
  }
}

save(results, file=paste0(path2, "/models/",
                          modelType[1], "/results.RData" ))

# Joint Hmsc ####

# SAM #####
modelType  <- models[5]
results <- expand.grid(species, sizes, replicates)
results <- results[order(results$Var1, results$Var2, results$Var3), ]
names(results) <- c("species", "size", "rep")

taus <- expand.grid(species, sizes, replicates)
taus <- taus[order(taus$Var1, taus$Var2, taus$Var3), ]
names(taus) <- c("species", "size", "rep")

results$auc <- rep(NA, length(results$species))
results$RMSEWeighted <- rep(NA, length(results$species))
results$TjursR2 <- rep(NA, length(results$species))
results$PC1_rankCor <- rep(NA, length(results$species))
results$PC3_rankCor <- rep(NA, length(results$species))
results$Rhat <- rep(NA, length(results$species))
results$pr_integral <- rep(NA, length(results$species))
results$varPC2Response <- rep(NA, length(results$species))

taus$Archetype1 <- rep(NA, length(taus$species))
taus$Archetype2 <- rep(NA, length(taus$species))
taus$Archetype3 <- rep(NA, length(taus$species))
taus$Archetype4 <- rep(NA, length(taus$species))
taus$Archetype5 <- rep(NA, length(taus$species))
taus$Archetype6 <- rep(NA, length(taus$species))

for (i in 1:length(results$species)) {
  if (file.exists(paste0(path2, "/models/", modelType[1], "/", results$species[i], "/", results$size[i], "/", "model_", results$rep[i], ".RData"))) {
    load(
      paste0(
        path2,
        "/models/",
        modelType[1],
        "/",
        results$species[i],
        "/",
        results$size[i],
        "/",
        "model_",
        results$rep[i],
        ".RData"
      )
    )

    results$auc[i] <- model$auc
    results$RMSEWeighted[i] <- model$RMSEWeighted
    results$TjursR2[i] <- model$TjursR2
    if (names(model)[12] =="P1_rankCor")  {results$PC1_rankCor[i] <- model$P1_rankCor} #fix typo from original
    if (names(model)[12] =="PC1_rankCor") {results$PC1_rankCor[i] <- model$PC1_rankCor}
    results$PC3_rankCor[i] <- model$PC3_rankCor
    results$pr_integral[i] <- model$pr_integral
    results$varPC2Response[i] <- model$varPC2Response
    if(!is.null(model$tau)) {taus[i,4:9] <- model$tau}
    }
}

save(results, file=paste0(path2, "/models/",
                          modelType[1], "/results.RData" ))
save(taus, file=paste0(path2, "/models/",
                       modelType[1], "/taus.RData"))

# Split by convergence #####

#**GLM #####
load(paste0(path2, "/models/",
            models[1], "/results.RData" ))
results_glm <- results
results_glm$model <- rep(models[1], length(results_glm$species))
rm(results)

#**ESM_simple #####
load(paste0(path2, "/models/",
            models[2], "/results.RData" ))
results_ESM_simple <- results
results_ESM_simple$model <- rep(models[2], length(results_ESM_simple$species))
rm(results)
results_ESM_simple_converged <- subset(results_ESM_simple, !is.na(results_ESM_simple$auc))

#**ESM_complex #####
load(paste0(path2, "/models/",
            models[3], "/results.RData" ))
results_ESM_complex <- results
results_ESM_complex$model <- rep(models[3], length(results_ESM_complex$species))
rm(results)
results_ESM_complex_converged <- subset(results_ESM_complex, !is.na(results_ESM_complex$auc))

#**Hmsc_joint #####
load(paste0(path2, "/models/",
            models[4], "/results.RData" ))
results_Hmsc_joint <- results
results_Hmsc_joint$model <- rep(models[4], length(results_Hmsc_joint$species))
rm(results)
results_Hmsc_joint_converged <- subset(results_Hmsc_joint, results_Hmsc_joint$Rhat < 1.2)

#**SAM#####
load(paste0(path2, "/models/",
            models[5], "/results.RData" ))
results_SAM <- results
results_SAM$model <- rep(models[5], length(results_SAM$species))
rm(results)
results_SAM_converged <- subset(results_SAM, !is.na(results_SAM$auc))

# COMBINE INTO ONE #####

full <- rbind(results_glm,
              results_ESM_simple,
              results_ESM_complex,
              results_Hmsc_joint,
              results_SAM)

converged <- rbind(results_glm,
                   results_ESM_simple_converged,
                   results_ESM_complex_converged,
                   results_Hmsc_joint_converged,
                   results_SAM_converged)

#RECODE SIZES #####
full$size <- full$size %>% recode("size2" = "2", "size4" = "4", "size8"="8", "size16"="16", "size32"="32",
                                  "size64"="64") %>%
  factor(levels=c("2", "4", "8", "16", "32","64"))
converged$size <- converged$size %>% recode("size2" = "2", "size4" = "4", "size8"="8", "size16"="16", "size32"="32",
                                            "size64"="64") %>%
  factor(levels=c("2", "4", "8", "16", "32","64"))

full$species <- factor(full$species)
full$model <- factor(full$model, levels=c("glm", "ESM_simple", "ESM_complex", "Hmsc_joint", "SAM"))

converged$species <- factor(converged$species)
converged$model <- factor(converged$model, levels=c("glm", "ESM_simple", "ESM_complex", "Hmsc_joint", "SAM"))

#Calculate plotting positions #####

for (i in 1:length(full$model)) {
  if(full$model[i] == "glm") {
    full$pos[i] = as.numeric(full$size[i]) - 0.16
  }
  if(full$model[i] == "ESM_simple") {
    full$pos[i] = as.numeric(full$size[i]) - 0.08
  }

  if(full$model[i] == "ESM_complex") {
    full$pos[i] = as.numeric(full$size[i])}

  if(full$model[i] == "Hmsc_joint") {
    full$pos[i] = as.numeric(full$size[i]) + 0.08
  }
  if(full$model[i] == "SAM") {
    full$pos[i] = as.numeric(full$size[i]) + 0.16
  }
}

for (i in 1:length(converged$model)) {
  if(converged$model[i] == "glm") {
    converged$pos[i] = as.numeric(converged$size[i]) - 0.16
  }
  if(converged$model[i] == "ESM_simple") {
    converged$pos[i] = as.numeric(converged$size[i]) - 0.08
  }

  if(converged$model[i] == "ESM_complex") {
    converged$pos[i] = as.numeric(converged$size[i])}

  if(converged$model[i] == "Hmsc_joint") {
    converged$pos[i] = as.numeric(converged$size[i]) + 0.08
  }
  if(converged$model[i] == "SAM") {
    converged$pos[i] = as.numeric(converged$size[i]) + 0.16
  }
}

# SEPARATE BY SPECIES #####
dat_broad_avg_full <- subset(full, full$species == "broad_avg")
dat_broad_avg_converged <- subset(converged, converged$species == "broad_avg")
dat_broad_ext_full <- subset(full, full$species == "broad_ext")
dat_broad_ext_converged <- subset(converged, converged$species == "broad_ext")
dat_narrow_avg_full <- subset(full, full$species == "narrow_avg")
dat_narrow_avg_converged <- subset(converged, converged$species == "narrow_avg")
dat_narrow_ext_full <- subset(full, full$species == "narrow_ext")
dat_narrow_ext_converged <- subset(converged, converged$species == "narrow_ext")

write.csv(dat_broad_avg_full, file="../data/dat_broad_avg_full.csv")
write.csv(dat_broad_avg_converged, file="../data/dat_broad_avg_converged.csv")
write.csv(dat_broad_ext_full, file="../data/dat_broad_ext_full.csv")
write.csv(dat_broad_ext_converged, file="../data/dat_broad_ext_converged.csv")
write.csv(dat_narrow_avg_full, file="../data/dat_narrow_avg_full.csv")
write.csv(dat_narrow_avg_converged, file="../data/dat_narrow_avg_converged.csv")
write.csv(dat_narrow_ext_full, file="../data/dat_narrow_ext_full.csv")
write.csv(dat_narrow_ext_converged, file="../data/dat_narrow_ext_converged.csv")

# CALCULATE SAMPLE SIZE #####
sample_size_broad_avg = dat_broad_avg_converged %>% group_by(model, size) %>% summarize(num=n())
sample_size_broad_ext = dat_broad_ext_converged %>% group_by(model, size) %>% summarize(num=n())
sample_size_narrow_avg = dat_narrow_avg_converged %>% group_by(model, size) %>% summarize(num = n())
sample_size_narrow_ext = dat_narrow_ext_converged %>% group_by(model, size) %>% summarize(num = n())

sample_size_broad_avg$pos <- rep(NA, length(sample_size_broad_avg$model))
#Assign pos value
# Broad Avg
for (i in 1:length(sample_size_broad_avg$model)) {
  if(sample_size_broad_avg$model[i] == "glm") {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) - 0.4
  }
  if(sample_size_broad_avg$model[i] == "ESM_simple") {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) - 0.2
  }

  if(sample_size_broad_avg$model[i] == "ESM_complex") {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i])}

  if(sample_size_broad_avg$model[i] == "Hmsc_joint") {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) + 0.2
  }
  if(sample_size_broad_avg$model[i] == "SAM") {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) + 0.4
  }
}
#Broad Ext
for (i in 1:length(sample_size_broad_ext$model)) {
  if(sample_size_broad_ext$model[i] == "glm") {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) - 0.4
  }
  if(sample_size_broad_ext$model[i] == "ESM_simple") {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) - 0.2
  }

  if(sample_size_broad_ext$model[i] == "ESM_complex") {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i])}

  if(sample_size_broad_ext$model[i] == "Hmsc_joint") {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) + 0.2
  }
  if(sample_size_broad_ext$model[i] == "SAM") {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) + 0.4
  }
}
#Narrow Avg
for (i in 1:length(sample_size_narrow_avg$model)) {
  if(sample_size_narrow_avg$model[i] == "glm") {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) - 0.4
  }
  if(sample_size_narrow_avg$model[i] == "ESM_simple") {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) - 0.2
  }

  if(sample_size_narrow_avg$model[i] == "ESM_complex") {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i])}

  if(sample_size_narrow_avg$model[i] == "Hmsc_joint") {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) + 0.2
  }
  if(sample_size_narrow_avg$model[i] == "SAM") {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) + 0.4
  }
}
#Narrow Ext
for (i in 1:length(sample_size_narrow_ext$model)) {
  if(sample_size_narrow_ext$model[i] == "glm") {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) - 0.4
  }
  if(sample_size_narrow_ext$model[i] == "ESM_simple") {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) - 0.2
  }

  if(sample_size_narrow_ext$model[i] == "ESM_complex") {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i])}

  if(sample_size_narrow_ext$model[i] == "Hmsc_joint") {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) + 0.2
  }
  if(sample_size_narrow_ext$model[i] == "SAM") {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) + 0.4
  }
}


write.csv(sample_size_broad_avg, file="../data/sample_size_broad_avg.csv")
write.csv(sample_size_broad_ext, file="../data/sample_size_broad_ext.csv")
write.csv(sample_size_narrow_avg, file="../data/sample_size_narrow_avg.csv")
write.csv(sample_size_narrow_ext, file="../data/sample_size_narrow_ext.csv")