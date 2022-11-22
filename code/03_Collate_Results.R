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
#As of 11/10 83 reps have run
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
converged$species <- factor(converged$species)


full$model <- full$model %>%
  recode("glm" = "GLM", "ESM_simple" = "ESM-L",
         "ESM_complex"="ESM-P",
         "Hmsc_joint"="HMSC", "SAM"="SAM") %>%
  factor(levels=c("GLM", "ESM-L", "ESM-P", "HMSC", "SAM"))

converged$model <- converged$model %>%
  recode("glm" = "GLM", "ESM_simple" = "ESM-L",
         "ESM_complex"="ESM-P",
         "Hmsc_joint"="HMSC", "SAM"="SAM") %>%
  factor(levels=c("GLM", "ESM-L", "ESM-P", "HMSC", "SAM"))



#Calculate plotting positions #####

for (i in 1:length(full$model)) {
  if(full$model[i] == levels(full$model)[1]) {
    full$pos[i] = as.numeric(full$size[i]) - 0.16
  }
  if(full$model[i] == levels(full$model)[2]) {
    full$pos[i] = as.numeric(full$size[i]) - 0.08
  }

  if(full$model[i] == levels(full$model)[3]) {
    full$pos[i] = as.numeric(full$size[i])}

  if(full$model[i] == levels(full$model)[4]) {
    full$pos[i] = as.numeric(full$size[i]) + 0.08
  }
  if(full$model[i] == levels(full$model)[5]) {
    full$pos[i] = as.numeric(full$size[i]) + 0.16
  }
}

for (i in 1:length(converged$model)) {
  if(converged$model[i] == levels(converged$model)[1]) {
    converged$pos[i] = as.numeric(converged$size[i]) - 0.16
  }
  if(converged$model[i] == levels(converged$model)[2]) {
    converged$pos[i] = as.numeric(converged$size[i]) - 0.08
  }

  if(converged$model[i] == levels(converged$model)[3]) {
    converged$pos[i] = as.numeric(converged$size[i])}

  if(converged$model[i] == levels(converged$model)[4]) {
    converged$pos[i] = as.numeric(converged$size[i]) + 0.08
  }
  if(converged$model[i] == levels(converged$model)[5]) {
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
  if(sample_size_broad_avg$model[i] == levels(full$model)[1]) {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) - 0.4
  }
  if(sample_size_broad_avg$model[i] == levels(full$model)[2]) {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) - 0.2
  }

  if(sample_size_broad_avg$model[i] == levels(full$model)[3]) {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i])}

  if(sample_size_broad_avg$model[i] == levels(full$model)[4]) {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) + 0.2
  }
  if(sample_size_broad_avg$model[i] == levels(full$model)[5]) {
    sample_size_broad_avg$pos[i] = as.numeric(sample_size_broad_avg$size[i]) + 0.4
  }
}

#Broad Ext
for (i in 1:length(sample_size_broad_ext$model)) {
  if(sample_size_broad_ext$model[i] == levels(full$model)[1]) {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) - 0.4
  }
  if(sample_size_broad_ext$model[i] == levels(full$model)[2]) {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) - 0.2
  }

  if(sample_size_broad_ext$model[i] == levels(full$model)[3]) {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i])}

  if(sample_size_broad_ext$model[i] == levels(full$model)[4]) {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) + 0.2
  }
  if(sample_size_broad_ext$model[i] == levels(full$model)[5]) {
    sample_size_broad_ext$pos[i] = as.numeric(sample_size_broad_ext$size[i]) + 0.4
  }
}
#Narrow Avg
for (i in 1:length(sample_size_narrow_avg$model)) {
  if(sample_size_narrow_avg$model[i] == levels(full$model)[1]) {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) - 0.4
  }
  if(sample_size_narrow_avg$model[i] == levels(full$model)[2]) {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) - 0.2
  }

  if(sample_size_narrow_avg$model[i] == levels(full$model)[3]) {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i])}

  if(sample_size_narrow_avg$model[i] == levels(full$model)[4]) {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) + 0.2
  }
  if(sample_size_narrow_avg$model[i] == levels(full$model)[5]) {
    sample_size_narrow_avg$pos[i] = as.numeric(sample_size_narrow_avg$size[i]) + 0.4
  }
}
#Narrow Ext
for (i in 1:length(sample_size_narrow_ext$model)) {
  if(sample_size_narrow_ext$model[i] == levels(full$model)[1]) {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) - 0.4
  }
  if(sample_size_narrow_ext$model[i] == levels(full$model)[2]) {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) - 0.2
  }

  if(sample_size_narrow_ext$model[i] == levels(full$model)[3]) {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i])}

  if(sample_size_narrow_ext$model[i] == levels(full$model)[4]) {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) + 0.2
  }
  if(sample_size_narrow_ext$model[i] == levels(full$model)[5]) {
    sample_size_narrow_ext$pos[i] = as.numeric(sample_size_narrow_ext$size[i]) + 0.4
  }
}


write.csv(sample_size_broad_avg, file="../data/sample_size_broad_avg.csv")
write.csv(sample_size_broad_ext, file="../data/sample_size_broad_ext.csv")
write.csv(sample_size_narrow_avg, file="../data/sample_size_narrow_avg.csv")
write.csv(sample_size_narrow_ext, file="../data/sample_size_narrow_ext.csv")

#Calculate correlation counts #####
cor <- expand.grid(models, c(2, 4, 8, 16, 32, 64))
names(cor) <- c("model", "size")
cor$size <- factor(cor$size)
cor$model <- cor$model %>%
  recode("glm" = "GLM", "ESM_simple" = "ESM-L",
         "ESM_complex"="ESM-P",
         "Hmsc_joint"="HMSC", "SAM"="SAM") %>%
  factor(levels=c("GLM", "ESM-L", "ESM-P", "HMSC", "SAM"))



cor$PC1_neg <- rep(NA, length(cor$size))
cor$PC1_neutral <- rep(NA, length(cor$size))
cor$PC1_pos <- rep(NA, length(cor$size))
cor$PC3_neg <- rep(NA, length(cor$size))
cor$PC3_neutral <- rep(NA, length(cor$size))
cor$PC3_pos <- rep(NA, length(cor$size))



cor_broad_avg <- cor
cor_broad_ext <- cor
cor_narrow_avg <- cor
cor_narrow_ext <- cor

tabulate_cor <- function(dat, cor) {
  for (i in 1:length(cor$size)) {
    temp <- subset(dat,
                   dat$model == cor$model[i] &
                     dat$size == cor$size[i])
    cor$PC1_neg[i] <- nrow(temp[temp$PC1_rankCor < - 0.33, ])
    cor$PC1_neutral[i] <-nrow(temp[temp$PC1_rankCor >= - 0.33 & temp$PC1_rankCor < 0.33, ])
    cor$PC1_pos[i] <- nrow(temp[temp$PC1_rankCor >= 0.33,])

    cor$PC3_neg[i] <- nrow(temp[temp$PC3_rankCor < - 0.33, ])
    cor$PC3_neutral[i] <-nrow(temp[temp$PC3_rankCor >= - 0.33 & temp$PC3_rankCor < 0.33, ])
    cor$PC3_pos[i] <- nrow(temp[temp$PC3_rankCor >=0.33,])

  }
  return(cor)
}

cor_broad_avg <- tabulate_cor(dat_broad_avg_converged, cor_broad_avg)
cor_broad_ext <- tabulate_cor(dat_broad_ext_converged, cor_broad_ext)
cor_narrow_avg <- tabulate_cor(dat_narrow_avg_converged, cor_narrow_avg)
cor_narrow_ext <- tabulate_cor(dat_narrow_ext_converged, cor_narrow_ext)

cor_broad_avg_long_PC1 <- gather(cor_broad_avg[,c(1,2, 3:5)], "PC1", "PC1_counts", PC1_neg, PC1_neutral, PC1_pos)
cor_broad_avg_long_PC3 <- gather(cor_broad_avg[,c(1,2, 6:8)], "PC3", "PC3_counts", PC3_neg, PC3_neutral, PC3_pos)
cor_broad_ext_long_PC1 <- gather(cor_broad_ext[,c(1,2, 3:5)], "PC1", "PC1_counts", PC1_neg, PC1_neutral, PC1_pos)
cor_broad_ext_long_PC3 <- gather(cor_broad_ext[,c(1,2, 6:8)], "PC3", "PC3_counts", PC3_neg, PC3_neutral, PC3_pos)
cor_narrow_avg_long_PC1 <- gather(cor_narrow_avg[,c(1,2, 3:5)], "PC1", "PC1_counts", PC1_neg, PC1_neutral, PC1_pos)
cor_narrow_avg_long_PC3 <- gather(cor_narrow_avg[,c(1,2, 6:8)], "PC3", "PC3_counts", PC3_neg, PC3_neutral, PC3_pos)
cor_narrow_ext_long_PC1 <- gather(cor_narrow_ext[,c(1,2, 3:5)], "PC1", "PC1_counts", PC1_neg, PC1_neutral, PC1_pos)
cor_narrow_ext_long_PC3 <- gather(cor_narrow_ext[,c(1,2, 6:8)], "PC3", "PC3_counts", PC3_neg, PC3_neutral, PC3_pos)

cor_broad_avg_long <- cbind(cor_broad_avg_long_PC1, cor_broad_avg_long_PC3[,3:4])
cor_broad_ext_long <- cbind(cor_broad_ext_long_PC1, cor_broad_ext_long_PC3[,3:4])
cor_narrow_avg_long <- cbind(cor_narrow_avg_long_PC1, cor_narrow_avg_long_PC3[,3:4])
cor_narrow_ext_long <- cbind(cor_narrow_ext_long_PC1, cor_narrow_ext_long_PC3[,3:4])

cor_broad_avg_long$PC1 <- factor(cor_broad_avg_long$PC1, levels=c("PC1_neg", "PC1_neutral", "PC1_pos"))
cor_broad_avg_long$PC3 <- factor(cor_broad_avg_long$PC3, levels=c("PC3_neg", "PC3_neutral", "PC3_pos"))
cor_broad_ext_long$PC1 <- factor(cor_broad_ext_long$PC1, levels=c("PC1_neg", "PC1_neutral", "PC1_pos"))
cor_broad_ext_long$PC3 <- factor(cor_broad_ext_long$PC3, levels=c("PC3_neg", "PC3_neutral",  "PC3_pos"))

cor_narrow_avg_long$PC1 <- factor(cor_narrow_avg_long$PC1, levels=c("PC1_neg", "PC1_neutral", "PC1_pos"))
cor_narrow_avg_long$PC3 <- factor(cor_narrow_avg_long$PC3, levels=c("PC3_neg", "PC3_neutral", "PC3_pos"))
cor_narrow_ext_long$PC1 <- factor(cor_narrow_ext_long$PC1, levels=c("PC1_neg", "PC1_neutral", "PC1_pos"))
cor_narrow_ext_long$PC3 <- factor(cor_narrow_ext_long$PC3, levels=c("PC3_neg", "PC3_neutral", "PC3_pos"))

write.csv(cor_broad_avg_long, file="../data/cor_broad_avg_long.csv")
write.csv(cor_broad_ext_long, file="../data/cor_broad_ext_long.csv")
write.csv(cor_narrow_avg_long, file="../data/cor_narrow_avg_long.csv")
write.csv(cor_narrow_ext_long, file="../data/cor_narrow_ext_long.csv")

#Correlation Percentages #####
cor_broad_avg_percent <- cor_broad_avg_long
cor_broad_ext_percent <- cor_broad_ext_long
cor_narrow_avg_percent <- cor_narrow_avg_long
cor_narrow_ext_percent <- cor_narrow_ext_long

#Denominator

getPercent <- function(cor, sample) {
  cor$denom <- rep(NA, length(cor$model))
  cor$PC1Frac <- rep(NA, length(cor$model))
  cor$PC3Frac <- rep(NA, length(cor$model))

  for(i in 1:length(cor$model)) {
    temp <- subset(sample, sample$size == cor$size[i])
    temp <- subset(temp, temp$model == cor$model[i])
    if(length(temp$num) == 0) {
      cor$denom[i] <- 0
      cor$PC1Frac[i] <- 0
      cor$PC3Frac[i] <- 0} else{
    cor$denom[i] <- temp$num
    cor$PC1Frac[i] <- cor$PC1_counts[i]/cor$denom[i]
    cor$PC3Frac[i] <- cor$PC3_counts[i]/cor$denom[i]
    }}
  return(cor)
}

cor_broad_avg_percent <- getPercent(cor_broad_avg_percent, sample_size_broad_avg)
cor_broad_ext_percent <- getPercent(cor_broad_ext_percent, sample_size_broad_ext)
cor_narrow_avg_percent <- getPercent(cor_narrow_avg_percent, sample_size_narrow_avg)
cor_narrow_ext_percent <- getPercent(cor_narrow_ext_percent, sample_size_narrow_ext)


# cor_broad_avg_percent$model <- recode_factor(cor_broad_avg_percent$model, glm = "glm",
#                                                        ESM_simple = "ESM \n Linear",
#                                                           ESM_complex = "ESM \n Quadratic",
#                                                           Hmsc_joint = "Hmsc",
#                                                           SAM = "SAM")
# cor_broad_ext_percent$model <- recode_factor(cor_broad_ext_percent$model, glm = "glm",
#                                              ESM_simple = "ESM \n Linear",
#                                              ESM_complex = "ESM \n Quadratic",
#                                              Hmsc_joint = "Hmsc",
#                                              SAM = "SAM")
# cor_narrow_avg_percent$model <- recode_factor(cor_narrow_avg_percent$model, glm = "glm",
#                                              ESM_simple = "ESM \n Linear",
#                                              ESM_complex = "ESM \n Quadratic",
#                                              Hmsc_joint = "Hmsc",
#                                              SAM = "SAM")
# cor_narrow_ext_percent$model <- recode_factor(cor_narrow_ext_percent$model, glm = "glm",
#                                              ESM_simple = "ESM \n Linear",
#                                              ESM_complex = "ESM \n Quadratic",
#                                              Hmsc_joint = "Hmsc",
#                                              SAM = "SAM")
cor_broad_avg_percent$PC1 <- recode_factor(cor_broad_avg_percent$PC1,
                                           PC1_neg = "Negative",
                                           PC1_neutral = "Neutral",
                                           PC1_pos = "Positive")
cor_broad_avg_percent$PC3 <- recode_factor(cor_broad_avg_percent$PC3,
                                           PC3_neg = "Negative",
                                           PC3_neutral = "Neutral",
                                           PC3_pos = "Positive")
cor_broad_ext_percent$PC1 <- recode_factor(cor_broad_ext_percent$PC1,
                                           PC1_neg = "Negative",
                                           PC1_neutral = "Neutral",
                                           PC1_pos = "Positive")
cor_broad_ext_percent$PC3 <- recode_factor(cor_broad_ext_percent$PC3,
                                           PC3_neg = "Negative",
                                           PC3_neutral = "Neutral",
                                           PC3_pos = "Positive")
cor_narrow_avg_percent$PC1 <- recode_factor(cor_narrow_avg_percent$PC1,
                                           PC1_neg = "Negative",
                                           PC1_neutral = "Neutral",
                                           PC1_pos = "Positive")
cor_narrow_avg_percent$PC3 <- recode_factor(cor_narrow_avg_percent$PC3,
                                           PC3_neg = "Negative",
                                           PC3_neutral = "Neutral",
                                           PC3_pos = "Positive")
cor_narrow_ext_percent$PC1 <- recode_factor(cor_narrow_ext_percent$PC1,
                                           PC1_neg = "Negative",
                                           PC1_neutral = "Neutral",
                                           PC1_pos = "Positive")
cor_narrow_ext_percent$PC3 <- recode_factor(cor_narrow_ext_percent$PC3,
                                           PC3_neg = "Negative",
                                           PC3_neutral = "Neutral",
                                           PC3_pos = "Positive")
write.csv(cor_broad_avg_percent, file="../data/cor_broad_avg_percent.csv")
write.csv(cor_broad_ext_percent, file="../data/cor_broad_ext_percent.csv")
write.csv(cor_narrow_avg_percent, file="../data/cor_narrow_avg_percent.csv")
write.csv(cor_narrow_ext_percent, file="../data/cor_narrow_ext_percent.csv")

# Count up simulated species
countSimSp <- function(species, size) {
  count <- rep(0, length(south$long))
  for (r in 1:100){
    count <- count + sims[[species]][[size]][[replicates[r]]]$YSim
  }
  return(count)
}

counts <- sims
for (s in 1:length(species)) {
  for (n in 1:length(sizes)) {
    count <- rep(0, length(south$long))
    count <- countSimSp(species[s], sizes[n])
    counts[[s]][[n]]$dat <- cbind(south, count)


  }
}

makeSimSpPlot <- function(species, size){
  dat2 <- data.frame(counts[[species]][[size]]$dat)
  dat2 <- subset(dat2, dat2$SimSp != 0)
  p1  <-
    ggplot() +
    geom_point(data=dat, aes(x=PC1, y=PC3, fill="SimSp"))
  }


p <- makeSimSpPlot("broad_avg","size2")

Y <- sims$broad_avg$size2$rep1$YSim
dat <- data.frame(cbind(south, Y))
dat2 <- subset(dat, dat$SimSp !=0)
p1 <- ggplot() +
  geom_point(data=dat2, aes(x=PC1, y=PC3,
                           color=SimSp)) +
  scale_color_viridis()