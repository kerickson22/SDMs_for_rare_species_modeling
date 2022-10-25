#01_Simulate_Species

# Script for simulating virtual species.

# Outputs of this script:
# *sims_data.RData

if(Sys.info()['sysname'] == "Darwin") {
path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
path <- "C:/Users/kerickson/Documents/GitHub"
}

source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))


## Virtual Species Niches #####

set.seed(12)
#?To set seed or not?

sims <- vector("list", length(species))
names(sims) <- species
for(s in 1:length(species)){
  sims[[s]] <- vector("list", length(numPresences))
  names(sims[[s]]) <- sizes
  for(n in 1:length(numPresences)) {
    sims[[s]][[n]] <- vector("list", length(replicates))
    names(sims[[s]][[n]]) <- replicates

  }}


for (s in 1:length(species)) {
  switch(species[s],
         "broad_avg"={
           sims[[s]]$L <-dnorm(south[,3], mu_avg[1],
                               sd=sd_broad[1,1]) * dnorm(south[,5],
                              mu_avg[2], sd=sd_broad[2,2])
          sims[[s]]$mu <- mu_avg
          sims[[s]]$sd <- sd_broad
         },
         "broad_ext" = {
           sims[[s]]$L <-dnorm(south[,3], mu_ext[1],
                    sd=sd_broad[1,1]) * dnorm(south[,5],
                    mu_ext[2], sd=sd_broad[2,2])
           sims[[s]]$mu <- mu_ext
           sims[[s]]$sd <- sd_broad
         },
         "narrow_avg" = {
           sims[[s]]$L <-dnorm(south[,3],
                mu_avg[1], sd=sd_narrow[1,1]) *
             dnorm(south[,5], mu_avg[2], sd=sd_narrow[2,2])
           sims[[s]]$mu <- mu_avg
           sims[[s]]$sd <- sd_narrow
           },
         "narrow_ext" = {
           sims[[s]]$L <-dnorm(south[,3],
                mu_ext[1],sd=sd_narrow[1,1]) *
             dnorm(south[,5], mu_ext[2], sd=sd_narrow[2,2])
           sims[[s]]$mu <- mu_ext
           sims[[s]]$sd <- sd_narrow
           })#end of switch

  for(n in 1:length(sizes)) {

    for(r in 1:length(replicates)){
      sims[[s]][[n]][[r]]$index <- matrix(ncol=(numPresences[n]+64))
      sims[[s]][[n]][[r]]$YSim <- matrix(0, nrow=NSites, ncol=1)
      sims[[s]][[n]][[r]]$YSim  <- data.frame(sims[[s]][[n]][[r]]$YSim )
      names(sims[[s]][[n]][[r]]$YSim ) <- "SimSp"

      sims[[s]][[n]][[r]]$index <- sample(1:NSites, (numPresences[n]+64), replace=F, prob=sims[[s]]$L/sum(sims[[s]]$L))

      for(site in 1:NSites) {
        if(site %in% sims[[s]][[n]][[r]]$index) {
          sims[[s]][[n]][[r]]$YSim[site,] <- 1
        }
      }

      sites <- 1:NSites
      sims[[s]][[n]][[r]]$presences_train <- sims[[s]][[n]][[r]]$index[1:numPresences[n]]
      sims[[s]][[n]][[r]]$presences_test <- sims[[s]][[n]][[r]]$index[(numPresences[n]+1):(numPresences[n]+64)]
      NAbsences_test <- 64
      NAbsences_train <- NSites - numPresences[n]-NAbsences_test-64
      absences <- sites[! sites %in% sims[[s]][[n]][[r]]$index]
      sims[[s]][[n]][[r]]$absences_train <- sample(absences, NAbsences_train, replace=F)
      sims[[s]][[n]][[r]]$absences_test <- absences[!absences %in% sims[[s]][[n]][[r]]$absences_train]
      sims[[s]][[n]][[r]]$train <- c(sims[[s]][[n]][[r]]$presences_train, sims[[s]][[n]][[r]]$absences_train)
      sims[[s]][[n]][[r]]$test <- c(sims[[s]][[n]][[r]]$presences_test, sims[[s]][[n]][[r]]$absences_test)


    }
  }
}


save(sims, file="../data/sims_data.RData")