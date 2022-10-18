#01_Simulate_Species

# Script for simulating virtual species.

# Outputs of this script:
# *south.csv

if(Sys.info()['sysname'] == "Darwin") {
path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
path <- "C:/Users/kerickson/Documents/GitHub"
}

source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

sims_new <- sims
#The originally generated species are stored in

load("../data/sims_data_original.RData")

sims_original <-sims

#Now go run 01_Simulate_Species
# Once that is done, come back and run
# the below code to fill in sizes 4:64 with
# the original values of species

      sims_new$broad_avg$size4  <- sims_original$broad_avg$size4
      sims_new$broad_avg$size8  <- sims_original$broad_avg$size8
      sims_new$broad_avg$size16 <- sims_original$broad_avg$size16
      sims_new$broad_avg$size32 <- sims_original$broad_avg$size32
      sims_new$broad_avg$size64 <- sims_original$broad_avg$size64

      sims_new$broad_ext$size4  <- sims_original$broad_ext$size4
      sims_new$broad_ext$size8  <- sims_original$broad_ext$size8
      sims_new$broad_ext$size16 <- sims_original$broad_ext$size16
      sims_new$broad_ext$size32 <- sims_original$broad_ext$size32
      sims_new$broad_ext$size64 <- sims_original$broad_ext$size64

      sims_new$narrow_avg$size4  <- sims_original$narrow_avg$size4
      sims_new$narrow_avg$size8  <- sims_original$narrow_avg$size8
      sims_new$narrow_avg$size16 <- sims_original$narrow_avg$size16
      sims_new$narrow_avg$size32 <- sims_original$narrow_avg$size32
      sims_new$narrow_avg$size64 <- sims_original$narrow_avg$size64

      sims_new$narrow_ext$size4  <- sims_original$narrow_ext$size4
      sims_new$narrow_ext$size8  <- sims_original$narrow_ext$size8
      sims_new$narrow_ext$size16 <- sims_original$narrow_ext$size16
      sims_new$narrow_ext$size32 <- sims_original$narrow_ext$size32
      sims_new$narrow_ext$size64 <- sims_original$narrow_ext$size64

sims<-sims_new
save(sims, file="../data/sims_data.RData")
