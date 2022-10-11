# Run Hmsc joint on computer 2
if(Sys.info()['sysname'] == "Darwin") {
path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
path <- "C:/Users/kerickson/Documents/GitHub"
}

source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

reps <- sample(1:100)

session <- sessionInfo()
save(session, file=paste0(path, "/SDMs_for_rare_species_modeling/data/models/Hmsc_joint/sessionInfo_computer2.RData"))



source(paste0(path, "/SDMS_for_rare_species_modeling/code/02_Model_5_Joint_Hmsc.R"))
if(!file.exists(paste0(path2, "/models/",
                       modelType[1], "/results.RData" ))) {
  source(paste0(path, "/SDMS_for_rare_species_modeling/code/03_Collate_Results.R"))
}
