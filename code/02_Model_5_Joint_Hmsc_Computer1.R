# Run Hmsc joint on computer 1
if(Sys.info()['sysname'] == "Darwin") {
path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
path <- "C:/Users/kerickson/Documents/GitHub"
}

source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

repStart <- 76
repEnd <- 100

session <- sessionInfo()
save(session, file=paste0(path2, "/models/Hmsc_joint/sessionInfo_computer1.RData"))



source(paste0(path, "/SDMS_for_rare_species_modeling/code/02_Model_5_Joint_Hmsc.R"))
