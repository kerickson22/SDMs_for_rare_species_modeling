# Run Hmsc joint on computer 2
if(Sys.info()['sysname'] == "Darwin") {
path <- "/Users/curculion/Documents/GitHub"
path2 <- path
}

if(Sys.info()['sysname'] == "Windows") {
path <- "C:/Users/kerickson/Documents/GitHub"
path2 <- "H:/Global Change Program/Research/ENMs - Modeling Methods for Rare Species (Kelley Erickson)/rare_species/data"
}

source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

reps <- sample(1:100)

session <- sessionInfo()
save(session, file=paste0(path2, "/models/Hmsc_joint/sessionInfo_computer2.RData"))



source(paste0(path, "/SDMS_for_rare_species_modeling/code/02_Model_5_Joint_Hmsc.R"))
