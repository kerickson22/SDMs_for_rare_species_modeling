# Run Hmsc joint on computer 1

#Computer 1 doesn't seem to work well with google

if(Sys.info()['sysname'] == "Darwin") {
path <- "/Users/curculion/Documents/GitHub"
}

if(Sys.info()['sysname'] == "Windows") {
path <- "C:/Users/kerickson/Documents/GitHub"
}

source(paste0(path, "/SDMs_for_rare_species_modeling/code/00b_Constants.R"))

comp1 <- 1:25
comp2 <- 26:50
comp3 <- 51:75
comp4 <- 76:100

reps <- c(comp1, comp2[25:1], comp3[25:1],
          comp4[25:1])

session <- sessionInfo()
save(session, file=paste0(path2, "/models/Hmsc_joint/sessionInfo_computer1.RData"))



source(paste0(path, "/SDMS_for_rare_species_modeling/code/02_Model_5_Joint_Hmsc_no_google.R"))
