# Run Hmsc joint on computer 3
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

reps <- c(comp3, comp4[25:1], comp2[25:1],
          comp1[25:1])

session <- sessionInfo()
save(session, file=paste0(path2, "/models/Hmsc_joint/sessionInfo_computer3.RData"))



source(paste0(path, "/SDMS_for_rare_species_modeling/code/02_Model_5_Joint_Hmsc.R"))

if(!file.exists(paste0(path2, "/models/",
                       modelType[1], "/results_start.RData" ))) {
  k <- "computer_3"
  save(k, file=paste0(path2, "/models/",
                       modelType[1], "/results_start.RData"))
  myfile <- paste0(path2, "/models/",
                   modelType[1], "/results_start.RData" ) %>%
    drive_upload(paste0("status_updates_for_Hmsc_joint/", k, ".RData"))
  source(paste0(path, "/SDMS_for_rare_species_modeling/code/03_Collate_Results_joint_Hmsc.R"))
}
