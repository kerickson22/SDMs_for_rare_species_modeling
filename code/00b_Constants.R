#00b_Constants
if(Sys.info()['sysname'] == "Darwin") {
  path <- "/Users/curculion/Documents/GitHub"
  path2 <- "/Users/curculion/Documents/GitHub/SDMs_for_rare_species_modeling/data"
}

if(Sys.info()['sysname'] == "Windows") {
  path <- "C:/Users/kerickson/Documents/GitHub"
  path2 <- "H:/Global Change Program/Research/ENMs - Modeling Methods for Rare Species (Kelley Erickson)/rare_species/data"
}

# This script contains universal constants that
# are shared across ALL scripts


# This script should be sourced by other scripts

# Load packages #####
library(sp)
#library(rstudioapi)
library(ggplot2)
library(RColorBrewer)
library(mixtools) #for drawing ellipses
library(Hmsc)
library(ecomix)
library(MCMCvis)
library(enmSdm)
library(arm) #for invlogit
library(PRROC) #for AUCpr
library(googledrive)
library(dplyr)
#library(ggpubr)
library(abind)

## Set working directory #####
# Set working directory to current file location
setwd(paste0(path, "/SDMS_for_rare_species_modeling/code"))       # Set working directory to source file location

# Define constants #####
numPresences <- c(2, 4, 8, 16, 32, 64)
sizes <- paste0("size", numPresences)# Desired sample size
numTestPresences <- 64 #How many test presences to use
species <- c("broad_avg", "broad_ext", "narrow_avg", "narrow_ext")
models <- c("glm", "ESM_simple", "ESM_complex", "Hmsc_joint","SAM")
numReplicates <- 100
replicates <- paste0("rep", 1:numReplicates)

x1 <- seq(from=-9, to = 13, length.out=100)
x2 <- seq(from=-9, to = 16, length.out=100)
x3 <- seq(from=-11, to = 4, length.out=100)

#Load real species data (object created in 00a_Run_Once.R)
if(file.exists(paste0(path, "/SDMs_for_rare_species_modeling/data/south.csv" ))){
  south <- read.csv(paste0(path, "/SDMs_for_rare_species_modeling/data/south.csv" ),
                    header=TRUE)}
if(file.exists(paste0(path, "/SDMs_for_rare_species_modeling/data/X_bar.csv" ))){
  X_bar <- read.csv(file=paste0(path, "/SDMs_for_rare_species_modeling/data/X_bar.csv" ))}

NSites <- length(south$long)
row.names(south) <- 1:NSites



## Load objects created from 01_Simulate_Species #####

if(file.exists(paste0(path, "/SDMs_for_rare_species_modeling/data/south.csv" ))){
  load(paste0(path, "/SDMs_for_rare_species_modeling/data/sims_data.RData" ))}

#Plotting Information #####
#us <- map_data('state')

nicheCols <- brewer.pal(6, "Dark2")

#Functions
computeRMSEWeighted = function(Y, predY) {
  RMSE <- sqrt((mean((Y[Y==1]-predY[Y==1])^2, na.rm=TRUE) +
                  mean((Y[Y==0]-predY[Y==0])^2, na.rm=TRUE)))
  return(RMSE)
}

computeTjurR2 = function(Y, predY) {


    R2 = mean(predY[which(Y == 1)]) - mean(predY[which(Y == 0)])

  return(R2)
}

# This is not the correct function to use
# computeR2 = function(Y, predY, method="pearson"){
#
#
#   co = cor(Y, predY, method=method, use='pairwise')
#   R2 = sign(co)*co^2
#
#   return(R2)
# }

#Create directories for storing models

for(i in 1:length(models)){
  if(!file.exists(paste0(path2, "/models/",
                         models[i]))){
    dir.create(paste0(path2, "/models/",
                      models[i]))
  }
  if(!file.exists(paste0(path2, "/models/",
                         models[i], "/status"))){
    dir.create(paste0(path2,"/models/",
                      models[i], "/status"))
  }
  for(j in 1:length(species)){
    if(!file.exists(paste0(path2,"/models/",
                           models[i], "/",
                           species[j]))){
      dir.create(paste0(path2,"/models/",
                        models[i], "/",
                        species[j]))
    }
    for(k in 1:length(sizes)){
      if(!file.exists(paste0(path2, "/models/",
                             models[i], "/",
                             species[j], "/",
                             sizes[k]))){
        dir.create(path=paste0(path2, "/models/",
                               models[i], "/",
                               species[j], "/",
                               sizes[k]))
      }
    }
  }
}


#These are defined in the script 00a_Run_Once.R using
# data from the real species community
var_narrow <- diag(c(0.9, 2.5))
sd_narrow <- sqrt(var_narrow)

mu_ext <- c(-3.5, 2)
mu_avg <- mu_avg <- c(mean(X_bar$mu_PC1), mean(X_bar$mu_PC3))
var_broad <- diag(x=c(quantile(X_bar$sd_PC1, probs=0.90, na.rm=T),
                     quantile(X_bar$sd_PC3, probs=0.90, na.rm=T)))^2
sd_broad <- sqrt(var_broad)


#ESM formulas #####
formulaMatrix <- read.csv("../data/formulaMatrix.csv",
                          header=T)

#Realized Niche Estimates #####
#These are created in #00a_Run_Once
load("../data/estSp.RData")

#Color Palette for plotting
#https://personal.sron.nl/~pault/#sec:qualitative
bright <- c("#4477AA", "#66CCEE",
            "#228833", "#CCBB44",
            "#EE6677", "#AA3377",
            "#BBBBBB")
my_bright <- bright[c(3,1, 2, 5, 6)]

#blue, cyan, green, yellow, red, purple, grey


