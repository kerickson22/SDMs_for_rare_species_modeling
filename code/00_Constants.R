#00_Constants
path <- "/Users/curculion/Documents/GitHub"
# This script contains universal constants that
# are shared across ALL scripts


# This script should be sourced by other scripts

# Load packages #####
library(sp)
library(rstudioapi)
library(ggplot2)

## Set working directory #####
# Set working directory to current file location
setwd(dirname(getActiveDocumentContext()$path))       # Set working directory to source file location

# Define constants #####
sizes <- c(2, 4, 8, 16, 32, 64)
numTestPresences <- 64 #How many test presences to use

## Load objects created from 01_Simulate_Species #####

if(file.exists(paste0(path, "/SDMs_for_rare_species_modeling/data/south.csv" ))){
  south <- read.csv(paste0(path, "/SDMs_for_rare_species_modeling/data/south.csv" ),
                  header=TRUE)}

us <- map_data('state')




