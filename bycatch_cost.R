rm(list = ls())

##################################
########### Data Setup ###########
##################################
require(readr)
require(dplyr)
require(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(pbapply)

###########################################
############ Read in the data #############
###########################################

# upsides_uncert <- read_csv("Data/upsides_uncert.csv", col_types = cols(regionfao = "c"))
upsides <- read_csv("Data/upsides.csv", col_types = cols(regionfao = "c"))

bycatch_df <- read_csv("Data/bycatch_species.csv")

target_df <- read_csv("Data/target_species.csv")


###################################
########## Load functions #########
###################################

source("bycatch_funcs_cost.R")


###############################
########### Results ###########
###############################

## Desired chance (%) that the bycatch reduction threshold is met
pctchance <- 95
## Sampling parameters
n1 <- 100#1000
n2 <- 100


turtle_species_samp <- c("Loggerhead turtle", "Olive ridley turtle (NEI)")#(filter(bycatch_df, grp=="turtle"))$species
## Run the bycatch function over all turtle species (takes 55s on my laptop)
turtles_samp <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
bycatchdistggplot(turtles_samp) +
  facet_wrap(~species, scales = "free")
costggplot(turtles_samp)

# humpback <- bind_rows(pblapply("Humpback dolphin", bycatch_func))
humpback <- bycatch_func("Humpback dolphin")
humpback
bycatchdistggplot(humpback) 
costggplot(humpback) 

