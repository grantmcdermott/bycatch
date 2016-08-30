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

source("bycatch_funcs.R")


###############################
########### Results ###########
###############################

## Desired chance (%) that the bycatch reduction threshold is met
pctchance <- 95
## Sampling parameters
n1 <- 10000
n2 <- 100


###############
### Turtles ###
###############

turtle_species <- (filter(bycatch_df, grp=="turtle"))$species

## Run the bycatch function over all turtle species (takes 55s on my laptop)
turtles <- bind_rows(pblapply(turtle_species, bycatch_func))

## Plot the data
bycatchdistggplot(turtles) +
  ggsave("TablesFigures/turtles1.png", width = 8, height = 6)
bycatchdistggplot(turtles) +
  facet_wrap(~species, scales = "free") + 
  ggsave("TablesFigures/turtles2.png", width = 8, height = 6)
# ## Note: To add axes to each plot while keeping scales="fixed"
# bycatchdistggplot(turtles) +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
#   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)


###############
### Mammals ###
###############

mammal_species <- (filter(bycatch_df, grp=="mammal"))$species

## Run the bycatch function over all turtle species (takes 47s on my laptop)
mammals <- bind_rows(pblapply(mammal_species, bycatch_func))

## Plot the data
bycatchdistggplot(mammals) 
## NZ sea lion is a huge outlier
bycatchdistggplot(filter(mammals, species != "NZ sea lion")) + 
  ggsave("TablesFigures/mammals1.png", width = 8, height = 6)
bycatchdistggplot(mammals) + 
  facet_wrap(~species, scales = "free") + 
  ggsave("TablesFigures/mammals2.png", width = 8, height = 6)


#############
### Birds ###
#############

bird_species <- (filter(bycatch_df, grp=="bird"))$species

## Run the bycatch function over all turtle species (takes 01m 14s on my laptop)
birds <- bind_rows(pblapply(bird_species, bycatch_func))

## Plot the data
bycatchdistggplot(birds) 
## Sooty shearwater is a huge outlier
bycatchdistggplot(filter(birds, species != "Sooty shearwater")) + 
  facet_wrap(~species, ncol = 2) + 
  ggsave("TablesFigures/birds1.png", width = 8, height = 6)
bycatchdistggplot(filter(birds, species != "Sooty shearwater")) + 
  facet_wrap(~species, scales = "free", ncol = 2) + 
  ggsave("TablesFigures/birds2.png", width = 8, height = 6)
