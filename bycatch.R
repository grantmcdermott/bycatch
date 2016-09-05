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

##################################################
##################################################
### Add pctred fvfmey columns to upsides table ###
##################################################
##################################################

ids <- c()
eqfvfmeys <- c()

# calculate equilibrium MEYs
for (i in 1:length(upsides$idorig)) {
  ids <- append(ids, upsides$idorig[i])
  eqfvfmeys <- append(eqfvfmeys, 
                      (1/(optim(par = 0.005, fn = function (x) {
                        - ((upsides$price[i] * upsides$g[i] * x * upsides$k[i] * ((1 - ((upsides$g[i] * x * upsides$phi[i])/(upsides$g[i] * (upsides$phi[i] + 1))))^(1/upsides$phi[i]))) - 
                             (upsides$marginalcost[i] * ((upsides$g[i] * x)^upsides$beta[i])))
                      }, method = "Brent", lower = 0, upper = 1.5)$par * (1/upsides$fvfmsy[i])))
  )
}

# add equilibrium MEY cols to upsides
test <- data_frame(idorig = ids, eqfvfmey = eqfvfmeys)
upsides <- left_join(upsides, test, by = "idorig")

# make pctred columns
upsides <- upsides %>%
  mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
  mutate(pctredfmey = 100 * (1 - (1/eqfvfmey))) # defines pctredmey in terms of eqfmey, but we can change this if we want.
# mutate(pctredfmey = 100 * (1 - (1/(fvfmsy/fmeyvfmsy)))) # defines pctredmey in terms of NPV fmey


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
n1 <- 1000
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
