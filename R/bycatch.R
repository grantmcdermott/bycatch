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
library(pbmcapply)

## Choose the number of process to run in parallel. Will choose the smaller of
## 4 or the number of CPUs your computer has. Edit this as needed. Also be advised
## that it is better to turn off all other applications while this code is
## runnning, to avoid overloading your computer.
num_cores <- min(4, detectCores())

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

test <- 
  pbmclapply(
    1:length(upsides$idorig),
    mc.cores = num_cores,
    function(i){
      
      ids <- upsides$idorig[i]
      
      eqfvfmeys <- 
        1/(optim(par = 0.005,  
                 fn = function(x){
                   - ((upsides$price[i] * upsides$g[i] * x * upsides$k[i] * 
                         ((1 - ((upsides$g[i] * x * upsides$phi[i])/(upsides$g[i] * (upsides$phi[i] + 1))))^(1/upsides$phi[i]))) - 
                        (upsides$marginalcost[i] * ((upsides$g[i] * x)^upsides$beta[i])))
                   }, 
                 method = "Brent", 
                 lower = 0, 
                 upper = 1.5)$par * 
             (1/upsides$fvfmsy[i]))
      
      df = data_frame(idorig = ids, eqfvfmey = eqfvfmeys)
      return(df)
         }) %>%
  bind_rows

# add equilibrium MEY cols to upsides
upsides <- left_join(upsides, test, by = "idorig")

# make pctred columns
upsides <- 
  upsides %>%
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

## Choose the percent increment for calculating (marginal) yield and profit
## losses when there is a shortfall.
incrmt <- 1 ## i.e. one tenth of a percent increment in bycatch reduction

###############
### Turtles ###
###############

turtle_species <- (filter(bycatch_df, grp=="turtle"))$species

## Run the bycatch function over all turtle species 
## Each progress bar represents the time needed for a particular turtle stock. 
## There are four turtle stocks, so there will be four progress bars in succession.
## This next line takes 4 and a bit minutes to run on my laptop (with a 1 percent 
## increment).
turtles <- bind_rows(pblapply(turtle_species, bycatch_func))

## % reduction plot
bycatchdistggplot(turtles) +
  ggsave("TablesFigures/turtles1.png", width = 8, height = 6)
bycatchdistggplot(turtles) +
  facet_wrap(~species, scales = "free") + 
  ggsave("TablesFigures/turtles2.png", width = 8, height = 6)
# ## Note: To add axes to each plot while keeping scales="fixed"
# bycatchdistggplot(turtles) +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
#   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

## Yield loss and cost plot
costggplot(turtles)


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
