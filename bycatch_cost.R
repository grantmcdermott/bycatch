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

