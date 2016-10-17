rm(list = ls()) #Clear environment

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
upsides_kobe <- read_csv("Kobe MEY data for chris.csv") %>%
  rename(idoriglumped = IdOrig,
         eqfvfmey = current_f_mey,
         eqfmeyvfmsy = f_mey) %>%
  select(idoriglumped, eqfvfmey, eqfmeyvfmsy)

upsides <- read_csv("Data/upsides.csv", col_types = cols(regionfao = "c"))
upsides <- left_join(upsides, upsides_kobe, by = 'idoriglumped') %>%
  mutate(curr_f = g * fvfmsy,
         f_mey = g * eqfmeyvfmsy,
         pctredfmsy = 100 * (1 - (1/fvfmsy)),
         pctredfmey = 100 * (1 - (1/eqfvfmey)))

bycatch_df <- read_csv("Data/bycatch_species.csv")

target_df <- read_csv("Data/target_species.csv")


###################################
########## Load functions #########
###################################

source("bycatch_funcs_cost.R")


###############################
########### Results ###########
###############################

## Sampling parameters
n1 <- 100
n2 <- 100

## Turtle results
turtle_species_samp <- (filter(bycatch_df, grp=="turtle"))$species #c("Loggerhead turtle", "Olive ridley turtle (NEI)")
## Run the bycatch function over all turtle species
turtles_samp <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtledistplots <- bycatchdistggplot(turtles_samp) +
  facet_wrap(~species, scales = "free")
turtlecostplots <- costggplot(turtles_samp)
#write_csv(turtles_samp, "turtles_test.csv")
turtledistplots
turtlecostplots

## Mammal results
mammal_species_samp <- (filter(bycatch_df, grp=="mammal"))$species 
## Run the bycatch function over all mammal species
mammals_samp <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammaldistplots <- bycatchdistggplot(mammals_samp) +
  facet_wrap(~species, scales = "free")
mammalcostplots <- costggplot(mammals_samp)
#write_csv(mammals_samp, "mammals_test.csv")
mammaldistplots
mammalcostplots

## Bird results
bird_species_samp <- (filter(bycatch_df, grp=="bird"))$species 
## Run the bycatch function over all bird species
bird_samp <- bind_rows(pblapply(bird_species_samp, bycatch_func))
birddistplots <- bycatchdistggplot(bird_samp) +
  facet_wrap(~species, scales = "free")
birdcostplots <- costggplot(bird_samp)
#write_csv(bird_samp, "bird_test.csv")
birddistplots
birdcostplots

# humpback <- bind_rows(pblapply("Humpback dolphin", bycatch_func))
humpback <- bycatch_func("Humpback dolphin")
humpback
bycatchdistggplot(humpback) 
costggplot(humpback) 

# humpback <- bind_rows(pblapply("Humpback dolphin", bycatch_func))
loghead <- bycatch_func("Loggerhead")
loghead
bycatchdistggplot(loghead) 
costggplot(loghead) 

###################################################
########### Fig. 2 - Loggerhead Example ###########
###################################################

### Extract distribution plots for specific sub-groups of target species

## Setup
# Function needed below to handle stocks with multiple fao regions
faofun <- function(x) paste("grepl(", toString(x), ",","dt$regionfao)", sep = "")

# Selecting stocks: option 1: list of species categories, fao regions specified
stockselectfig2 <- function(dt,spcat,faoreg){
  # Construct command to filter FAO regions
  a <- lapply(faoreg,faofun)
  b <- paste(a, collapse = "|")
  dtuse <- 
    dt %>%
    filter((eval(parse(text = b)))) %>% 
    filter(speciescat %in% spcat) %>%
    filter(fmeyvfmsy > 0) %>%
    mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
    select(idorig,pctredfmsy,pctredfmey,wt)
    #select(pctredfmsy,pctredfmey)
  return(dtuse)
}

# Plotting distribution
fig2distplot <- function(bdist) {
  bd <- bdist %>%
    rename(MSY = pctredfmsy,
           MEY = pctredfmey)
  # den <- apply(bd, 2, density)
  plot(bd$MSY, bd$wt, xlim=c(-50,100), col=4, pch = 19,
       xlab = "% Reduction in Mortality", ylab = "2012 Effort ($)")
  par(new = TRUE)
  plot(bd$MEY, bd$wt, xlim=c(0,100), col=2, pch = 19,
       xlab = "% Reduction in Mortality", ylab = "2012 Effort ($)")
  #legend("topleft", legend=names(den), fill=c(4,6))
}

## shrimp
shrimpfaoreg <- c(21,31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
shrimpspcat <- c(45)
dtshrimp <- stockselectfig2(upsides, shrimpspcat, shrimpfaoreg)
f2shrimp <- fig2distplot(dtshrimp)
f2shrimp

## other demersal
demfaoreg <- c(21,31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
demspcat <- c(31,33,34) # Manual match, i.e. parameters from literature. Same for all.
dtdem <- stockselectfig2(upsides, demspcat, demfaoreg)
f2dem <- fig2distplot(dtdem)
f2dem

## tuna, billfish, shark
tunafaoreg <- c(21,31,27,37) #Pelagic impacts assumed to also include NE Atlantic and Med. (based on Wallace et al. 2010 RMUs)
tunaspcat <- c(36,38)
dttuna <- stockselectfig2(upsides, tunaspcat, tunafaoreg)
f2tuna <- fig2distplot(dttuna)
f2tuna
