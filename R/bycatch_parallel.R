## Clear environment
rm(list = ls())

##################################
########### Data Setup ###########
##################################
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(pbapply)
library(magrittr)
library(stringr)
library(pbmcapply)
library(snow)
library(snowfall)

# Date
run_date <- "20170123" #gsub("-", "", Sys.Date())

# Make directory to store the results
run_dir = paste0('results/', run_date)

if(dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
}

# Start time
start.time <- Sys.time()

###################################
########## Load functions #########
###################################

source("R/bycatch_funcs_cost.R")


##############################
########## Load data #########
##############################

## Load bycatch data
bycatch_df <- read_csv("Data/bycatch_species.csv")
target_df <- read_csv("Data/target_species.csv")

### Load target data, derived from the "upsides" model of Costello et al. (PNAS, 2016)
## Uncertainty, no nei stocks version
upsides <- read_csv("Data/bycatch-upuncert-input.csv", col_types = cols(regionfao = "c"))
# upsides_fread <- data.table::fread("Data/bycatch-upuncert-input.csv") %>% as_data_frame() ## quicker
#
## No uncertainty, nei stocks version
upsides <- read_csv("Data/bycatch-nouncert-input.csv", col_types = cols(regionfao = "c"))

###############################
########### Results ###########
###############################

## Sampling parameters
n1 <- 100 # Run n = 1000 in 10 chunks of 100
n2 <- 100
sens_exp <- 1 # 1 is the normal value (i.e. same as main analysis). 
# Please also run the no-uncert analysis once each with sens_exp <- 0.5 and sens_exp <- 2. Thanks!

# n1 <- 10 
# n2 <- 10

######## TIME CONSUMING PART ######################
all_species_samp <- bycatch_df$species 

# Track errors
options(error = traceback)

# Initialize cluster
sfInit(parallel = do.parallel, cpus = NumCPUs)

# Functions and parameters needed
sfExportAll()

# Source functions
sfSource("R/bycatch_funcs.R")
# sfSource("bycatch_funcs_cost_uncert - Copy.R")

# Load packages on all cores
sfLibrary(tidyr)
sfLibrary(dplyr)

# All species results (run 10 times)
all_samp1 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp2 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp3 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp4 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp5 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp6 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp7 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp8 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp9 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp10 <- bind_rows(tryCatch(sfLapply(all_species_samp, bycatch_func), error = function(e) NULL))
all_samp <- bind_rows(all_samp1, all_samp2, all_samp3, 
                      all_samp4, all_samp5, all_samp6, 
                      all_samp7, all_samp8, all_samp9, 
                      all_samp10)
                       
# Stop parallel processing
sfStop()

### ORIGINAL CODE #####
# # All species results
# all_species_samp <- bycatch_df$species 
# all_samp1 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp2 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp3 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp4 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp5 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp6 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp7 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp8 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp9 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp10 <- bind_rows(pblapply(all_species_samp, bycatch_func))
# all_samp <- bind_rows(all_samp1, all_samp2, all_samp3, 
#                       all_samp4, all_samp5, all_samp6, 
#                       all_samp7, all_samp8, all_samp9, 
#                       all_samp10
# )

# End time
print(Sys.time() - start.time)

# Save bycatch results
write_csv(all_samp, path = paste0(run_dir, "/bycatch_results_", run_date, ".csv"))

############################################################################################################################

# Remove samples that are no longer needed
rm(all_samp1,all_samp2,all_samp3,all_samp4,all_samp5,all_samp6,all_samp7,all_samp8,all_samp9,all_samp10)

# Read in bycatch results file from above
all_dt <- read_csv(paste0(run_dir, "/bycatch_results_", run_date, ".csv"))

alldistplots <- bycatchdistggplot(all_dt) +
  facet_wrap(~species, ncol = 3, scales = "free")
allcostplots <- costggplot(all_dt) +
  facet_wrap(~species, ncol = 3, scales = "free")
alldistplots
allcostplots


### Get median, 2.5th, 25th, 75th, 97.5th percentiles from runs for Figs. 3a,b
results_summary <- all_dt %>%
  group_by(species) %>%
  summarise(pctredmsy025 = quantile(pctredmsy, probs = 0.025),
            pctredmsy25 = quantile(pctredmsy, probs = 0.25),
            pctredmsy50 = median(pctredmsy),
            pctredmsy75 = quantile(pctredmsy, probs = 0.75),
            pctredmsy975 = quantile(pctredmsy, probs = 0.975),
            pctredmey025 = quantile(pctredmey, probs = 0.025),
            pctredmey25 = quantile(pctredmey, probs = 0.25),
            pctredmey50 = median(pctredmey),
            pctredmey75 = quantile(pctredmey, probs = 0.75),
            pctredmey975 = quantile(pctredmey, probs = 0.975),
            ycostmsy025 = quantile(ycostmsy, probs = 0.025),
            ycostmsy25 = quantile(ycostmsy, probs = 0.25),
            ycostmsy50 = median(ycostmsy),
            ycostmsy75 = quantile(ycostmsy, probs = 0.75),
            ycostmsy975 = quantile(ycostmsy, probs = 0.975),
            pcostmey025 = quantile(pcostmey, probs = 0.025),
            pcostmey25 = quantile(pcostmey, probs = 0.25),
            pcostmey50 = median(pcostmey),
            pcostmey75 = quantile(pcostmey, probs = 0.75),
            pcostmey975 = quantile(pcostmey, probs = 0.975))
results_summary <- left_join(results_summary, bycatch_df, by = 'species')

# Save results
write_csv(results_summary, paste("bycatch_results_summary_", run_date, ".csv", sep = ""))

print(Sys.time() - start.time)


################################################################################################################


## Turtle results
turtle_species_samp <- (filter(bycatch_df, grp=="turtle"))$species #c("Loggerhead turtle", "Olive ridley turtle (NEI)")
## Run the bycatch function over all turtle species
turtles_samp1 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp2 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp3 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp4 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp5 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp6 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp7 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp8 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp9 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp10 <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtles_samp <- bind_rows(turtles_samp1, turtles_samp2, turtles_samp3, 
                          turtles_samp4, turtles_samp5#, turtles_samp6, 
                          #turtles_samp7, turtles_samp8, turtles_samp9, 
                          #turtles_samp10
                          )
#write_csv(turtles_samp, "turtles_results.csv") 
turtledistplots <- bycatchdistggplot(turtles_samp1) +
  facet_wrap(~species, scales = "free")
turtlecostplots <- costggplot(turtles_samp1) +
  facet_wrap(~species, scales = "free")
turtledistplots
turtlecostplots

## Mammal results
mammal_species_samp <- (filter(bycatch_df, grp=="mammal"))$species 
## Run the bycatch function over all mammal species
mammals_samp1 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp2 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp3 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp4 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp5 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp6 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp7 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp8 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp9 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp10 <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammals_samp <- bind_rows(mammals_samp1, mammals_samp2, mammals_samp3, 
                          mammals_samp4, mammals_samp5#, mammals_samp6, 
                          #mammals_samp7, mammals_samp8, mammals_samp9, 
                          #mammals_samp10
                          )
#write_csv(mammals_samp, "mammals_results.csv")
mammaldistplots <- bycatchdistggplot(mammals_samp1) +
  facet_wrap(~species, scales = "free")
mammalcostplots <- costggplot(mammals_samp1) +
  facet_wrap(~species, scales = "free")
mammaldistplots
mammalcostplots

## Bird results
bird_species_samp <- (filter(bycatch_df, grp=="bird"))$species 
## Run the bycatch function over all bird species
bird_samp1 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp2 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp3 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp4 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp5 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp6 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp7 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp8 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp9 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp10 <- bind_rows(pblapply(bird_species_samp, bycatch_func))
bird_samp <- bind_rows(bird_samp1, bird_samp2, bird_samp3, 
                       bird_samp4, bird_samp5, bird_samp6, 
                       bird_samp7, bird_samp8, bird_samp9, 
                       bird_samp10)
#write_csv(bird_samp, "bird_results.csv")
birddistplots <- bycatchdistggplot(bird_samp1) +
  facet_wrap(~species, scales = "free")
birdcostplots <- costggplot(bird_samp1) +
  facet_wrap(~species, scales = "free")
birddistplots
birdcostplots

# vaquita test
vaquit <- bycatch_func("Vaquita porpoise")
vext <- extract_func("Vaquita porpoise")
upsamp <- upsides_subset_func(vext$Totoaba)

bycatchdistggplot(vaquit) 
costggplot(vaquit) 

# Tristan albatross
talb <- bycatch_func("Tristan albatross")

bycatchdistggplot(vaquit) 
costggplot(vaquit) 

# Loggerhead test
lh <- bycatch_func("Loggerhead turtle")
lhext <- extract_func("Loggerhead turtle")

lh <- all_dt %>%
  filter(species == "Loggerhead turtle (NW Atlantic)")
bycatchdistggplot(lh) 
costggplot(lh) 

# Leatherback test
lb <- bycatch_func("Leatherback turtle")
bycatchdistggplot(lb) 
costggplot(lb) 

# Australian sea lion test
as <- bycatch_func("Australian sea lion")
bycatchdistggplot(as) 
costggplot(as) 

# NZ sea lion test
ns <- bycatch_func("NZ sea lion")
bycatchdistggplot(ns) 
costggplot(ns) 


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
  plot(bd$MEY, bd$wt, xlim=c(-50,100), col=2, pch = 19,
       xlab = "% Reduction in Mortality", ylab = "2012 Effort ($)")
  #legend("topleft", legend=names(den), fill=c(4,6))
}

## shrimp
shrimpfaoreg <- c(31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
shrimpspcat <- c(45)
dtshrimp <- stockselectfig2(upsides, shrimpspcat, shrimpfaoreg)
f2shrimp <- fig2distplot(dtshrimp)
f2shrimp

## other demersal
demfaoreg <- c(31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
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


########################################################
########### Overall Reductions (MEY and MSY) ###########
########################################################

ovrred <- upsides %>%
  group_by(idoriglumped) %>%
  summarise(margc = mean(marginalcost),
            bet = mean(beta),
            g = mean(g),
            fvfmey = mean(eqfvfmey),
            fvfmsy = mean(fvfmsy),
            pctmey = mean(pctredfmey),
            pctmsy = mean(pctredfmsy)) %>%
  mutate(wt = margc * ((g * fvfmsy)^bet),
         cstcurr = wt,
         cstmey = margc * (((g * fvfmsy)/fvfmey)^bet),
         cstmsy = margc * ((g)^bet)) %>%
  ungroup() %>%
  mutate(wt = wt/sum(wt, na.rm = T)) %>%
  mutate(wtpctmey = wt * pctmey,
         wtpctmsy = wt * pctmsy) 

avpctmey <- sum(ovrred$wtpctmey, na.rm = T)
avpctmsy <- sum(ovrred$wtpctmsy, na.rm = T)

avpctmey <- 100 * (1 - (sum(ovrred$cstmey, na.rm = T)/sum(ovrred$cstcurr, na.rm = T)))
avpctmsy <- 100 * (1 - (sum(ovrred$cstmsy, na.rm = T)/sum(ovrred$cstcurr, na.rm = T)))

## List of species categories
# list("Shads" = 24, "Flounders, halibuts, soles" = 31, 
#   "Cods, hakes, haddocks" = 32,"Miscellaneous coastal fishes" = 33,
#  "Miscellaneous demersal fishes" = 34,"Herrings, sardines,anchovies" = 35,
# "Tunas,bonitos,billfishes" = 36,"Miscellaneous pelagic fishes" = 37,
#"Sharks, rays, chimeras" = 38,"Shrimps, prawns" = 45,
#"Carps, barbels and other cyprinids" = 11,"Sturgeons, paddlefishes" = 21,
#"Salmons, trouts, smelts" = 23,"Miscellaneous diadromous fishes" = 25,
#"Crabs, sea-spiders" = 42,"Lobsters, spiny rock lobsers" = 43,
#"King crabs, squat lobsters" = 44,"Miscellaneous marine crustaceans" = 47,
#"Abalones, winkles, conchs" = 52,"Oysters" = 53,"Mussels" = 54,
#"Scallops, pectens" = 55,"Clams, cockles, arkshells" = 56,
#"Squids, cuttlefishes, octopuses" = 57,"Horseshoe crabs and other arachnoids" = 75,
#"Sea-urchins and other echinoderms" = 76,"Miscellaneous aquatic invertebrates" = 56)

##########################################################
########### Checking species-level assumptions ###########
##########################################################

# Loggerhead
lhext <- extract_func("Loggerhead turtle")
lhshrimp <- upsides %>%
  filter(regionfao %in% lhext$Shrimp$faoreg) %>% 
  filter(speciescat %in% lhext$Shrimp$spcat) %>%
  mutate(maxmp = mprofitf(0,price,marginalcost,g,k,phi,beta))
lhshrimpstocks <- lhshrimp %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

lhdem <- upsides %>%
  filter(regionfao %in% lhext$Demersals$faoreg) %>% 
  filter(speciescat %in% lhext$Demersals$spcat) %>%
  mutate(maxmp = mprofitf(0,price,marginalcost,g,k,phi,beta))
lhdemstocks <- lhdem %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

lhtuna <- upsides %>%
  filter(regionfao %in% lhext$Tuna$faoreg) %>% 
  filter(speciescat %in% lhext$Tuna$spcat) %>%
  mutate(maxmp = mprofitf(0,price,marginalcost,g,k,phi,beta)) %>%
  filter(marginalcost > 0)
lhtunastocks <- lhtuna %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

# Leatherback
lbext <- extract_func("Leatherback turtle")

lbdem <- upsides %>%
  filter(regionfao %in% lbext$Demersals$faoreg) %>% 
  filter(speciescat %in% lbext$Demersals$spcat) %>%
  filter(country != "USA")
lbdemstocks <- lbdem %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

lbtuna <- upsides %>%
  filter(regionfao %in% lbext$Tuna$faoreg) %>% 
  filter(speciescat %in% lbext$Tuna$spcat) %>%
  filter(country != "USA")
lbtunastocks <- lbtuna %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

# Leatherback
lbext <- extract_func("Leatherback turtle")

lbdem <- upsides %>%
  filter(regionfao %in% lbext$Demersals$faoreg) %>% 
  filter(speciescat %in% lbext$Demersals$spcat)
lbdemstocks <- lbdem %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

lbtuna <- upsides %>%
  filter(regionfao %in% lbext$Tuna$faoreg) %>% 
  filter(speciescat %in% lbext$Tuna$spcat)
lbtunastocks <- lbtuna %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

# Olive ridley - NE Indian Ocean
orext <- extract_func("Olive ridley turtle (NEI)")

ordem <- upsides %>%
  filter(regionfao %in% orext$target_species$faoreg) %>% 
  filter(speciescat %in% orext$target_species$spcat) %>%
  filter(country %in% orext$target_species$countries)
ordemstocks <- ordem %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

# Olive ridley - W Indian Ocean
orext <- extract_func("Olive ridley turtle (WI)")

ordem <- upsides %>%
  filter(regionfao %in% orext$target_species$faoreg) %>% 
  filter(speciescat %in% orext$target_species$spcat)
ordemstocks <- ordem %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

# Australian sea lion
asext <- extract_func("Australian sea lion")

assh <- upsides %>%
  filter(regionfao %in% asext$Sharks$faoreg) %>% 
  filter(speciescat %in% asext$Sharks$spcat) %>%
  filter(country %in% asext$Sharks$countries)
ordemstocks <- ordem %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

# New Zealand sea lion
nsext <- extract_func("NZ sea lion")

nsoth <- upsides %>%
  filter(regionfao %in% nsext$Other$faoreg) %>% 
  filter(speciescat %in% nsext$Other$spcat) %>%
  filter(country %in% nsext$Other$countries)
nsssqstocks <- nsoth %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

# vaquita
vext <- extract_func("Vaquita")
upsamp <- upsides_subset_func(vext$Totoaba)
upsampd <- upsides_subset_func(vext$Demersals)
upsampd <- upsides %>%
  filter(regionfao %in% vext$Demersals$faoreg) %>% 
  filter(speciescat %in% vext$Demersals$spcat) %>%
  filter(country %in% vext$Demersals$countries)

#### Testing subset function for bycatch species
 
test <- extract_func(c('Tristan albatross'))
testdt <- upsides_subset_func(test$All)
teststocks <- testdt %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))
