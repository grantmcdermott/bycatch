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

###################################
########## Load functions #########
###################################

source("bycatch_funcs_cost_2.R")

#################################################
############ Read in/clean the data #############
#################################################

## Load upsides (i.e. target species projection) data from Costello et al. 2016
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
         pctredfmey = 100 * (1 - (1/eqfvfmey))) %>%
  select(idorig,idoriglumped,commname,sciname,country,speciescat,speciescatname,regionfao,
         k,fvfmsy,g,beta,phi,price,eqfvfmey,fmeyvfmsy,curr_f,f_mey,pctredfmsy,pctredfmey)

# Add Totoaba row to upsides (for vaquita)
totoab <- data_frame(idorig = "toto",
                     idoriglumped = "totoaba",
                     commname = "Totoaba",
                     sciname = "Totoaba macdonaldi",
                     country = "Mexico",
                     speciescat = 0, # N/A
                     speciescatname = "N/A",
                     regionfao = "77",
                     k = 15825,
                     fvfmsy = 0.52631579,
                     g = 0.057,
                     beta = 1.3,
                     phi = 0.188,
                     price = 50000,
                     eqfvfmey = 0.52631579,
                     fmeyvfmsy = 0.966133,
                     curr_f = 0.03,
                     f_mey = 0.055069581
) %>%
  mutate(pctredfmsy = 100 * (1-(g/curr_f)),
         pctredfmey = 100 * (1-(f_mey/curr_f)))

upsides <- bind_rows(upsides,totoab)
## recalculate unlumped marginalcosts based on f_mey estimates from Kobe file

# calculation function
margcost <- function(fmey,price,g,k,phi,beta) {
  mc <- uniroot((function (x) mprofitf(fmey,price,x,g,k,phi,beta)), 
                lower = -1000000000, upper = 100000000000)[1]$root
}

# calculation
upsides2 <- upsides %>%
  filter(f_mey > 0)
mcup <- c()
idsup <- c()
for (i in 1:length(upsides2$idorig)) {
  idsup <- append(idsup, upsides2$idorig[i])
  mcup <- append(mcup, margcost(upsides2$f_mey[i],upsides2$price[i],upsides2$g[i],upsides2$k[i],upsides2$phi[i],upsides2$beta[i]))
}
dtup <- data_frame(idorig = idsup, marginalcost = mcup) %>%
  filter(marginalcost > 0)
upsides <- left_join(upsides, dtup, by = 'idorig')
## end marginal cost calculation 

## clean up
rm(upsides2,dtup,mcup,idsup,upsides_kobe)
## end clean up


## Load bycatch data
bycatch_df <- read_csv("Data/bycatch_species.csv")

target_df <- read_csv("Data/target_species.csv")


###############################
########### Results ###########
###############################

## Sampling parameters
n1 <- 50
n2 <- 100

## Turtle results
turtle_species_samp <- (filter(bycatch_df, grp=="turtle"))$species #c("Loggerhead turtle", "Olive ridley turtle (NEI)")
## Run the bycatch function over all turtle species
turtles_samp <- bind_rows(pblapply(turtle_species_samp, bycatch_func))
turtledistplots <- bycatchdistggplot(turtles_samp) +
  facet_wrap(~species, scales = "free")
turtlecostplots <- costggplot(turtles_samp) +
  facet_wrap(~species, scales = "free")
#write_csv(turtles_samp, "turtles_test.csv")
turtledistplots
turtlecostplots

## Mammal results
mammal_species_samp <- (filter(bycatch_df, grp=="mammal"))$species 
## Run the bycatch function over all mammal species
mammals_samp <- bind_rows(pblapply(mammal_species_samp, bycatch_func))
mammaldistplots <- bycatchdistggplot(mammals_samp) +
  facet_wrap(~species, scales = "free")
mammalcostplots <- costggplot(mammals_samp) +
  facet_wrap(~species, scales = "free")
#write_csv(mammals_samp, "mammals_test.csv")
mammaldistplots
mammalcostplots

## Bird results
bird_species_samp <- (filter(bycatch_df, grp=="bird"))$species 
## Run the bycatch function over all bird species
bird_samp <- bind_rows(pblapply(bird_species_samp, bycatch_func))
birddistplots <- bycatchdistggplot(bird_samp) +
  facet_wrap(~species, scales = "free")
birdcostplots <- costggplot(bird_samp) +
  facet_wrap(~species, scales = "free")
#write_csv(bird_samp, "bird_test.csv")
birddistplots
birdcostplots

# vaquita test
vaquit <- bycatch_func("Vaquita")
vext <- extract_func("Vaquita")
upsamp <- upsides_subset_func(vext$Totoaba)

bycatchdistggplot(vaquit) 
costggplot(vaquit) 


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
  filter(speciescat %in% lhext$Shrimp$spcat)
lhshrimpstocks <- lhshrimp %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

lhdem <- upsides %>%
  filter(regionfao %in% lhext$Demersals$faoreg) %>% 
  filter(speciescat %in% lhext$Demersals$spcat)
lhdemstocks <- lhdem %>%
  group_by(idoriglumped) %>%
  summarize(meanfvfmsy = mean(fvfmsy))

lhtuna <- upsides %>%
  filter(regionfao %in% lhext$Tuna$faoreg) %>% 
  filter(speciescat %in% lhext$Tuna$spcat)
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
  filter(speciescat %in% orext$target_species$spcat)
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

# vaquita
vext <- extract_func("Vaquita")
upsamp <- upsides_subset_func(vext$Totoaba)
upsampd <- upsides_subset_func(vext$Demersals)
upsampd <- upsides %>%
  filter(regionfao %in% vext$Demersals$faoreg) %>% 
  filter(speciescat %in% vext$Demersals$spcat) %>%
  filter(country %in% vext$Demersals$countries)

  