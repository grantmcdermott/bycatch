rm(list = ls())

library(dplyr)
library(readr)

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

###############
###############
### Turtles ###
###############
###############

##########################################
## Loggerhead sea turtles - NW Atlantic ##
##########################################

loggerhead_list <- 
  list(
    Demersals =
      list(type = 1,
           wt = 0.25,
           faoreg = c(21, 31), # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
           spcat = c(31, 33, 34),
           countries = NA,
           lumpedstocks = NA), 
    Tuna =
      list(type = 1,
           wt = 0.07,
           faoreg = c(21, 27, 31, 37), #Pelagic impacts assumed to also include NE Atlantic and Med. (based on Wallace et al. 2010 RMUs)
           spcat = c(36, 38),
           countries = NA,
           lumpedstocks = NA), 
    Shrimp =
      list(type = 1,
           wt = 0.68,
           faoreg = c(21, 31), # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
           spcat = 45,
           countries = NA,
           lumpedstocks = NA)
  )

## Convert to data frame (and remove list)
loggerhead_df <-
  lapply(loggerhead_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Loggerhead turtle") %>%
  select(species, everything())
rm(loggerhead_list)


################################
## Leatherback - East Pacific ##
################################

leatherback_list <- 
  list(
    Demersals =
      list(type = 1,
           wt = 0.84,
           faoreg = c(77, 87), 
           spcat = c(31, 33, 34),
           countries = NA,
           lumpedstocks = NA), 
    Tuna =
      list(type = 1,
                wt = 0.16,
                faoreg = c(77, 87),
                spcat = c(36, 38),
                countries = NA,
                lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
leatherback_df <-
  lapply(leatherback_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Leatherback turtle") %>%
  select(species, everything())
rm(leatherback_list)


####################################
## Olive ridley - NE Indian Ocean ##
####################################

olive_nei_list <- 
  list(
    target_species = ## info? 
      list(type = 1, 
           wt = 1,
           faoreg = 57, 
           spcat = c(31, 32, 33, 34, 38, 45),
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
olive_nei_df <-
  lapply(olive_nei_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Olive ridley turtle (NEI)") %>%
  select(species, everything())
rm(olive_nei_list)

####################################
## Olive ridley - W Indian Ocean ##
####################################

olive_wi_list <- 
  list(
    target_species = ## info?
      list(type = 1, 
           wt = 1,
           faoreg = 51, 
           spcat = c(31, 32, 33, 34, 38, 45),
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
olive_wi_df <-
  lapply(olive_wi_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Olive ridley turtle (WI)") %>%
  select(species, everything())
rm(olive_wi_list)


###################
##### Mammals #####
###################

#########################
## Australian sea lion ##
#########################

aus_sealion_list <- 
  list(
    Demersals =
      list(type = 2, 
           wt = 0.25,
           faoreg = NA, 
           spcat = c(31, 32, 33, 34, 45), #Trawl - secondary threat (assume 25%)
           countries = "Australia",
           lumpedstocks = NA),
    Sharks =
      list(type = 2, 
           wt = 0.75,
           faoreg = NA, 
           spcat = 38, #Shark gillnet - major threat (assume 75%))
           countries = "Australia",
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
aus_sealion_df <-
  lapply(aus_sealion_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Australian sea lion") %>%
  select(species, everything())
rm(aus_sealion_list)

#################
## NZ Sea lion ##
#################

nz_sealion_list <- 
  list(
    target_species = ## info?
      list(type = 1, 
           wt = 1,
           faoreg = c(81, 88), 
           spcat = 57, 
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
nz_sealion_df <-
  lapply(nz_sealion_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "NZ sea lion") %>%
  select(species, everything())
rm(nz_sealion_list)


###################################
## Finless porpoise - NW Pacific ##
###################################

finless_list <- 
  list(
    target_species =
      list(type = 1, 
           wt = 1,
           faoreg = 61, 
           spcat = c(31, 32, 33, 34), 
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
finless_df <-
  lapply(finless_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Finless porpoise") %>%
  select(species, everything())
rm(finless_list)


###################################
## Humpback dolphin - SW Pacific ##
###################################

humpback_list <- 
  list(
    target_species = ## info?
      list(type = 5, 
           wt = 1,
           faoreg = 61, 
           spcat = c(31, 32, 33, 34, 38), 
           countries = c("China", "Taiwan Province of China"), 
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
humpback_df <-
  lapply(humpback_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Humpback dolphin") %>%
  select(species, everything())
rm(humpback_list)


###########################################
## Vaquita - Gulf of California - Mexico ##
###########################################

vaquita_list <- 
  list(
    target_species = ## info?
      list(type = 1, 
           wt = 1,
           faoreg = 77, 
           spcat = c(31, 32, 33, 34),
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
vaquita_df <-
  lapply(vaquita_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Vaquita") %>%
  select(species, everything())
rm(vaquita_list)



#################
##### Birds #####
#################

#################################################
## Amsterdam albatross - Southern Indian Ocean ##
#################################################

amsterdam_list <- 
  list(
    target_species = ## info?
      list(type = 1, 
           wt = 1,
           faoreg = 57, 
           spcat = 36,
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
amsterdam_df <-
  lapply(amsterdam_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Amsterdam albatross") %>%
  select(species, everything())
rm(amsterdam_list)


##########################################################################
## Sooty shearwater - Atlantic, Pacific, S Indian Ocean, Southern Ocean ##
##########################################################################

# 600% reduction needed

sooty_list <- 
  list(
    target_species = ## info?
      list(type = 1, 
           wt = 1,
           faoreg = c(61, 71, 81, 67, 77, 87, 21, 27, 31, 34, 41, 47, 48, 58, 88), 
           spcat = 36,
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
sooty_df <-
  lapply(sooty_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Sooty shearwater") %>%
  select(species, everything())
rm(sooty_list)


####################################
## Tristan albatross - S Atlantic ##
####################################

tristan_list <- 
  list(
    target_species = ## info?
      list(type = 1, 
           wt = 1,
           faoreg = c(41, 47), 
           spcat = 36,
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
tristan_df <-
  lapply(tristan_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "Tristan albatross") %>%
  select(species, everything())
rm(tristan_list)


############################################################################
## White-chinned petrel - S Pacific, S Atlantic, S Indian, Southern Ocean ##
############################################################################

petrel_list <- 
  list(
    target_species_1 = ## info?
      list(type = 1, 
           wt = 0.167,
           faoreg = c(48, 58), 
           spcat = 32,
           countries = NA,
           lumpedstocks = NA),
    target_species_2 = ## info?
      list(type = 1, 
           wt = 0.167,
           faoreg = 47, 
           spcat = 32,
           countries = NA,
           lumpedstocks = NA),
    target_species_3 = ## info?
      list(type = 1, 
           wt = 0.167,
           faoreg = 41, 
           spcat = 32,
           countries = NA,
           lumpedstocks = NA),
    target_species_4 = ## info?
      list(type = 1, 
           wt = 0.167,
           faoreg = 47, 
           spcat = 36,
           countries = NA,
           lumpedstocks = NA),
    target_species_5 = ## info?
      list(type = 1, 
           wt = 0.167,
           faoreg = 41, 
           spcat = c(36, 38),
           countries = NA,
           lumpedstocks = NA),
    target_species_6 = ## info?
      list(type = 1, 
           wt = 0.167,
           faoreg = c(71, 81), 
           spcat = 36,
           countries = NA,
           lumpedstocks = NA)
    )

## Convert to data frame (and remove list)
petrel_df <-
  lapply(petrel_list, expand.grid) %>%
  bind_rows(.id = "target") %>%
  mutate(species = "White-chinned petrel") %>%
  select(species, everything())
rm(petrel_list)



#######################
#######################
##### Combine all #####
#######################
#######################

target_species <- bind_rows(mget(ls()))
write_csv(target_species, "Data/target_species.csv")
