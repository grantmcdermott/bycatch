##################################
########### Data Setup ###########
##################################

#rm(list = ls())
require(readr)
require(dplyr)
require(tidyr)
require(stringr)
require(grid)
library(tibble)
library(magrittr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtools)
library(graphics)
#Shrink full datasets to needed columns/rows
load("ProjectionData Data.rdata")
rm(OriginalProjectionData)
upsmallu <- as_data_frame(UnlumpedProjectionData) %>% 
  select(Year,IdOrig,CommName,SciName,Country,g,phi,k,MSY,MarginalCost,SpeciesCat,SpeciesCatName,RegionFAO,Policy,beta,FvFmsy,BvBmsy,Catch,Price) %>%
  filter(Year %in% c(2012,2050))
upsmalll <- as_data_frame(ProjectionData) %>% 
  select(Year,IdOrig,CommName,SciName,Country,g,phi,k,MSY,MarginalCost,SpeciesCat,SpeciesCatName,RegionFAO,Policy,beta,FvFmsy,BvBmsy,Catch,Price) %>%
  filter(Year %in% c(2012,2050))
upuncert <- as_data_frame(read.csv("upsidessamples.csv",stringsAsFactors=F))

#Make small files for IdOrig lookup table
upidslumped <- select(upsmalll,Year,IdOrig,CommName,SciName,Country,SpeciesCat,SpeciesCatName,RegionFAO,Catch) %>%
  filter(Year %in% c(2012))
upidsunlumped <- select(upsmallu,Year,IdOrig,CommName,SciName,Country,SpeciesCat,SpeciesCatName,RegionFAO,Catch) %>%
  filter(Year %in% c(2012))
write.csv(upidslumped,"UpIDsLumped.csv",row.names = F)
write.csv(upidsunlumped,"UpIDsUnlumped.csv",row.names = F)

ullookup <- as_data_frame(read.csv("upsidesidlookup.csv",stringsAsFactors=F))

#Ensure all files are data frames and convert column names to lowercase 
names(upsmalll) %<>% tolower
names(upsmallu) %<>% tolower
names(upuncert) %<>% tolower
names(ullookup) %<>% tolower

## Convert upsides columns to bycatch-relevant ones (e.g. cost) in shrunken dataset and uncertainty samples

#Rename columns to make lumped/unlumped tables compatible
upuncert <- upuncert %>%
  rename(idoriglumped = idorig) %>%
  rename(gsamp = g) %>%
  rename(ksamp = k) %>%
  rename(msysamp = msy) %>%
  rename(phisamp = phi) %>%
  select(-x)

upsmallu <- 
  left_join(
    upsmallu,
    ullookup
  )

# Remove unnecessary files
rm(ProjectionData,UnlumpedProjectionData,upsmalll)

# Make new unlumped table with MSY, MEY and Historic split out 
uphist <-
  upsmallu %>%
  ungroup() %>%
  filter(year %in% c(2012)) %>%
  filter(policy %in% c("Historic")) %>%
  select(-year)
head(uphist)

upmey <-
  upsmallu %>%
  ungroup() %>%
  filter(year %in% c(2050)) %>%
  filter(policy %in% c("Opt")) %>%
  select(idorig,fvfmsy,bvbmsy) %>%
  rename(fmeyvfmsy = fvfmsy) %>%
  rename(bmeyvbmsy = bvbmsy)
head(upmey)

upcatchshare <-
  upsmallu %>%
  ungroup() %>%
  filter(year %in% c(2050)) %>%
  filter(policy %in% c("CatchShare")) %>%
  select(idorig,fvfmsy,bvbmsy) %>%
  rename(fcsvfmsy = fvfmsy) %>%
  rename(bcsvbmsy = bvbmsy)
head(upcatchshare)

restabnouncert <- 
  left_join(
  left_join(uphist,
            upmey),
  upcatchshare)
head(restabnouncert)

# Remove unnecessary files
rm(ullookup,upsmallu,uphist,upmey,upcatchshare)

#Combine uncertainty and unlumped table
upcomb <- 
  left_join(
    restabnouncert %>%
      ungroup(),
    upuncert
  )

#Create new columns for results table
upcomb <- upcomb %>%
  mutate(g = ifelse(!is.na(gsamp), gsamp, g)) %>%
  mutate(k = ifelse(!is.na(ksamp), ksamp, k)) %>%
  mutate(phi = ifelse(!is.na(phisamp), phisamp, phi)) %>%
  mutate(msy = ifelse(!is.na(msysamp), msysamp, msy)) %>%
  mutate(fvfmsy = ifelse(!is.na(finalfvfmsy), finalfvfmsy, fvfmsy)) %>%
  mutate(bvbmsy = ifelse(!is.na(finalbvbmsy), finalbvbmsy, bvbmsy)) %>%
  mutate(costmsy = marginalcost*(g^beta)*(catch/lumpedcatch2012)) %>%
  mutate(costmey = marginalcost*(fmeyvfmsy^beta)*(g^beta)*(catch/lumpedcatch2012)) %>%
  mutate(costcs = 0.77*marginalcost*(fcsvfmsy^beta)*(g^beta)*(catch/lumpedcatch2012)) %>%
  mutate(costhist = marginalcost*(fvfmsy^beta)*(g^beta)*(catch/lumpedcatch2012)) %>%
  mutate(pricecs = 1.31*price) %>%
  mutate(marginalcostcs = 0.77*marginalcost) %>%
  mutate(fvfmey = fvfmsy/fmeyvfmsy) %>%
  mutate(fvfcs = fvfmsy/fcsvfmsy) %>%
  mutate(pctchangemsy = 100*(1 - (1/fvfmsy))) %>%
  mutate(pctchangemey = 100*(1 - (1/fvfmey))) %>%
  mutate(pctchangecs = 100*(1 - (1/fvfcs)))
head(upcomb)

upcomb <- 
  upcomb %>%
  group_by(idorig) %>%
  mutate(cnt = n()) %>%
  mutate(repmult = 1001 - cnt)
head(upcomb)

# Remove unnecessary files
rm(upuncert)

# Create results table
restab <- 
  upcomb %>%
  select(-gsamp,-ksamp,-msysamp,-phisamp,-finalfvfmsy,-finalbvbmsy,-policy,-repmult)
head(restab)

# Remove unnecessary files
rm(upcomb)

# Export final results input file
write.csv(restab,"BycatchResultsInputFile2016.csv",row.names = F)
write.csv(restabnouncert,"BycatchResultsInputFileNoUncert2016.csv",row.names = F)
# restab <- as_data_frame(read.csv("BycatchResultsInputFile2016.csv",stringsAsFactors=F))
# restabnouncert <- as_data_frame(read.csv("BycatchResultsInputFileNoUncert2016.csv",stringsAsFactors=F))

##################################
############ Results #############
##################################


###################################################################
########## Make table samples, functions that Grant needs #########
###################################################################

tabsampf1 <- 
  restabnouncert %>%
  mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
  mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
  mutate(pctredfmey = 100 * (1 - (fmeyvfmsy/fvfmsy))) %>%
  head(20) %>%
  select(idorig,pctredfmsy,pctredfmey,wt)

wtsum <- sum(tabsampf1$wt,na.rm = TRUE)

tabsampf1 <- 
  tabsampf1 %>%
  mutate(weight = wt/wtsum) %>%
  select(-wt)

write.csv(tabsampf1,"sample-table-f1.csv",row.names = F)

tabsampf2 <- 
  restabnouncert %>%
  mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
  mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
  mutate(pctredfmey = 100 * (1 - (fmeyvfmsy/fvfmsy))) %>%
  head(20) %>%
  select(idorig,pctredfmsy,pctredfmey,wt,fvfmsy,g,beta,phi,price,marginalcost)

wtsum <- sum(tabsampf2$wt,na.rm = TRUE)

tabsampf2 <- 
  tabsampf2 %>%
  mutate(weight = wt/wtsum) %>%
  select(-wt)

write.csv(tabsampf2,"sample-table-f2.csv",row.names = F)

eqprofit <- function(pctred,price,marginalcost,fvfmsy,g,k,phi,beta) {
  eqp <-  (0.01 * price * fvfmsy * g * k * (100 - pctred) * (((1 + phi - (fvfmsy * phi) + (0.01 * fvfmsy * pctred * phi))/(1 + phi))^(1/phi))) - 
    (marginalcost * ((fvfmsy * g * (1 - (pctred/100)))^beta)) + 
   return(eqp)
}

eqyield <- function(pctred,fvfmsy,g,k,phi) {
  eqy <- 0.01 * fvfmsy * g * k * (100 - pctred) * (((1 + phi - (fvfmsy * phi) + (0.01 * fvfmsy * pctred * phi))/(1 + phi))^(1/phi))
  return(eqy)
}

bycatchbenefit <- function(pctred,weight) {
  benefit <- pctred * weight
  return(benefit)
}

##############################################
########## Results functions #########
##############################################

##### pseudocode
#1. Subset full table into small table with weights, stock ids and pctred values (stockselect functions below)
#2. For each table, create a sampled version with n rows--which samples from pctred 
#     (i.e. a vector of pctred samples of length n)--combine MSY and MEY vectors into data frame
#3. Take weighted averages of resulting vectors for different species groups

##### end pseudocode

### Function 1: Subset full table into table with weights (4 versions) "stockselect"

# Function needed below to handle stocks with multiple fao regions
faofun <- function(x) paste("grepl(", toString(x), ",","dt$regionfao)", sep = "")

# Selecting stocks: option 1: list of species categories, fao regions specified
stockselect1 <- function(dt,spcat,faoreg){
  # Construct command to filter FAO regions
  a <- lapply(faoreg,faofun)
  b <- paste(a, collapse = "|")
  dtuse <- 
    dt %>%
    filter((eval(parse(text = b)))) %>% 
    filter(speciescat %in% spcat) %>%
    filter(fmeyvfmsy > 0) %>%
    mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
    mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
    mutate(pctredfmey = 100 * (1 - (fmeyvfmsy/fvfmsy))) %>%
    select(idorig,pctredfmsy,pctredfmey,wt,fvfmsy,g,beta,phi,price,marginalcost)
  return(dtuse)
}

# Selecting stocks: option 2: list of species categories, countries specified
stockselect2 <- function(dt,spcat,countries){
  # Construct command to filter FAO regions
  dtuse <- 
    dt %>%
    filter(speciescat %in% spcat) %>%
    filter(country %in% countries) %>%
    filter(fmeyvfmsy > 0) %>%
    mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
    mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
    mutate(pctredfmey = 100 * (1 - (fmeyvfmsy/fvfmsy))) %>%
    select(idorig,pctredfmsy,pctredfmey,wt,fvfmsy,g,beta,phi,price,marginalcost)
  return(dtuse)
}

# Selecting stocks: option 3: list of lumped stocks, fao regions specified
stockselect3 <- function(dt,lumpedstocks,faoreg){
  # Construct command to filter FAO regions
  a <- lapply(faoreg,faofun)
  b <- paste(a, collapse = "|")
  dtuse <- 
    dt %>%
    filter((eval(parse(text = b)))) %>% 
    filter(idoriglumped %in% lumpedstocks) %>%
    filter(fmeyvfmsy > 0) %>%
    mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
    mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
    mutate(pctredfmey = 100 * (1 - (fmeyvfmsy/fvfmsy))) %>%
    select(idorig,pctredfmsy,pctredfmey,wt,fvfmsy,g,beta,phi,price,marginalcost)
  return(dtuse)
} 

# Selecting stocks: option 4: list of lumped stocks, countries specified
stockselect4 <- function(dt,lumpedstocks,countries){
  # Construct command to filter FAO regions
  dtuse <- 
    dt %>%
    filter(idoriglumped %in% lumpedstocks) %>%
    filter(country %in% countries) %>%
    filter(fmeyvfmsy > 0) %>%
    mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
    mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
    mutate(pctredfmey = 100 * (1 - (fmeyvfmsy/fvfmsy))) %>%
    select(idorig,pctredfmsy,pctredfmey,wt,fvfmsy,g,beta,phi,price,marginalcost)
  return(dtuse)
} 


### Function 2: Turns subsettted table with pctred and wt into sample of pctred

bycdist <- function(dt,n1,n2) {
  sumwt <- sum(dt$wt)
  dt1 <- dt %>%
    mutate(wt = wt/sumwt)
  sampmsy <- c()
  sampmey <- c()
  for (i in 1:n1) {
    sampmsy <- append(sampmsy, mean(sample(dt1$pctredfmsy, size = n2, replace = T, prob = dt1$wt)))
    sampmey <- append(sampmey, mean(sample(dt1$pctredfmey, size = n2, replace = T, prob = dt1$wt)))
  }
  samp <- data.frame(pctredmsy=sampmsy,pctredmey=sampmey)
  return(samp)
} 

### Function 3: Plot that asks: Will rebuilding be enough to stop decline? 

bycatchdistplot <- function(bdist,pctrdbpt,pctrdbl,pctrdbu) {
  bd <- bdist %>%
    rename(MSY = pctredmsy,
           MEY = pctredmey)
  den <- apply(bd, 2, density)
  plot(NA, xlim=range(sapply(den, "[", "x")), 
       ylim=range(sapply(den, "[", "y")), 
       xlab = "% Reduction in Mortality", ylab = "Density")
  rect(pctrdbl,-1,pctrdbu,1,col = 'lightgrey',border = NA)
  mapply(lines, den, col=c(4,6))
  abline(v=pctrdbpt,col = 1, lty = 2)
  legend("topleft", legend=names(den), fill=c(4,6))
}

#Version 2 fixes the range at zero to 100
bycatchdistplot2 <- function(bdist,pctrdbpt,pctrdbl,pctrdbu) {
  bd <- bdist %>%
    rename(MSY = pctredmsy,
           MEY = pctredmey)
  den <- apply(bd, 2, density)
  plot(NA, xlim=c(0,100), 
       ylim=range(sapply(den, "[", "y")), 
       xlab = "% Reduction in Mortality", ylab = "Density")
  rect(pctrdbl,-1,pctrdbu,1,col = 'lightgrey',border = NA)
  mapply(lines, den, col=c(4,6))
  abline(v=pctrdbpt,col = 1, lty = 2)
  legend("topleft", legend=names(den), fill=c(4,6))
}

#Version 3 fixes the range at -500 to 100 (for NZ sea lion)
bycatchdistplot3 <- function(bdist,pctrdbpt,pctrdbl,pctrdbu) {
  bd <- bdist %>%
    rename(MSY = pctredmsy,
           MEY = pctredmey)
  den <- apply(bd, 2, density)
  plot(NA, xlim=c(-500,100), 
       ylim=range(sapply(den, "[", "y")), 
       xlab = "% Reduction in Mortality", ylab = "Density")
  rect(pctrdbl,-1,pctrdbu,1,col = 'lightgrey',border = NA)
  mapply(lines, den, col=c(4,6))
  abline(v=pctrdbpt,col = 1, lty = 2)
  legend("topleft", legend=names(den), fill=c(4,6))
}

#Version 4 fixes the range at -20 to 100 (for Finless porpoise)
bycatchdistplot4 <- function(bdist,pctrdbpt,pctrdbl,pctrdbu) {
  bd <- bdist %>%
    rename(MSY = pctredmsy,
           MEY = pctredmey)
  den <- apply(bd, 2, density)
  plot(NA, xlim=c(-50,100), 
       ylim=range(sapply(den, "[", "y")), 
       xlab = "% Reduction in Mortality", ylab = "Density")
  rect(pctrdbl,-1,pctrdbu,1,col = 'lightgrey',border = NA)
  mapply(lines, den, col=c(4,6))
  abline(v=pctrdbpt,col = 1, lty = 2)
  legend("topleft", legend=names(den), fill=c(4,6))
}

#Version 5 puts legend on right (for Au sea lion)
bycatchdistplot5 <- function(bdist,pctrdbpt,pctrdbl,pctrdbu) {
  bd <- bdist %>%
    rename(MSY = pctredmsy,
           MEY = pctredmey)
  den <- apply(bd, 2, density)
  plot(NA, xlim=c(0,100), 
       ylim=range(sapply(den, "[", "y")), 
       xlab = "% Reduction in Mortality", ylab = "Density")
  rect(pctrdbl,-1,pctrdbu,1,col = 'lightgrey',border = NA)
  mapply(lines, den, col=c(4,6))
  abline(v=pctrdbpt,col = 1, lty = 2)
  legend("topright", legend=names(den), fill=c(4,6))
}

#Version 6 allows range to be specified
bycatchdistplot6 <- function(bdist,pctrdbpt,pctrdbl,pctrdbu,rl,ru) {
  bd <- bdist %>%
    rename(MSY = pctredmsy,
           MEY = pctredmey)
  den <- apply(bd, 2, density)
  plot(NA, xlim=c(rl,ru), 
       ylim=range(sapply(den, "[", "y")), 
       xlab = "% Reduction in Mortality", ylab = "Density")
  rect(pctrdbl,-1,pctrdbu,1,col = 'lightgrey',border = NA)
  mapply(lines, den, col=c(4,6))
  abline(v=pctrdbpt,col = 1, lty = 2)
  legend("topleft", legend=names(den), fill=c(4,6))
}

### Function 4: How much does it cost to reduce mortality enough? (minimized) 

## Inputs: 
# 1: Filtered dataframe, "datafilt"
# 2: "pctredbpt", "pctredbl", "pctredbu"
# 3: A desired percent chance that the bycatch reduction threshold is met, "pctchance"
# 4: A function ("eqprofit") describing the relationship for a target stock between the equilibrium profit 
#     and the percent reduction in F (pctredf; a variable) in terms of columns in the input table.
#     The function will be the same for all stocks.

# 5: A function ("eqyield") describing the relationship for a target stock between the equilibrium yield (i.e. catch) 
#     and the percent reduction in F (pctredf; a variable) in terms of columns in the input table.
#     The function will be the same for all stocks.

# 6: A function ("bycatch benefit") describing the relationship between bycatch benefit (in terms of percent reduction in bycatch mortality)
#     and the percent reduction in F for a target stock (pctredf; a variable)


## Outputs: 
# 1. 2 values (not a data frame): (i) minimum profit cost of achieving the point/median percent change needed for bycatch species
#                                 (ii) minimum yield cost of achieving the point/median percent change needed for bycatch species
# 2,3. Same as 1. but for each of the lower and upper percent change needed.

# mydf %>%
# ungroup() %>% ## just in case
#   arrange(margcost)

# mydf %>% mutate(cum_margcost = cumsum(margcost))

######################
### Function tests ###
######################

## Inputs ##

# Estimates of percent change needed for the bycatch species of interest
pctredbpt <- 56 # Point estimae
pctredbl  <- 35 # Lower estimate
pctredbu  <- 78 # Upper estimate

# Desired percent chance that the bycatch reduction threshold is met
pctchance <- 95

# Relevant taxonomic groups
stocklist <- c("Lumped-Atlantic mackerel-FaoRegion37","Lumped-Atlantic saury-FaoRegion37") # Relevant lumped stocks, if given
cntrylist <- c("Greece","Italy","France") # Relevant countries, if given
faoreglist <- c(21,31) # Relevant FAO fishing areas, if given
spcatlist <- c(31,33,34) # Relevant species categories, if given

# Test stockselect functions
stockselect1(restabnouncert,spcatlist,faoreglist)
stockselect2(restabnouncert,spcatlist,cntrylist)
stockselect3(restabnouncert,stocklist,faoreglist)
stockselect4(restabnouncert,stocklist,cntrylist)

### Tests with sample table 
t1 <- stockselect1(restabnouncert,spcatlist,faoreglist)
t2 <- bycdist(t1,100,100)
bycatchdistplot(t2,pctredbpt,pctredbl,pctredbu)

# Combining multiple sub-groups with given weights (e.g. shrimps, tunas, etc. for NW Atl. loggerhead)
n1 <- 100
n2 <- 100
demfaoreg <- c(21,31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
demspcat <- c(31,33,34)
demwt <- 0.25
demres <- stockselect1(restabnouncert,demspcat,demfaoreg)
demdist <- bycdist(demres,n1,n2)

tunafaoreg <- c(21,31,27,37) #Pelagic impacts assumed to also include NE Atlantic and Med. (based on Wallace et al. 2010 RMUs)
tunaspcat <- c(36,38)
tunawt <- 0.07
tunares <- stockselect1(restabnouncert,tunaspcat,tunafaoreg)
tunadist <- bycdist(tunares,n1,n2)

shrimpfaoreg <- c(21,31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
shrimpspcat <- c(45)
shrimpwt <- 0.68
shrimpres <- stockselect1(restabnouncert,shrimpspcat,shrimpfaoreg) 
shrimpdist <- bycdist(shrimpres,n1,n2)

loggerheadnwadist <- (demwt * demdist) + (shrimpwt * shrimpdist) + (tunawt * tunadist)

bycatchdistplot(loggerheadnwadist,pctredbpt,pctredbl,pctredbu)  
  
###############################
########### Results ###########
###############################

# Desired percent chance that the bycatch reduction threshold is met
pctchance <- 95
n1 <- 10000
n2 <- 100

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



##########################################
## Loggerhead sea trutles - NW Atlantic ##
##########################################

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.98
fe <- 0.067
  
pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller

demfaoreg <- c(21,31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
demspcat <- c(31,33,34)
demwt <- 0.25
demres <- stockselect1(restabnouncert,demspcat,demfaoreg)
demdist <- bycdist(demres,n1,n2)

tunafaoreg <- c(21,31,27,37) #Pelagic impacts assumed to also include NE Atlantic and Med. (based on Wallace et al. 2010 RMUs)
tunaspcat <- c(36,38)
tunawt <- 0.07
tunares <- stockselect1(restabnouncert,tunaspcat,tunafaoreg)
tunadist <- bycdist(tunares,n1,n2)

shrimpfaoreg <- c(21,31) # Demersal impacts assumed to be restricted to US and Caribbean (NMFS 2008)
shrimpspcat <- c(45)
shrimpwt <- 0.68
shrimpres <- stockselect1(restabnouncert,shrimpspcat,shrimpfaoreg) 
shrimpdist <- bycdist(shrimpres,n1,n2)

loggerheadnwadist <- (demwt * demdist) + (shrimpwt * shrimpdist) + (tunawt * tunadist)
rm(demdist,tunadist,shrimpdist)

#bycatchdistplot(loggerheadnwadist,pctredbpt,pctredbl,pctredbu)
bycatchdistplot2(loggerheadnwadist,pctredbpt,pctredbl,pctredbu)

loggerheadnwaplot1 <- bycatchdistplot(loggerheadnwadist,pctredbpt,pctredbl,pctredbu)

################################
## Leatherback - East Pacific ##
################################

# 83% reduction needed - Pelagic longlines (16%), Coastal driftnets (84%) (Kaplan et al. 2005 for breakdown)

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.97
fe <- 0.036
lambda2 <- 0.846
fe2 <- 0.334

pctredbpt <- 100 * ((1 - lambda2)/fe2) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda2))/(fe2 * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda2))/(fe2 * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller

faoreg <- c(77,87) 
spgr <- c(31,33,34)
demwt <- 0.84
demresLBT <- stockselect1(restabnouncert,spgr,faoreg)
demdistLBT <- bycdist(demresLBT,n1,n2)

faoreg <- c(77,87)  
spgr <- c(36,38) #Longline
tunawt <- 0.16
tunaresLBT <- stockselect1(restabnouncert,spgr,faoreg)
tunadistLBT <- bycdist(tunaresLBT,n1,n2)

LBTdist <- (demwt * demdistLBT) + (tunawt * tunadistLBT)
rm(demdistLBT,tunadistLBT)

#bycatchdistplot(LBTdist,pctredbpt,pctredbl,pctredbu)
bycatchdistplot2(LBTdist,pctredbpt,pctredbl,pctredbu)

LBTplot1 <- bycatchdistplot(LBTdist,pctredbpt,pctredbl,pctredbu)

####################################
## Olive ridley - NE Indian Ocean ##
####################################

# 23% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.97
fe <- 0.0871

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller

faoreg <- c(57) 
spgr <- c(31,32,33,34,38,45)
ORTEIres <- stockselect1(restabnouncert,spgr,faoreg)

ORTEIdist <- bycdist(ORTEIres,n1,n2)

#bycatchdistplot(ORTEIdist,pctredbpt,pctredbl,pctredbu)
bycatchdistplot2(ORTEIdist,pctredbpt,pctredbl,pctredbu)

ORTEIplot1 <- bycatchdistplot(ORTEIres,pctredbpt,pctredbl,pctredbu)

####################################
## Olive ridley - W Indian Ocean ##
####################################

# 129% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.82
fe <- 0.14

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller

faoreg <- c(51) 
spgr <- c(31,32,33,34,38,45)
ORTWIres <- stockselect1(restabnouncert,spgr,faoreg)

ORTWIdist <- bycdist(ORTWIres,n1,n2)

#bycatchdistplot(ORTWIdist,pctredbpt,pctredbl,pctredbu)
bycatchdistplot2(ORTWIdist,pctredbpt,pctredbl,pctredbu)

###################
##### Mammals #####
###################

#########################
## Australian sea lion ##
#########################

# 6.5% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.9923
fe <- 0.1177

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller


cntry <- c("Australia") 
spgr <- c(31,32,33,34,45) #Trawl - secondary threat (assume 25%)
demwt <- 0.25
demresAuslM <- stockselect2(restabnouncert,spgr,cntry)
demdistAuslM <- bycdist(demresAuslM,n1,n2)

spgr <- c(38) #Shark gillnet - major threat (assume 75%)
sharkwt <- 0.75
sharkresAuslM <- stockselect2(restabnouncert,spgr,cntry)
sharkdistAuslM <- bycdist(sharkresAuslM,n1,n2)

AuslMdist <- (demwt * demdistAuslM) + (sharkwt * sharkdistAuslM)
rm(demdistAuslM,sharkdistAuslM)

#bycatchdistplot(AuslMdist,pctredbpt,pctredbl,pctredbu)
bycatchdistplot5(AuslMdist,pctredbpt,pctredbl,pctredbu)

#################
## NZ Sea lion ##
#################

# 59% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.96
fe <- 0.06734

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller

faoreg <- c(81,88) 
spgr <- c(57)
NZslMres <- stockselect1(restabnouncert,spgr,faoreg)
NZslMdist <- bycdist(NZslMres,n1,n2)

bycatchdistplot3(NZslMdist,pctredbpt,pctredbl,pctredbu)

###################################
## Finless porpoise - NW Pacific ##
###################################

# 57% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.947
fe <- 0.093

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller

#faoreg <- c(61) # Change to China/Taiwan? Need to update
cntry <- c("China","Taiwan Province of China") 
spgr <- c(31,32,33,34)
FPMres <- stockselect2(restabnouncert,spgr,cntry)

FPMdist <- bycdist(FPMres,n1,n2)

bycatchdistplot4(FPMdist,pctredbpt,pctredbl,pctredbu)

###################################
## Humpback dolphin - SW Pacific ##
###################################

# 38% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.9754
fe <- 0.0646

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller

faoreg <- c(61) 
spgr <- c(31,32,33,34,11)
HDMres <- stockselect1(restabnouncert,spgr,faoreg)

HDMdist <- bycdist(HDMres,n1,n2)

bycatchdistplot2(HDMdist,pctredbpt,pctredbl,pctredbu)

###########################################
## Vaquita - Gulf of California - Mexico ##
###########################################

# 76% reduction needed
# Totoaba - Cr. EN

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.87
fe <- 0.17

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller


faoreg <- c(77) 
spgr <- c(31,32,33,34)
VaMres <- stockselect1(restabnouncert,spgr,faoreg)

VaMdist <- bycdist(VaMres,n1,n2)

bycatchdistplot2(VaMdist,pctredbpt,pctredbl,pctredbu)


#################
##### Birds #####
#################


#################################################
## Amsterdam albatross - Southern Indian Ocean ##
#################################################

# 53% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.953
fe <- 0.088

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller


faoreg <- c(57) 
spgr <- c(36)
AABres <- stockselect1(restabnouncert,spgr,faoreg)

AABdist <- bycdist(AABres,n1,n2)

bycatchdistplot6(AABdist,pctredbpt,pctredbl,pctredbu,-100,100)

##########################################################################
## Sooty shearwater - Atlantic, Pacific, S Indian Ocean, Southern Ocean ##
##########################################################################

# 600% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.982
fe <- 0.003

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller


faoreg <- c(61,71,81,67,77,87,21,27,31,34,41,47,48,58,88) 
spgr <- c(36)
SSBres <- stockselect1(restabnouncert,spgr,faoreg)

SSBdist <- bycdist(SSBres,n1,n2)

bycatchdistplot6(SSBdist,pctredbpt,pctredbl,pctredbu,-100,100)

####################################
## Tristan albatross - S Atlantic ##
####################################

# 91% reduction needed

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.971
fe <- 0.032

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller


faoreg <- c(41,47) 
spgr <- c(36)
TABres <- stockselect1(restabnouncert,spgr,faoreg)

TABdist <- bycdist(TABres,n1,n2)

bycatchdistplot6(TABdist,pctredbpt,pctredbl,pctredbu,-100,100)

############################################################################
## White-chinned petrel - S Pacific, S Atlantic, S Indian, Southern Ocean ##
############################################################################

# 28% reduction needed 

# Estimates of percent change needed for the bycatch species of interest
lambda <- 0.964
fe <- 0.13

pctredbpt <- 100 * ((1 - lambda)/fe) # Point estimae
pctredbl  <- 100 * ((0.75 * (1 - lambda))/(fe * 1.25)) # Lower estimate, assumes decline is 25% smaller, fe 25% larger
pctredbu  <- 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper estimate, assumes decline is 25% larger, fe 25% smaller


faoreg <- c(48,58) 
spgr <- c(32)
wt1 <- 0.167
combresWCPB1 <- stockselect1(restabnouncert,spgr,faoreg)
WCPBdist1 <- bycdist(combresWCPB1,n1,n2)

faoreg <- c(47) 
spgr <- c(32)
wt2 <- 0.167
combresWCPB2 <- stockselect1(restabnouncert,spgr,faoreg)
WCPBdist2 <- bycdist(combresWCPB2,n1,n2)

faoreg <- c(41) 
spgr <- c(32)
wt3 <- 0.167
combresWCPB3 <- stockselect1(restabnouncert,spgr,faoreg)
WCPBdist3 <- bycdist(combresWCPB3,n1,n2)

faoreg <- c(47) 
spgr <- c(36)
wt4 <- 0.167
combresWCPB4 <- stockselect1(restabnouncert,spgr,faoreg)
WCPBdist4 <- bycdist(combresWCPB4,n1,n2)

faoreg <- c(41) 
spgr <- c(36,38)
wt5 <- 0.167
combresWCPB5 <- stockselect1(restabnouncert,spgr,faoreg)
WCPBdist5 <- bycdist(combresWCPB5,n1,n2)

faoreg <- c(71,81) 
spgr <- c(36)
wt6 <- 0.167
combresWCPB6 <- stockselect1(restabnouncert,spgr,faoreg)
WCPBdist6 <- bycdist(combresWCPB6,n1,n2)

WCPBdist <- (WCPBdist1 * wt1) + (WCPBdist2 * wt2) + (WCPBdist3 * wt3) + (WCPBdist4 * wt4) + (WCPBdist5 * wt5) + (WCPBdist6 * wt6)
rm(WCPBdist1,WCPBdist2,WCPBdist3,WCPBdist4,WCPBdist5,WCPBdist6)

bycatchdistplot5(WCPBdist,pctredbpt,pctredbl,pctredbu)

#####################
##### Figure 1 ######
#####################

turtleplot1 <- grid.arrange(loggerheadnwaplot1,LBTplot1,ORTWIplot1,ORTEIplot1,ncol = 2)
mammalplot1 <- grid.arrange(AuslMplot1,FPMplot1,HDMplot1,VaMplot1,ncol = 2)
birdplot1 <- grid.arrange(AABplot1,SSBplot1,TABplot1,WCPBplot1,ncol = 2)

grid.arrange(loggerheadnwaplot1,LBTplot1,ORTWIplot1,ORTEIplot1,ncol = 2)
grid.arrange(AuslMplot1,FPMplot1,HDMplot1,VaMplot1,ncol = 2)
grid.arrange(AABplot1,SSBplot1,TABplot1,WCPBplot1,ncol = 2)


#####################
######## End ########
#####################