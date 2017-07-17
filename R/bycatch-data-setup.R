rm(list = ls())
library(magrittr)
library(splitstackshape)
library(tidyverse)
library(stringr)
library(cowplot)
library(pbapply)
library(pbmcapply)

###########################################################
############ Load Costello et al. projections #############
###########################################################

upsidesunlumped <- read_csv("Data/TBD_IF_NEEDED/UnlumpedProjectionData.csv", col_types = cols(RegionFAO = "c")) # UnlumpedProjectionData is by stock and country
names(upsidesunlumped) %<>% tolower
upsmallu <- 
  upsidesunlumped %>% 
  select(year,idorig,commname,sciname,country,g,phi,k,msy,c,speciescat,speciescatname,regionfao,policy,fvfmsy,bvbmsy,catch,price) %>%
  filter(year %in% c(2012)) %>%
  rename(marginalcost = c) %>%
  mutate(beta = 1.3)


################################################################
############ Set up data file for bycatch analysis #############
################################################################

upsides_kobe <- read_csv("Data/TBD_IF_NEEDED/Kobe_MEY.csv") %>%
  rename(idoriglumped = IdOrig,
         eqfvfmey = current_f_mey,
         eqfmeyvfmsy = f_mey) %>%
  select(idoriglumped, eqfvfmey, eqfmeyvfmsy)

## First read in a lookup table to map idoriglumped column to upsmallu DF
idlookup <- read_csv("Data/TBD_IF_NEEDED/upsidesidlookup.csv", col_types = "ccd")
names(idlookup) %<>% tolower

# load unlumped table
upsides <- 
  upsmallu %>%
  left_join(idlookup, by = 'idorig')
rm(upsmallu)

upsides <- left_join(upsides, upsides_kobe, by = 'idoriglumped') %>%
  mutate(curr_f = g * fvfmsy,
         f_mey = g * eqfmeyvfmsy,
         pctredfmsy = 100 * (1 - (1/fvfmsy)),
         pctredfmey = 100 * (1 - (1/eqfvfmey))) %>%
  rename(fmeyvfmsy = eqfmeyvfmsy) %>%
  select(idorig,idoriglumped,commname,sciname,country,speciescat,speciescatname,regionfao,
         k,fvfmsy,bvbmsy,g,beta,phi,price,eqfvfmey,fmeyvfmsy,curr_f,f_mey,pctredfmsy,pctredfmey)


## 3-YEAR MOVING AVERAGES FOR F (skip if using 2012 only as base year)
gm_mean = function(x, na_rm = TRUE, zero_propagate = TRUE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero_propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na_rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na_rm) / length(x))
  }
}

up3yr <- 
  upsidesunlumped %>% 
  select(idorig, year, fvfmsy) %>% 
  filter(year %in% 2010:2012) %>% 
  group_by(idorig) %>% 
  summarise(fvfmsy3yr = gm_mean(fvfmsy))

upsides <- 
  upsides %>%
  left_join(up3yr, by = 'idorig') %>%
  mutate(
    fvfmsy = fvfmsy3yr,
    curr_f = fvfmsy * g,
    eqfvfmey = curr_f/f_mey,
    pctredfmsy = 100 * (1-(g/curr_f)),
    pctredfmey = 100 * (1-(f_mey/curr_f))
    ) %>%
  select(-fvfmsy3yr)
rm(up3yr)
# END 3 yr moving average

# Add conservation concern scenario columns: fconmsy, fconmey
upsides <- upsides %>%
  mutate(bmsy = k * ((1/(1 + phi))^(1/phi)),
         curr_b = bvbmsy * bmsy, # calculate bvbmey, which is needed to calculate fconmey below
         bmey = (k * ((1-((f_mey * phi)/(g * (phi + 1))))^(1/phi))), 
         bvbmey = curr_b/bmey,
         fconmsy = if_else((fvfmsy < 1) & (bvbmsy > 1), curr_f, g), # calculate fconmsy
         fconmey = if_else((eqfvfmey < 1) & (bvbmey > 1), curr_f, f_mey), # calculate fconmey
         pctredfmsycon = 100 * (1-(fconmsy/curr_f)),
         pctredfmeycon = 100 * (1-(fconmey/curr_f))) %>%
  select(-bvbmsy, -bvbmey, -curr_b, -bmey, -bmsy)

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
                     f_mey = 0.055069581,
                     fconmsy = 0.057
) %>%
  mutate(pctredfmsy = 100 * (1-(g/curr_f)),
         pctredfmey = 100 * (1-(f_mey/curr_f)),
         pctredfmsycon = 100 * (1-(fconmsy/curr_f)))

upsides <- bind_rows(upsides,totoab)

## recalculate unlumped marginalcosts based on f_mey estimates from Kobe file
source("R/bycatch_funcs.R")

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
rm(upsides2,dtup,mcup,idsup,upsides_kobe,totoab)
## end clean up

## Add a database ID (i.e. if stock assessment source was FAO or RAM legacy)
dbaseid <- read_csv("Data/TBD_IF_NEEDED/ram_stock_lkup.csv") ## already loaded
names(dbaseid) %<>% tolower
upsides <-
  upsides %>%
  left_join(dbaseid %>% select(-lumpedid)) %>%
  select(dbase, everything()) %>%
  mutate(dbase = ifelse(is.na(dbase), "FAO", dbase)) %>%
  mutate(dbase = ifelse(idorig=="toto", "Totoaba", dbase)) ## Totoaba was added manually (not part of upsides)

rm(dbaseid,idlookup,upsidesunlumped)

# Write final main input file
write_csv(upsides, "Data/upsides_main.csv")

# If using only 2012 as base year for F (instead of 3-yr moving averages): Write final input file
write_csv(upsides, "Data/upsides_2012only.csv")

## OPTIONAL: Remove 'nei' stocks
upsides <-
  upsides %>%
  mutate(nei = grepl(' nei',commname)) %>%
  filter(nei == FALSE) %>%
  select(-nei) 

# Write final input file without nei stocks
write_csv(upsides, "Data/upsides_nonei.csv")

############ Diagnostics #############
max((upsides %>% filter(fmeyvfmsy > 0))$fmeyvfmsy) ## no fmeyvfmsy
hist((upsides %>%
        filter(fmeyvfmsy > 0))$fmeyvfmsy,
     breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,30000),
     xlim = c(0,2))
######## end diagnostics
