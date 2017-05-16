#rm(list = ls())
# library(readr)
# library(dplyr)
# library(tidyr)
# library(tibble)
# library(ggplot2)
library(splitstackshape)
library(tidyverse)
library(stringr)
library(magrittr)
library(cowplot)
library(pbapply)
library(pbmcapply)

#################################################################
############ Uncertainty & no-nei stocks data setup #############
#################################################################

## Load Dan Ovando's C-MSY samples for lumped, non-NEI, non-RAM stocks in Costello et al. 2016.
upuncert <- read_csv("Data/TBD_IF_NEEDED/upsidessamples.csv")
names(upuncert) %<>% tolower

# Load in upsides data (from Costello et al. 2016) and clean it.
upsideslumped <- read_csv("Data/TBD_IF_NEEDED/ProjectionData.csv", col_types = cols(RegionFAO = "c")) # ProjectionData is by stock
upsidesunlumped <- read_csv("Data/TBD_IF_NEEDED/UnlumpedProjectionData.csv", col_types = cols(RegionFAO = "c")) # UnlumpedProjectionData is by stock and country
upsmalll <- 
  upsideslumped %>% 
  select(Year,IdOrig,CommName,SciName,Country,g,phi,k,MSY,c,SpeciesCat,SpeciesCatName,RegionFAO,Policy,FvFmsy,BvBmsy,Catch,Price) %>%
  filter(Year %in% c(2012)) %>%
  rename(marginalcost = c) %>%
  mutate(beta = 1.3)
upsmallu <- 
  upsidesunlumped %>% 
  select(Year,IdOrig,CommName,SciName,Country,g,phi,k,MSY,c,SpeciesCat,SpeciesCatName,RegionFAO,Policy,FvFmsy,BvBmsy,Catch,Price) %>%
  filter(Year %in% c(2012)) %>%
  rename(marginalcost = c) %>%
  mutate(beta = 1.3)
names(upsmallu) %<>% tolower
upidslumped <- 
  select(upsmalll,Year,IdOrig,CommName,marginalcost,Price,g,phi,k,beta) %>%
  filter(Year %in% c(2012))
names(upidslumped) %<>% tolower
rm(upsmalll,upsidesunlumped,upsideslumped)

# calculate 'b_oa', 'f_oa', remove marginalcost from upidslumped
#         Costello et al. SI pg 8: MSY = ((g * k)/((phi + 1)^(1/phi)))
#         Costello et al. SI eq. 10: marginalcost = (price * f_oa * b_oa * MSY)/((g *f_oa)^beta) 
#         From Costello et al. SI eq. 5, at equilibrium: f_oa =  (((phi + 1)/phi) * (1 - ((b_oa^phi)/(phi + 1))))
#         Therefore: marginalcost = (price * f_oa * b_oa * MSY)/((g *f_oa)^beta)
# calculate b_oa from equations above:
num_cores <- min(24, detectCores()) # The specified number or the number of CPUs your computer has, whichever is smaller.
test <- 
  # pbmclapply(
  pblapply(
    1:length(upidslumped$idorig),
    # mc.cores = num_cores,
    function(i){
      
      ids <- upidslumped$idorig[i]
      
      boas <- uniroot(f = function(x){
        upidslumped$marginalcost[i] - 
          ((upidslumped$price[i] * 
              (((upidslumped$phi[i] + 1)/upidslumped$phi[i]) * 
                 (1 - ((x^upidslumped$phi[i])/(upidslumped$phi[i] + 1)))) * 
              x * 
              ((upidslumped$g[i] * upidslumped$k[i])/((upidslumped$phi[i] + 1)^(1/upidslumped$phi[i]))))/
          ((upidslumped$g[i] *(((upidslumped$phi[i] + 1)/upidslumped$phi[i]) * 
                              (1 - ((x^upidslumped$phi[i])/
                                      (upidslumped$phi[i] + 1)))))^upidslumped$beta[i]))
      },
                     interval = c(0,1),
      extendInt="yes"
      )$root
      df = data_frame(idorig = ids, b_oa = boas)
      return(df)
    },
    cl = num_cores) %>%
  bind_rows
# merge 'b_oa' column into upidslumped
upidslumped <- left_join(upidslumped,test,by = 'idorig')
rm(test)

# calculate 'f_oa' from 'b_oa'
upidslumped <- 
  upidslumped %>%
  mutate(f_oa = (((phi + 1)/phi) * (1 - ((b_oa^phi)/(phi + 1)))))
upidslumped <- 
  upidslumped %>%
  select(idorig,commname,price,b_oa,f_oa,k) %>%
  rename(klumped = k)

# join price, 'b_oa', 'f_oa' columns to upuncert (the bootstrapped dataset)
upuncert <- left_join(upuncert,upidslumped,by='idorig')

# Add beta column and rename fvfmsy column
upuncerttest <- 
  upuncert %>%
  mutate(beta = 1.3) %>%
  rename(
    fvfmsy = finalfvfmsy,
    bvbmsy = finalbvbmsy
    )
upuncert <- upuncerttest
rm(upuncerttest)

# calculate marginal costs for bootstraps
upuncerttest <- 
  upuncert %>%
  mutate(marginalcost = (price * f_oa * b_oa * ((g * k)/((phi + 1)^(1/phi))))/((g *f_oa)^beta))
upuncert <- upuncerttest
rm(upuncerttest)

# Save intermediate file (so don't have to calculate fmeys again)
# write.csv(upuncert,"bycatch-upuncert-input.csv",row.names = F)

# load intermediate file if starting new session
# upuncert <- read_csv("bycatch-upuncert-input.csv")

# Calculate fmey's for each bootstrap
num_cores <- min(24, detectCores()) # The specified number or the number of CPUs your computer has, whichever is smaller.
test <- 
  pblapply(
  # pbmclapply(
    1:length(upuncert$idorig),
    # mc.cores = num_cores,
    function(i){
      
      ids <- upuncert$idorig[i]
      
      eqfvfmeys <- 
        1/(optim(par = 0.005,  
                 fn = function(x){
                   - ((upuncert$price[i] * upuncert$g[i] * x * upuncert$k[i] * 
                         ((1 - ((upuncert$g[i] * x * upuncert$phi[i])/(upuncert$g[i] * (upuncert$phi[i] + 1))))^(1/upuncert$phi[i]))) - 
                        (upuncert$marginalcost[i] * ((upuncert$g[i] * x)^upuncert$beta[i])))
                 }, 
                 method = "Brent", 
                 lower = 0, 
                 upper = 6.32)$par * # Upper bound is equal to (phi + 1)/phi, which is the maximum growth rate as a fraction of fmsy
             (1/upuncert$fvfmsy[i]))
      
      df = data_frame(idorig = ids, eqfvfmey = eqfvfmeys)
      return(df)
    },
    cl = num_cores) %>%
  bind_rows
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 01h 34m 42s

# Save intermediate file (so don't have to calculate fmeys again)
# write_csv(test, "Data/TBD_IF_NEEDED/fmeys-uncert.csv")

# load 'test' if starting new session
# test <- read_csv("Data/TBD_IF_NEEDED/fmeys-uncert.csv")

# Merge fmeys into original table
upuncert$eqfvfmey <- test$eqfvfmey
rm(test)

upuncert <- 
  upuncert %>%
  mutate(
    fmeyvfmsy = fvfmsy/eqfvfmey,
    curr_f = g * fvfmsy,
    f_mey = g * fmeyvfmsy
    ) %>%
  mutate(
    pctredfmsy = 100 * (1 - (1/fvfmsy)),
    pctredfmey = 100 * (1 - (1/eqfvfmey))
    )

# Save intermediate file (so don't have to calculate fmeys again)
# write_csv(upuncert, "Data/TBD_IF_NEEDED/bycatch-upuncert-input.csv")

# load 'upuncert' if starting new session
# upuncert <- read_csv("Data/TBD_IF_NEEDED/bycatch-upuncert-input.csv")

upuncert <- 
  upuncert %>%
  rename(idoriglumped = idorig) 

# Create unlumped version of table
# create repeated rows in the country table

## First read in a lookup table to map idoriglumped column to upsmallu DF
idlookup <- read_csv("Data/TBD_IF_NEEDED/upsidesidlookup.csv", col_types = "ccd")
names(idlookup) %<>% tolower

# load unlumped table
upsides <- 
  upsmallu %>%
  left_join(idlookup) %>%
  select(idorig,idoriglumped,commname,sciname,country,speciescat,speciescatname,regionfao,g,k) %>%
  rename(
    g2 = g,
    k2 = k
    )

# Remove all nei stocks
test <- 
  upidslumped %>%
  ungroup() %>%
  mutate(nei = grepl(' nei',commname)) %>%
  filter(nei == FALSE) %>%
  rename(idoriglumped = idorig) %>%
  select(idoriglumped, nei) 

upsides <- left_join(upsides, test, by = 'idoriglumped')
upsides <- 
  upsides %>%
  filter(nei == FALSE) # 2781 non-nei lumped stocks
rm(test, idlookup)

############ Diagnostics #############
median((upuncert %>% filter(fmeyvfmsy > 0))$g)
max((upuncert %>% filter(fmeyvfmsy > 0))$fmeyvfmsy)
max((upsides %>% filter(fmeyvfmsy > 0))$fmeyvfmsy) ## no fmeyvfmsy
head((upuncert %>% filter(fmeyvfmsy > 0))$fmeyvfmsy)
hist((upuncert %>%
        filter(fmeyvfmsy > 0))$fmeyvfmsy,
     breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,30000),
     xlim = c(0,7))
hist((upsides %>%
        filter(fmeyvfmsy > 0))$fmeyvfmsy,
     breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,30000),
     xlim = c(0,2))
length((upuncert %>% filter(fmeyvfmsy > 0))$g)
######## end diagnostics

# Number upuncert rows and create stock-by-replicate column
repno <- c(rep(seq(1:1000),3415)) 
upuncert$repnumber <- repno
rm(repno)
test <- upuncert %>%
  mutate(idlumped_n_rep = paste(idoriglumped,repnumber,sep = "-"))
upuncert <- test
rm(test)

# Save intermediate file (so don't have to calculate fmeys again)
# write_csv(upuncert, "Data/TBD_IF_NEEDED/bycatch-upuncert-input.csv")

# load 'upuncert' if starting new session
# upuncert <- read_csv("Data/TBD_IF_NEEDED/bycatch-upuncert-input.csv")

# Replicate upsides rows 1000 times and create stock-by-replicate column
test <- upsides
for (i in 1:10) {
  test <- rbind(test,test)
}
upsides <- test
rm(test)
repno <- sort(c(rep(seq(1:1024),5400)))
upsides$repnumber <- repno
upsides <- upsides %>%
  mutate(idlumped_n_rep = paste(idoriglumped,repnumber,sep = "-")) %>%
  filter(repnumber < 1001)
rm(repno)

# Merge upuncert columns into upsides
upsides <- left_join(upsides,upuncert,by = 'idlumped_n_rep')
rm(upuncert)

# Save intermediate file (so don't have to calculate fmeys again)
# write_csv(upsides,"Data/TBD_IF_NEEDED/bycatch-upuncert-input.csv")

# load 'upuncert' if starting new session
# upsides <- read_csv("Data/TBD_IF_NEEDED/bycatch-upuncert-input.csv")

# Clean column names so that they match
upsides <- upsides %>%
  rename(idoriglumped = idoriglumped.x,
         repnumber = repnumber.x,
         commname = commname.x) %>%
  select(idorig,idoriglumped,commname,sciname,
         country,speciescat,speciescatname,regionfao,
         k,k2,fvfmsy,g,g2,beta,phi,price,eqfvfmey,fmeyvfmsy,
         curr_f,f_mey,repnumber)
upsides <- upsides %>%
  mutate(gdiff = g - g2)
test <- sum(upsides$gdiff,na.rm = TRUE)/5400000
rm(test)
upsides <- upsides %>%
  select(-g2)

# Adjust ks (from lumped to unlumped)

# temporary code (skip if running from the beginning. this was just added due to an omission in the original code)
upsides2 <- read_csv("Data/upsides.csv", col_types = cols(regionfao = "c")) %>%
  select(idorig,k) %>%
  rename(k2 = k)
upsides <- left_join(upsides,upsides2,by = 'idorig')
rm(upsides2)
# end temporary

upidslumped <- upidslumped %>%
  rename(idoriglumped = idorig,
         klumped = k) %>%
  select(idoriglumped,klumped)
upsides <- left_join(upsides,upidslumped,by = 'idoriglumped')
rm(upidslumped)
upsides <- upsides %>%
  mutate(k_adj_factor = k2/klumped) %>%
  mutate(k_adj = k * k_adj_factor) %>%
  select(-k) %>%
  rename(k = k_adj) %>%
  select(idorig,idoriglumped,commname,sciname,
         country,speciescat,speciescatname,regionfao,
         k,fvfmsy,g,beta,phi,price,eqfvfmey,fmeyvfmsy,
         curr_f,f_mey,repnumber)

# Add totoaba
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
                     repnumber = 0 # 0 denotes the fact that there is only one replication of totoaba
) 
upsides <- bind_rows(upsides,totoab)
rm(totoab)

# Add pctred columns
upsides <- upsides %>%
  mutate(pctredfmsy = 100 * (1-(g/curr_f)),
         pctredfmey = 100 * (1-(f_mey/curr_f)))
upsides <- upsides %>%
  mutate(id_n_rep = paste(idorig,repnumber,sep = "-")) # we add this column so that there is a column in which each row is unique, 
                                                       #  for the marginal cost calculation below

# Save intermediate file (so don't have to calculate fmeys again)
# write.csv(upsides,"bycatch-upuncert-input.csv",row.names = F)

# load 'upsides' if starting new session
# upsides <- read_csv("bycatch-upuncert-input.csv")

## Recalculate marginal costs for unlumped stocks
# calculation function
source("R/bycatch_funcs.R")

margcost <- function(fmey,price,g,k,phi,beta) {
  mc <- uniroot((function (x) mprofitf(fmey,price,x,g,k,phi,beta)), 
                lower = -1000000000, upper = 100000000000)[1]$root
}

# parallelized calculation
upsides2 <- upsides %>%
  filter(f_mey > 0)
num_cores <- min(24, detectCores()) # The specified number or the number of CPUs your computer has, whichever is smaller.
dtup <- 
  pbmclapply(
    1:length(upsides2$id_n_rep),
    mc.cores = num_cores,
    function(i){
      
      idsup <- upsides2$id_n_rep[i]
      
      mcup <- 
        margcost(upsides2$f_mey[i],upsides2$price[i],
                 upsides2$g[i],upsides2$k[i],
                 upsides2$phi[i],upsides2$beta[i])
      
      df = data_frame(id_n_rep = idsup, marginalcost = mcup)
      return(df)
    }) %>%
  bind_rows %>%
  filter(marginalcost > 0)

upsides <- left_join(upsides, dtup, by = 'id_n_rep')
rm(upsides2)
## end marginal cost calculation 

# save marginal costs
# write_csv(dtup,"Data/TBD_IF_NEEDED/uncert-marg-costs-dtup.csv")

# load marginal costs if starting new session
# dtup <- read_csv("uncert-marg-costs-dtup.csv")

## clean up
rm(dtup)
## end clean up

## Add a database ID (i.e. if stock assessment source was FAO or RAM legacy)
dbaseid <- read_csv("Data/TBD_IF_NEEDED/ram_stock_lkup.csv")
names(dbaseid) %<>% tolower
upsides <-
  upsides %>%
  left_join(dbaseid %>% select(-lumpedid)) %>%
  select(dbase, everything()) %>%
  mutate(dbase = ifelse(is.na(dbase), "FAO", dbase)) %>%
  mutate(dbase = ifelse(idorig=="toto", "Totoaba", dbase)) ## Totoaba was added manually (not part of upsides)

# Export final file
write_csv(upsides, "Data/upsides_uncert.csv")



####################################################
############ No uncertainty data setup #############
####################################################

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
  left_join(idlookup)
rm(upsmallu)

upsides <- left_join(upsides, upsides_kobe, by = 'idoriglumped') %>%
  mutate(curr_f = g * fvfmsy,
         f_mey = g * eqfmeyvfmsy,
         pctredfmsy = 100 * (1 - (1/fvfmsy)),
         pctredfmey = 100 * (1 - (1/eqfvfmey))) %>%
  rename(fmeyvfmsy = eqfmeyvfmsy) %>%
  select(idorig,idoriglumped,commname,sciname,country,speciescat,speciescatname,regionfao,
         k,fvfmsy,g,beta,phi,price,eqfvfmey,fmeyvfmsy,curr_f,f_mey,pctredfmsy,pctredfmey)


# OPTIONAL: change Fs and Bs to 3-year geometric mean (2010-2012) in upsides stocks
upsidesunlumped <- read_csv("Data/TBD_IF_NEEDED/UnlumpedProjectionData.csv", col_types = cols(RegionFAO = "c")) # UnlumpedProjectionData is by stock and country
names(upsidesunlumped) %<>% tolower
up2012 <- 
  upsidesunlumped %>% 
  filter(year %in% c(2012)) %>%
  select(idorig,fvfmsy) %>%
  rename(fvfmsy2012 = fvfmsy)
up2011 <- 
  upsidesunlumped %>% 
  filter(year %in% c(2011)) %>%
  select(idorig,fvfmsy) %>%
  rename(fvfmsy2011 = fvfmsy)
up2010 <- 
  upsidesunlumped %>% 
  filter(year %in% c(2010)) %>%
  select(idorig,fvfmsy) %>%
  rename(fvfmsy2010 = fvfmsy)
up3yr <- up2012 %>%
  left_join(up2011, by = 'idorig') %>%
  left_join(up2010, by = 'idorig') %>%
  mutate(fvfmsy3yr = (fvfmsy2010 * fvfmsy2011 * fvfmsy2012)^(1/3)) %>%
  select(idorig,fvfmsy3yr)
upsides <- 
  upsides %>%
  left_join(up3yr, by = 'idorig') %>%
  mutate(fvfmsy = fvfmsy3yr,
         curr_f = fvfmsy * g,
         eqfvfmey = curr_f/f_mey) %>%
  select(-fvfmsy3yr)
rm(upsidesunlumped,up2010,up2011,up2012,up3yr)
# END OPTIONAL

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
# dbaseid <- read_csv("Data/TBD_IF_NEEDED/ram_stock_lkup.csv") ## already loaded
# names(dbaseid) %<>% tolower
upsides <-
  upsides %>%
  left_join(dbaseid %>% select(-lumpedid)) %>%
  select(dbase, everything()) %>%
  mutate(dbase = ifelse(is.na(dbase), "FAO", dbase)) %>%
  mutate(dbase = ifelse(idorig=="toto", "Totoaba", dbase)) ## Totoaba was added manually (not part of upsides)

# rm(dbaseid,idlookup)

# Write final input file
write_csv(upsides, "Data/upsides_nouncert.csv")


############ Diagnostics #############
max((upsides %>% filter(fmeyvfmsy > 0))$fmeyvfmsy) ## no fmeyvfmsy
hist((upsides %>%
        filter(fmeyvfmsy > 0))$fmeyvfmsy,
     breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,30000),
     xlim = c(0,2))
######## end diagnostics
