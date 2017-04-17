require(readr)
require(dplyr)
require(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(pbapply)
library(magrittr)
require(stringr)
library(pbmcapply)

#################################################
############ Read in/clean the data #############
#################################################

## Load Dan's Boostraps for lumped, non-NEI stocks
upuncert <- read_csv("upsidessamples.csv")
names(upuncert) %<>% tolower

# Add in price and marginal cost columns (which are not bootstrapped)
load("ProjectionData.rdata")
rm(OriginalProjectionData)
upsmalll <- as_data_frame(ProjectionData) %>% 
  select(Year,IdOrig,CommName,SciName,Country,g,phi,k,MSY,MarginalCost,SpeciesCat,SpeciesCatName,RegionFAO,Policy,beta,FvFmsy,BvBmsy,Catch,Price) %>%
  filter(Year %in% c(2012,2050))
rm(ProjectionData,UnlumpedProjectionData)
upidslumped <- select(upsmalll,Year,IdOrig,CommName,MarginalCost,Price) %>%
  filter(Year %in% c(2012))
names(upidslumped) %<>% tolower
rm(upsmalll)
upuncert <- left_join(upuncert,upidslumped,by='idorig')

# Add beta column and rename fvfmsy column
upuncerttest <- upuncert %>%
  mutate(beta = 1.3) %>%
  rename(fvfmsy = finalfvfmsy,
         bvbmsy = finalbvbmsy)
upuncert <- upuncerttest
rm(upuncerttest)

# Save intermediate file (so don't have to calculate fmeys again)
write.csv(upuncert,"bycatch-upuncert-input.csv",row.names = F)

# load 'test' if starting new session
upuncert <- read_csv("bycatch-upuncert-input.csv")

# Calculate fmey's for each bootstrap
num_cores <- min(4, detectCores()) # The specified number or the number of CPUs your computer has, whichever is smaller.
test <- 
  pbmclapply(
    1:length(upuncert$idorig),
    mc.cores = num_cores,
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
                 upper = 1.5)$par * 
             (1/upuncert$fvfmsy[i]))
      
      df = data_frame(idorig = ids, eqfvfmey = eqfvfmeys)
      return(df)
    }) %>%
  bind_rows

# Save intermediate file (so don't have to calculate fmeys again)
# write.csv(test,"fmeys-uncert.csv",row.names = F)

# load 'test' if starting new session
test <- read_csv("fmeys-uncert.csv")

# Merge fmeys into original table
upuncert$eqfvfmey <- test$eqfvfmey
rm(test)

upuncert <- upuncert %>%
  rename(fmeyvfmsy = eqfvfmey) %>%
  mutate(curr_f = g * fvfmsy,
         f_mey = g * fmeyvfmsy,
         eqfvfmey = (fvfmsy * g)/f_mey) %>%
  mutate(pctredfmsy = 100 * (1 - (1/fvfmsy)),
         pctredfmey = 100 * (1 - (1/eqfvfmey)))

# Save intermediate file (so don't have to calculate fmeys again)
# write.csv(upuncert,"bycatch-upuncert-input.csv",row.names = F)

# load 'upuncert' if starting new session
upuncert <- read_csv("bycatch-upuncert-input.csv")

upuncert <- upuncert %>%
  rename(idoriglumped = idorig)

# Create unlumped version of table
# create repeated rows in the country table

# load unlumped table
upsides <- read_csv("Data/upsides.csv", col_types = cols(regionfao = "c")) %>%
  select(idorig,idoriglumped,commname,sciname,country,speciescat,speciescatname,regionfao,g) %>%
  rename(g2 = g)

test <- upidslumped %>%
  ungroup() %>%
  mutate(nei = grepl(' nei',commname)) %>%
  filter(nei == FALSE) %>%
  rename(idoriglumped = idorig) %>%
  select(idoriglumped, nei) 

upsides <- left_join(upsides, test, by = 'idoriglumped')
upsides <- upsides %>%
  filter(nei == FALSE)

# diagnostics
median((upuncert %>%
          filter(fmeyvfmsy > 0))$g)
length((upuncert %>%
          filter(fmeyvfmsy > 0))$g)
# end diagnostics



  
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

#################################################
############ Write final input file #############
#################################################
write.csv(upuncert,"bycatch-upuncert-input.csv",row.names = F)
