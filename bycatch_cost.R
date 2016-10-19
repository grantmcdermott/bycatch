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
  select(-marginalcost) #drop original marginalcost row 

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

## Load bycatch data
bycatch_df <- read_csv("Data/bycatch_species.csv")

target_df <- read_csv("Data/target_species.csv")


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


###########################
### Test cost functions ###
###########################
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
    filter(fmeyvfmsy > 0,
           marginalcost > 0) %>%
    mutate(wgt = marginalcost * ((curr_f)^beta)) %>%
    select(idorig,pctredfmsy,pctredfmey,wgt,fvfmsy,g,beta,phi,price,marginalcost,eqfvfmey,curr_f,f_mey,k) %>%
    mutate(wgt = wgt/sum(wgt, na.rm = T))
  return(dtuse)
}

# Test target stocks
testdf <- stockselect1(upsides, c(31,33,34), c(21,31))

# Test target sample
testsamp <- testdf %>%
  sample_n(100, replace = T, weight = wgt) %>%
  mutate(wgt = wgt/sum(wgt, na.rm = T)) %>%
  mutate(pctredwt = pctredfmsy * wgt) %>%
  mutate(pctredwtp = pctredfmey * wgt)
testsamp$gen_id <- 1:nrow(testsamp)

fmytest <- c()
idstest <- c()
for (i in 1:length(testsamp$idorig)) {
  idstest <- append(idstest, testsamp$gen_id[i])
  fmytest <- append(fmytest, inv_marg_yield(20,testsamp$g[i],testsamp$k[i],testsamp$phi[i]))
}
dttest <- data_frame(gen_id = idstest, f_my = fmytest)
dttest <- left_join(testsamp, dttest, by = 'gen_id')

# debug mp_calc
pctredbtest <- 70
dttest <- testsamp %>%
    mutate(maxmp = mprofitf(0,price,marginalcost,g,k,phi,beta))
mxmptest <- max(dttest$maxmp)
imptest <- inverse(function (mp) redncost_giv_mp(dttest, mp)$pctred, 0.0001, mxmptest)
ifelse(
    redncost_giv_mp(dttest, mxmptest)$pctred < pctredbtest,
    outputtest <- mxmptest,
    outputtest <- imptest(pctredbtest)
  ) 

source("bycatch_funcs_cost_2.R")

meanpct <- sum(testsamp$pctredwt)
meanpctp <- sum(testsamp$pctredwtp)

cost_yield(testsamp, 80, meanpct)

cost_profit(testsamp, 99.99999, meanpctp)

utest <- upsides %>%
  mutate(maxmc = mprofitf(f_mey,price,0,g,k,phi,beta)) %>%
  filter(maxmc < 0)
rm(utest)

# Test component functions
inv_marg_profit(100000,testsamp$f_mey[1],testsamp$price[1],testsamp$marginalcost[1],testsamp$g[1],testsamp$k[1],testsamp$phi[1],testsamp$beta[1])

mprofitf(0.1,testsamp$price[1],testsamp$marginalcost[1],testsamp$g[1],testsamp$k[1],testsamp$phi[1],testsamp$beta[1])

myieldf(0,testsamp$g[1],testsamp$k[1],testsamp$phi[1])

f_imyl <- inv_marg_yield(20,testsamp$g[1],testsamp$k[1],testsamp$phi[1])

redncost_giv_mp(testsamp, 100000) 

redncost_giv_my(testsamp, 20)

my_calc(testsamp, 60)

mp_calc(testsamp, 70)
