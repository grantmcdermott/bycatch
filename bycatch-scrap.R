###################
### Final tests ###
###################

source("bycatch_funcs_cost_2.R")

# Test extract function 
lhtest <- extract_func(c('Loggerhead turtle'))

# Test upsides_subset_func  
testdt <- lapply(lhtest,
                 upsides_subset_func) %>%
  bind_rows()

# Test pctredbpt extraction
pctredbtest <- (bycatch_df %>% 
               filter(species == testdt$bycsp[1]))$pctredbpt[1]
pctredbtest <- 70

# Test single state of the world
swtest <- single_worldstate_outputs(lhtest, 100, pctredbtest, testdt)

# Test whole function
testdt2 <- disb_func(lhtest, n1 = 500, n2 = 100)

# Test mp_calc
stest <- samp_func(lhtest$Demersals, n2, testdt)
meanpct <- mean(stest$pctredfmsy)
meanpctp <- mean(stest$pctredfmey)
mp_calc(stest, 70)

##################
### Scrap code ###
##################


#1.4 Function that takes the sample dataframe of target stocks, and a marginal yield or profit, 
#      and computes total cost (as a percentage of total MSY or MEY), given that marginal cost.

pcost_giv_mp <- function(df, mp) {
  dt <- df %>%
    mutate(f_mp = inv_marg_profit(mp,f_mey,price,marginalcost,g,k,phi,beta)) %>%
    mutate(mey = eqprofitf(f_mey,price,marginalcost,g,k,phi,beta)) %>%
    mutate(pcost = profitcostf(f_mp,price,marginalcost,f_mey,g,k,phi,beta))
  return(100 * (sum(dt$pcost)/sum(dt$mey)))
}

ycost_giv_my <- function(df, my) {
  dt <- df %>%
    mutate(f_my = inv_marg_yield(my,g,k,phi)) %>%
    mutate(msy = eqyieldf(g,g,k,phi)) %>%
    mutate(ycost = yieldcostf(f_my,g,k,phi))
  return(100 * (sum(dt$ycost)/sum(dt$msy)))
}

### Check marginal cost parameter values ###
upsides <- upsides %>%
  mutate(fbar = ((phi + 1)/phi) * (1 - ((bvbmsy^phi)/(phi + 1)))) %>%
  mutate(mu = fbar - fvfmsy) %>%
  mutate(msy = ((g * k)/((phi + 1)^(1/phi)))) %>%
  mutate(bbar = ((1 + phi - (fbar * phi))^(1/phi))) %>%
  mutate(cpar = ((price * fbar * bbar * msy)/((g * fbar)^beta)))

### End marginal cost parameter value check ###

###########################
### Test cost functions ###
###########################
source("bycatch_funcs_cost_2.R")

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

cost_yield(testsamp, 97, meanpct)

cost_profit(testsamp, 60.6, meanpctp)

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

my_calc(testsamp, 46.50297754)

mp_calc(testsamp, 60.6)

# Looking at pieces of Grant's code.
source("bycatch_funcs_cost_2.R")

test <- target_df %>% 
         filter(species == 'Loggerhead turtle') %>%
         split(.$target)
       
test2 <- bind_rows(pblapply(turtle_species_samp, extract_func))

lhtest <- bycatch_func(c('Loggerhead turtle')) 


testlist <- c()
testlist <- append(testlist,head(upsides))

rm(testlist)
testlist[1] <- head(upsides)

test <- 3
c(1:test)

testdt <- disb_func(lhtest, n1 = 1, n2 = 100)

testsw <- single_worldstate_outputs(lhtest, 100, pctredb, reltdf)
