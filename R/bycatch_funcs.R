##################################################################
##################################################################
### Equilibrium profit and yield functions (for cost analysis) ###
##################################################################
##################################################################

# Using literal interpretations of eqn 8 in Costello et al. SI
eqprofitf <- 
  function(f,price,marginalcost,g,k,phi,beta) {
    ymsy <- ((g * k)/((phi + 1)^(1/phi)))
    bvbm <- ((1 + phi - ((f/g) * phi))^(1/phi))
    eqp <- 
      (price * ymsy * (f/g) * bvbm) -
      (marginalcost * ((f)^beta))
    return(eqp)
  }

# In terms of f:
eqyieldf <- 
  function(f,g,k,phi) {
    eqy <- 
      f *
      (k * ((1 - ((f * phi)/(g * (phi + 1))))^(1/phi)))
    return(eqy)
  }

# Old version
# eqprofitf <- 
#   function(f,price,marginalcost,g,k,phi,beta) {
#     eqp <-
#       (price * f *
#          eqyieldf(f,g,k,phi)) - 
#       (marginalcost * ((f)^beta))
#     return(eqp)
#   }


# Cost of f, relative to fmey (profit) or msy (yield):
profitcostf <- 
  function(f,price,marginalcost,fmey,g,k,phi,beta) {
    pc <-
      eqprofitf(fmey,price,marginalcost,g,k,phi,beta) -
      eqprofitf(f,price,marginalcost,g,k,phi,beta)
    return(pc)
  }

yieldcostf <- 
  function(f,g,k,phi) {
    yc <- 
      eqyieldf(g,g,k,phi) -
      eqyieldf(f,g,k,phi)
    return(yc)
  }

# Marginal profit/yield, in terms of f:
myieldf <- 
  function(f,g,k,phi) {
    my <- 
      (((g - f) * k *
          ((1 - ((f * phi)/(g * (phi + 1))))^((1/phi) - 1)))
       /g)
    return(my)
  }

# Using literal interpretations of eqn 8 in Costello et al. SI
mprofitf <- 
  function(f,price,marginalcost,g,k,phi,beta) {
    mp <-
      ((price * ((1 + phi - ((f * phi)/g))^(1/phi)) * (f - g) * k * ((1 + phi)^((phi - 1)/phi)))/
         ((f * phi) - (g * (1 + phi)))) -
      (beta * marginalcost * (f^(beta - 1)))
    return(mp)
  }

# Old version
# mprofitf <- 
#   function(f,price,marginalcost,g,k,phi,beta) {
#     mp <-
#       (-beta * (f^(beta - 1)) * marginalcost) +
#       (price * myieldf(f,g,k,phi))
#     return(mp)
#   }


#######################################
#######################################
### Chris' analytical cost analysis ###
#######################################
#######################################

##1. Given dataframe and marginal cost, return total cost

#1.1 Define numerical inverse function
inverse = function (f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]$root
}
# inverse <- possibly(inverse, NA_real)

#1.2 Inverse marginal cost (in eq. yield or profit) of reducing f for one target stock:
#     Input is marginal cost of reducing f (denoted 'mp', equal to marginal profit from increasing f ('my' for yield))
#     Output is f that gives you that marginal cost (i.e. the inverse of the eq. marginal profit function)
#     If f cannot be found, 0 is returned (because it implies that fishing must have shut down on stock of interest)

inv_marg_profit <- 
  function(mp,fmey,price,marginalcost,g,k,phi,beta) {
    imp <-
      inverse(function (f) mprofitf(f,price,marginalcost,g,k,phi,beta), 0, fmey)
    ifelse(
      # check for corner
      mprofitf(0,price,marginalcost,g,k,phi,beta) < mp,
      
      # corner
      output <- 0,
      
      # interior
      output <- imp(mp)
    ) 
    return(output)
  }

inv_marg_yield <- 
  function(my,g,k,phi) {
    imy <- 
      inverse(function (f) myieldf(f,g,k,phi), 0, g)
    ifelse(
      # check for corner
      myieldf(0,g,k,phi) < my,
      
      # corner
      output <- 0,
      
      # interior
      output <- imy(my)
    ) 
    return(output)
  }

#1.3 Function that takes the sample dataframe of target stocks, and a marginal yield or profit, 
#      and computes the pct reduction for the bycatch species, and the total cost 
#      (as a percentage of total MSY or MEY), given that marginal cost.

redncost_giv_mp <- function(df, mp) {
  # Add row ids to the input dataframe ('df')
  df$gen_id <- 1:nrow(df)
  
  # Create empty lists
  ids <- c()
  fmp <- c()
  for (i in 1:length(df$idorig)) {
    # Generate ids matching df row ids created above. 
    # These will be used to  merge the new column 'fmp' created below into df.
    ids <- append(ids, df$gen_id[i])
    # Calculate fishing mortality, fmp, for each stock (i) in the sample, df, 
    # given the marginal profit, mp, given as input.
    fmp <- append(fmp, inv_marg_profit(mp,df$f_mey[i],df$price[i],
                                df$marginalcost[i],df$g[i],
                                df$k[i],df$phi[i],df$beta[i]))
  }
  
  # Merge the fmp's calculated above into the dataframe df, 
  # in new column called 'f_mp'.
  dt <- data_frame(gen_id = ids, f_mp = fmp)
  dt <- left_join(df, dt, by = 'gen_id')
  
  dt <- 
    dt %>%
    # Create new column 'pctred_b' which specifies the % reduction in fishing
    # mortality for the bycatch stock resulting from the % reduction
    # for each target stock (relative to 2012) at f_mp.
    mutate(pctred_b = wgt * 100 * (1 - ((f_mp/curr_f)^alpha_exp))) %>%
    # Create new column 'mey' which specifies the mey for
    # each target stock.
    mutate(mey = eqprofitf(f_mey,price,marginalcost,g,k,phi,beta)) %>%
    # Create new column 'pcost' which specifies the mey for
    # each target stock.
    mutate(pcost = profitcostf(f_mp,price,marginalcost,f_mey,g,k,phi,beta))
  
  # Create empty data frame to contain the results.
  output <- data_frame(pctred = 0, cost = 0)
  
  # % reduction in mortality for bycatch species (pctred), given marginal profit (mp),
  # is the sum of pctred_b across all target stocks in df.
  output$pctred <- sum(dt$pctred_b)
  
  # 'cost' is the difference between the total profit at the f's having marginal cost, mp,
  # (and therefore required for % bycatch reduction 'pctred') and mey, as a fraction of mey,
  # cumulatively across all target stocks in df.
  output$cost <- 100 * (sum(dt$pcost)/sum(dt$mey))
  return(output)
}

redncost_giv_my <- function(df, my) { # function is analogous to 'redncost_giv_mp' above
  df$gen_id <- 1:nrow(df)
  ids <- c()
  fmy <- c()
  for (i in 1:length(df$idorig)) {
    ids <- append(ids, df$gen_id[i])
    fmy <- append(fmy, inv_marg_yield(my,df$g[i],df$k[i],df$phi[i]))
  }
  dt <- data_frame(gen_id = ids, f_my = fmy)
  dt <- left_join(df, dt, by = 'gen_id')
  
  dt <- 
    dt %>%
    mutate(pctred_b = wgt * 100 * (1 - ((f_my/curr_f)^alpha_exp))) %>%
    mutate(msy = eqyieldf(g,g,k,phi)) %>%
    mutate(ycost = yieldcostf(f_my,g,k,phi))
  
  output <- data_frame(pctred = 0, cost = 0)
  output$pctred <- sum(dt$pctred_b)
  output$cost <- 100 * (sum(dt$ycost)/sum(dt$msy))
  return(output)
}

##2. Given dataframe ('df') and pctredpt (the bycatch reduction target, called 'pctredb' in the function), 
#     calculate marginal costs (y and p).

my_calc <- function(df, pctredb) {
  dt <- 
    df %>%
    # Generate new column in df of maximum marginal yield (i.e. my at 0 fishing)
    mutate(maxmy = myieldf(0,g,k,phi))
  
  # store largest 'maxmy' among stocks in df as 'mxmy'
  mxmy <- max(dt$maxmy)
  
  # Create 'imy' function which calculates the marginal yield, given pctredb
  imy <- inverse(function (my) redncost_giv_my(dt, my)$pctred, 0, mxmy)
  ifelse(
    # Check if the maximum marginal yield ('mxmy') corresponds 
    # to a still-insufficient reduction in bycatch.
    redncost_giv_my(dt, mxmy)$pctred < pctredb,
    # If yes, then require the max marginal yield.
    output <- mxmy,
    # If no, use 'imy' function to calculate required marginal yield
    output <- imy(pctredb)
  ) 
  return(output)
}

mp_calc <- function(df, pctredb) { # function is analogous to 'my_calc' above
  dt <- df %>%
    mutate(maxmp = mprofitf(0,price,marginalcost,g,k,phi,beta))
  mxmp <- max(dt$maxmp)
  imp <- inverse(function (mp) redncost_giv_mp(dt, mp)$pctred, 1, mxmp)
  ifelse(
    redncost_giv_mp(dt, mxmp)$pctred < pctredb,
    output <- mxmp,
    output <- imp(pctredb)
  ) 
  return(output)
}

##3 Main function: Takes as input a dataframe of relevant stocks 
#                  and relevent columns from these, as well as total reduction needed
#                  for bycatch species, and percent reduction from MSY/MEY
#                  Returns 0 if there is no shortfall. 
#                  Returns 100 if pctredpt > 100 (i.e. if stopping all fishing isn't enough to stop decline). 
#                  Returns yield/profit cost as fraction of total MEY/MSY profit/yield otherwise.

cost_yield <- function(df, pctredb, meanpctredmsy) {
  
  ifelse(
    # Check if stopping all fishing is enough for bycatch target.
    # If not, cost = 100 (percent of fishing yield)
    round(pctredb) >= 100, ## GRM: Changed from `pctredb > 100
    ycost <- 100,
    ifelse(
      # Check if there is a shortfall.
      # If not, cost = 0 (because MSY reduces bycatch enough to stop decline)
      round(pctredb) < round(meanpctredmsy), ## GRM: Added rounding
      ycost <- 0,
      
      # There is a shortfall, and it is attainable through additional target species F reductions.
      # Code below calculates costs. It is assumed that df includes weights of each target stock (which sum to 1)
      
      #1. Calculate marginal yield cost (with function 'my_calc'), 
      #     given pct reduction needed for bycatch species.
      #2. Calculate total yield cost, given marginal cost calculated in 1. 
      # ycost <- redncost_giv_my(df, my_calc(df, pctredb))$cost
      ycost <-
        tryCatch(
          withTimeout(redncost_giv_my(df, my_calc(df, pctredb))$cost, timeout=5), 
          TimeoutException=function(ex) NA
          )
    )
  )
  return(ycost)
}

cost_profit <- function(df, pctredb, meanpctredmey) {
  ifelse(
    # Check if stopping all fishing is enough for bycatch target.
    # If not, cost = 100 (percent of fishing profits)
    round(pctredb) >= 100, ## GRM: Changed from `pctredb > 100`
    pcost <- 100,
    ifelse(
      # Check if there is a shortfall.
      # If not, cost = 0 (because MEY reduces bycatch enough to stop decline)
      round(pctredb) < (meanpctredmey), ## GRM: Changed from `pctredb < meanpctredmey + 1.2` (Rounds to the nearest whole percent)
      pcost <- 0,
      # There is a shortfall, and it is attainable through additional target species F reductions.
      # Code below calculates costs. It is assumed that df includes weights of each target stock (which sum to 1)
      #1. Calculate marginal profit cost (with function 'mp_calc'), 
      #     given pct reduction needed for bycatch species.
      #2. Calculate total profit cost, given marginal cost calculated in 1. 
      # pcost <- redncost_giv_mp(df, mp_calc(df, pctredb))$cost
      pcost <-
        tryCatch(
          withTimeout(redncost_giv_mp(df, mp_calc(df, pctredb))$cost, timeout=5), 
          TimeoutException=function(ex) NA
          )
    )
  )
  return(pcost)
}


##################################################################################
##################################################################################
### Bycatch function: Obtain distribution of sampled reductions in MSY and MEY ###
##################################################################################
##################################################################################

################################################################################################
## Step 1: Extract relevant list from target species data frame (needed for disb_func below) ###
################################################################################################

extract_func <-
  function(s){
    lapply(target_df %>% 
             filter(species == s) %>%
             split(.$target), 
           function(y){
             lapply(y, function(x) unique(x))
           }
    ) 
  }

###########################################################
## Step 2: Subset upsides data to relevant target stocks ##
###########################################################

upsides_subset_func <- 
  function(dt) {
    
    ## GRM: NEW TESTING
    upsides <-
      upsides %>% 
      separate_rows(regionfao)
    
    stocks_df <- 
      if(dt$type == 1){
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          # filter(grepl(paste0("\\b",dt$faoreg,"\\b"), regionfao)) %>% ## GRM CHANGED: To handle cases with a multi-FAO region string
          filter(speciescat %in% dt$spcat)
          # filter(grepl(paste0("\\b",dt$spcat,"\\b"), speciescat)) ## GRM CHANGED: To handle cases with a multi-species cat
      } else if (dt$type == 2) {
        upsides %>%
          filter(speciescat %in% dt$spcat) %>%
          filter(country %in% dt$countries)  
      } else if (dt$type == 3) {
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>% 
          filter(idoriglumped %in% dt$lumpedstocks)  
      } else if (dt$type == 4) {
        upsides %>%
          filter(idoriglumped %in% dt$lumpedstocks) %>%
          filter(country %in% dt$countries)   
      } else if (dt$type == 5) {
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>% 
          filter(speciescat %in% dt$spcat) %>%
          filter(country %in% dt$countries)  
      } else if (dt$type == 6) { # Leatherback demersal - remove U.S. stocks
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter(speciescat %in% dt$spcat) %>%
          filter(country != "USA")
      } else if (dt$type == 7) { # Indo-Pacific finless porpoise - add countries
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter(speciescat %in% dt$spcat) %>%
          filter(country %in% c("China",
                                "Taiwan Province of China",
                                "China Hong Kong SAR",
                                "China Macao SAR",
                                "Viet Nam",
                                "Cambodia",
                                "Thailand",
                                "Malaysia",
                                "Indonesia",
                                "Singapore",
                                "Brunei",
                                "Philippines",
                                "Myanmar",
                                "Bangladesh",
                                "India",
                                "Sri Lanka",
                                "Pakistan",
                                "Saudi Arabia",
                                "Oman",
                                "Iran"))
      } else if (dt$type == 8) { # Indo-Pacific humpack dolphin - remove stocks from Russia, Japan, Koreas
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter(speciescat %in% dt$spcat) %>%
          filter(country != "Russian Federation") %>%
          filter(country != "Japan") %>%
          filter(country != "Republic of Korea") %>%
          filter(country != "Democratic People's Republic of Korea")
      } else if (dt$type == 9) { # W Pacific leatherback - add countries
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter(speciescat %in% dt$spcat) %>%
          filter(country %in% c("Malaysia",
                                "Indonesia",
                                "Papua New Guinea",
                                "Solomon Islands"))
      } else if (dt$type == 10) { # NW Indian Ocean loggerhead turtle - add countries
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter(speciescat %in% dt$spcat) %>%
          filter(country %in% c("India",
                                "Sri Lanka",
                                "Pakistan",
                                "Saudi Arabia",
                                "Oman",
                                "Iran",
                                "Iraq",
                                "Kuwait",
                                "Bahrain",
                                "Qatar",
                                "United Arab Emirates",
                                "Yemen",
                                "Egypt",
                                "Sudan",
                                "Eritrea",
                                "Djibouti",
                                "Somalia"))
      } else if (dt$type == 11) { # Stocks for white-chinned petrel
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter((speciescat == '36')|
                   (sciname == "Dissostichus eleginoides")|
                   (grepl("Merluccius", sciname) == TRUE)
                 )    # Dissostichus eleginoides, Merluccius
      } else if (dt$type == 12) { # Hector's dolphin and Maui's dolphin
        upsides %>%
          filter(regionfao == 81) %>%
          filter(country == "New Zealand") %>%
          filter(sciname %in% c("Notodarus sloanii", # List from Table 1 in Hickford et al. 1997
                                "Thyrsites atun",
                                "Parapercis colias",
                                "Paristiopterus labiosus",
                                "Caesioperca lepidoplera",
                                "Pseudophycis breviuscula",
                                "Colistium guntheri",
                                "Scymnorhinus licha",
                                "Odax pullus",
                                "Cephaloscyllium isabellum",
                                "Latridopsis forsteri",
                                "Conger verreauxi",
                                "Myliobatis tenuicaudatus",
                                "Callorhinchus milii",
                                "Scomber australasicus",
                                "Peltorhamphus novaezeelandiae",
                                "Lepidopus caudatus",
                                "Aplodactylus arctidens",
                                "Chelidonichthys kumu",
                                "Polyprion oxygeneios",
                                "Macruronus novaezelandiae",
                                "Zeus faber",
                                "Trachurus declivis",
                                "Trachurus novaezelandiae",
                                "Arripis trutta",
                                "Seriola lalandi",
                                "Parika scaber",
                                "Cenypterus blacodes",
                                "Pelotretis flavilatus",
                                "Scorpis violaceus",
                                "Latridopsis ciliaris",
                                "Octopus maorum",
                                "Macrouridae",
                                "Plagiogeneion nibiginosus",
                                "Pseudophycis bachus",
                                "Paratrachichthys trailli",
                                "Cheilodactylus spectabilis",
                                "Upeneichthys lineatus",
                                "Scorpaena cardinalis",
                                "Raja nasuta",
                                "Galeorhinus galeus",
                                "Rhombosolea plebeia",
                                "Rexea solandri",
                                "Chrysophrys auratus",
                                "Squalus acanthias",
                                "Helicolenus percoides",
                                "Pseudolabrus miles",
                                "Mustelus lenticulatus",
                                "Genyagnus monopterygius",
                                "Raja innominata",
                                "Scorpis lineolatus",
                                "Nemadactylus macropterus",
                                "Pseudocaranx dentex",
                                "Seriolella brama",
                                "Arnoglossus scapha"))
      }
    
    stocks_df <-
      stocks_df %>%
      filter(
        curr_f > 0,
        fmeyvfmsy > 0,
        marginalcost > 0,
        fconmsy > 0
        ) %>%
      mutate(wgt = marginalcost * ((curr_f)^beta)) %>%
      select(
        dbase,idorig,idoriglumped,pctredfmsy,pctredfmey,pctredfmey,pctredfmsycon,pctredfmeycon,
        wgt,speciescat,speciescatname,fmeyvfmsy,k,fvfmsy,g,beta,phi,price,marginalcost,eqfvfmey,
        curr_f,f_mey,fconmsy,fconmey
        ) %>%
      mutate(trgcat = dt$target) %>%
      mutate(wt = dt$wt) %>%
      mutate(bycsp = dt$species) %>%
      mutate(wgt = wgt/sum(wgt, na.rm = T))
    
    return(stocks_df)
  }

########################################################################
## Step 3: Produce a data frame of the sample distributions and costs ##
########################################################################

samp_func <- 
  function(dt, n2, reltdf) {
    smple <- reltdf %>%
      filter(trgcat == dt$target) %>%
      sample_n(n2, replace = T, weight = wgt)
    return(smple)
  }

## Input is a list of 'dt's' -> output from 'extract_func'
single_worldstate_outputs <- 
  function(dt2, n2, pctredb, pctredbl, pctredbu, reltdf, sensrangept25) {
   
    if (scenario == "All stocks") {
    samp <- 
      lapply(dt2, function (x) samp_func(x, n2, reltdf)) %>%
      bind_rows() %>%
      mutate(
        wgt = wt * (1/n2),
        pctrmsywt = wgt * 100 * (1 - ((1 - (pctredfmsy/100))^alpha_exp)),
        pctrmeywt = wgt * 100 * (1 - ((1 - (pctredfmey/100))^alpha_exp))
        )
    } else {
      samp <- 
        lapply(dt2, function (x) samp_func(x, n2, reltdf)) %>%
        bind_rows() %>%
        mutate(
          wgt = wt * (1/n2),
          pctrmsywt = wgt * 100 * (1 - ((1 - (pctredfmsycon/100))^alpha_exp)),
          pctrmeywt = wgt * 100 * (1 - ((1 - (pctredfmeycon/100))^alpha_exp))
        )
    }
    
    mpctmsy <- sum(samp$pctrmsywt)
    mpctmey <- sum(samp$pctrmeywt)
    
    if (sensrangept25 == 1) {
    pctb <- runif(1, min = pctredbl, max = pctredbu) # if 25% uncertainty in Fe and delta is on
    } else {
      pctb <- pctredb
    }
    stwld <- 
      data_frame(
        pctredmsy = mpctmsy, 
        pctredmey = mpctmey,
        ycostmsy = cost_yield(samp, pctb, mpctmsy),
        pcostmey = cost_profit(samp, pctb, mpctmey)
        ) # need pctredpt from somewhere
    return(stwld)
  }

## Input is list of 'dt's'
disb_func <-
  function(dt2, n1, n2){
    
    ##########################################################
    ## Step 3.1: Create dataframe of relevant target stocks ##
    ##########################################################
    rel_targets <- 
      lapply(dt2, upsides_subset_func) %>%
      bind_rows()
    
    pctredbt <- (filter(bycatch_df, species==rel_targets$bycsp[1]))$pctredbpt[1]
    pctredbtl <- (filter(bycatch_df, species==rel_targets$bycsp[1]))$pctredbl[1]
    pctredbtu <- (filter(bycatch_df, species==rel_targets$bycsp[1]))$pctredbu[1]

    #########################################################################################################
    ## Step 3.2: Repeatedly sample from stocks data frame to create distribution of both pctreds and costs ##
    #########################################################################################################

    dists <-
      pblapply(1:n1, 
               possibly(function(i) {
                 evalWithTimeout(
                   single_worldstate_outputs(dt2, n2, pctredbt, pctredbtl, pctredbtu, rel_targets, sensrange25), 
                   timeout = 20, ## i.e. Time out after 20 seconds if can't resolve 
                   TimeoutException = function(ex) "TimedOut"
                   )
                 }, 
                 # otherwise = NULL
                 otherwise = data_frame(pctredmsy=NA, pctredmey=NA, ycostmsy=NA, pcostmey=NA) ## To catch failed uniroot cases
                 ),
               cl=num_cores
               ) %>%
      bind_rows()
    
    return(dists)
  }

###############################################################################
## Final step: Convenience wrapper to apply over all target stocks affecting ##
## a single bycatch species. This is what we'll call in the actual analysis. ##
###############################################################################

bycatch_func <- 
  function(z){
    print(z)
    disb_func(extract_func(z), n1, n2) %>%
      as_data_frame() %>%
      mutate(species = z)
  }

## E.g. bycatch_func("Loggerhead_turtle")




######################
######################
### PLOT FUNCTIONS ###
######################
######################


### Function to extract legend (to serve as common legend at bottom of composite figures) 
g_legend <- 
  function(a_ggplot){ 
    tmp <- ggplot_gtable(ggplot_build(a_ggplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)
  } 

##################################################################
### Plot that asks: Will rebuilding be enough to stop decline? ###
##################################################################

bycatchdist_plot <-
  function(bdist){
    
    df <- 
      left_join(bdist, bycatch_df) %>% 
      group_by(species) %>%
      select(-ycostmsy, -pcostmey) %>%
      mutate_if(is.double, funs(. / 100)) 
    
    df %>%
      rename(MSY = pctredmsy,
             MEY = pctredmey) %>%
      select(MSY:species) %>%
      gather(key, pctred, -species) %>%
      mutate(key = factor(key, levels = c("MSY", "MEY"))) %>% 
      ggplot() +
      geom_rect(data = filter(df, row_number() == 1),
                aes(ymin = -Inf, ymax = Inf, xmin = pctredbl, xmax = pctredbu),
                alpha = .25) +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = pctred, y = ..scaled.., col = key, fill = key), alpha = .5) +
      geom_vline(data = filter(df, row_number() == 1),
                 aes(xintercept = pctredbpt), lty = 2) +
      labs(x = "Reduction in mortality", y = "Density") + 
      scale_x_continuous(labels = percent) + 
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      facet_wrap(~species) +
      theme(legend.position = "bottom") 
  } 

##############################################################
### Plot that asks: How much will it cost to stop decline? ###
##############################################################
cost_plot <-
  function(bdist){
    
    df <- 
      left_join(bdist, bycatch_df) %>% 
      group_by(species) %>%
      select(-pctredmsy, -pctredmey) %>%
      mutate_if(is.double, funs(. / 100)) 
    
    df %>%
      rename(MSY = ycostmsy,MEY = pcostmey) %>%
      select(MSY:species) %>%
      gather(key, pctcost, -species) %>%
      mutate(key = factor(key, levels = c("MSY", "MEY"))) %>% 
      ggplot() +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = pctcost, y = ..scaled.., col = key, fill = key), alpha = .5, adjust = 0.5) +
      labs(x = "Cost (fraction of MSY or MEY)", y = "Density") + 
      # xlim(0, 100) +
      scale_x_continuous(limits=c(0,1), oob = rescale_none, labels = percent) + 
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      facet_wrap(~species) +
      theme(legend.position = "bottom") 
  } 


#################################################################################
### Fig. 2 Logic of analysis, as illustrated by NW Atlantic loggerhead turtle ###
#################################################################################

### Panel 2.c
## First an intermediate function that selects the relevant target stocks
stockselect_func <-
  function(speciesname){
    lapply(extract_func(speciesname), upsides_subset_func) %>%
      bind_rows() %>%
      mutate(wt_usd = marginalcost * ((g * fvfmsy)^beta)) %>%
      rename(wt_Fe = wt) %>% 
      select(idorig, speciescat, speciescatname, pctredfmsy,pctredfmey,trgcat, 
             wt_Fe,wt_usd)
  }
  
## Second, the function that actually plots the scatter figures
samples_plot <- 
  function(bdist) {

    lbl_df <- 
      data_frame(
        trgcat=c("Shrimp", "Demersals", "Tuna"),
        target_description=c("Shrimp~trawls", "Other~demersal~fisheries", "Pelagic~longline~fisheries")
      )
    
    bd <- 
      bdist %>%
      left_join(lbl_df) %>%
      rename(MSY = pctredfmsy, MEY = pctredfmey) %>%
      ## Special chars ('(', '%', etc) req. diff quotation marks and spacing tokens (~, *) for plotmath parsing to work
      mutate(target_lab = paste0("atop(", target_description, paste0(",~'('*", wt_Fe*100,"*'%'~of~italic(F)[e]*')')"))) %>%
      group_by(trgcat) %>%
      gather(key, pctred, c(MSY, MEY)) %>%
      mutate(key = factor(key, levels = c("MSY", "MEY"))) %>%
      mutate(pctred = pctred/100)
    
    bd %>% 
      ggplot(aes(x = pctred, y = wt_usd, col = key)) +
      geom_point(alpha=0.5, size = 3) +
      scale_color_brewer(palette = "Set1") +
      scale_x_continuous(limits = c(-1, 1), labels = percent) +
      scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
      labs(x = "Reduction in mortality", y = "2012 Effort (USD)") +
      facet_grid( ~ forcats::fct_reorder(target_lab, wt_Fe, .desc = T), 
                  labeller = label_parsed) +
      theme(
        # strip.text = element_text(size=14),
        legend.position = "none"
      )
}


#########################
### Fig. 3 Trade-offs ###
#########################

#########################
### Results summaries ###
#########################

summ_func <-
  function(df) {
    df %>%
      gather(key, value, -species) %>%
      group_by(species, key) %>%
      # do(q=(quantile(.$value))) %>%
      summarise(
        q025 = quantile(value, .025, na.rm = T),
        q50 = quantile(value, .5, na.rm = T),
        q975 = quantile(value, .975, na.rm = T)
      ) %>%
      left_join(bycatch_df, by = 'species') %>%
      select(species, grp, clade, region, contains("pctred"), everything())
  }


tradeoffs_plot <-
  function(summ_df, scenario) {
    
    if(toupper(scenario)=="MEY"){
      lvl0_a <- "pctredmey"
      lvl0_b <- "pcostmey"
      lvl1_a <- "Projected reduction in \nmortality (MEY scenario)"
      lvl1_b <- "Projected profit cost \n(as fraction of MEY)"
    }else{
      lvl0_a <- "pctredmsy"
      lvl0_b <- "ycostmsy"
      lvl1_a <- "Projected reduction in \nmortality (MSY scenario)"
      lvl1_b <- "Projected profit cost \n(as fraction of MSY)"
    }
    
    df <- 
      summ_df %>%
      mutate_if(is.double, funs(./100)) %>%
      filter(grepl(scenario, key, ignore.case=T)) %>%
      filter(q50>=0) %>% 
      filter(pctredbpt<=1) %>%
      mutate(pctredbu_dash = pctredbu) %>% ## This and next two lines for visualization (see geom_errorbarh)
      mutate(pctredbu = ifelse(pctredbu<=1, pctredbu, 1)) %>%
      mutate(pctredbl_dash = pctredbu) %>%
      mutate(q025_dash = q025) %>% ## This and next two lines for visualization (see geom_errorbar)
      mutate(q025 = ifelse(q025>=0, q025, 0)) %>%
      mutate(q975_dash = q025) %>%
      filter(n()==2) ## Make sure we only keep cases that can be depicted in both facets
    df$key <- factor(df$key, levels = c(lvl0_a, lvl0_b), labels = c(lvl1_a, lvl1_b))
    
    excld_species <- anti_join(distinct(summ_df, species), distinct(df, species))$species
    
    print(noquote(paste0("Note: The following (outlier) species has been excluded from the plot: ", excld_species)))

    df %>%
      mutate(clade = stringr::str_to_title(clade)) %>%
      ggplot(aes(x = pctredbpt, y = q50, col = clade, group = species)) +
      geom_segment(
        data = data.frame(x1=0, x2=1, y1=0, y2=1, key=lvl1_a),
        inherit.aes = F,
        aes(x = x1, y = y1, xend = x2, yend = y2), col="black", lty=2
        ) +
      geom_errorbarh(aes(xmin = pctredbl_dash, xmax = pctredbu_dash, col = clade), height = 0, linetype = "21", alpha = 0.7) +
      geom_errorbarh(aes(xmin = pctredbl, xmax = pctredbu, col = clade), height = 0, alpha = 0.7) +
      geom_errorbar(aes(ymin = q025_dash, ymax = q975_dash, col = clade), width = 0, linetype = "21", alpha = 0.7) +
      geom_errorbar(aes(ymin = q025, ymax = q975, col = clade), width = 0, alpha = 0.7) +
      # geom_point(stroke = 0.25, size = 3, alpha = 0.7) +
      # geom_point(stroke = 0.25, size = 3, shape = 1) +
      geom_point(aes(shape=clade), fill="white", size = 3.5, stroke = 0) + ## to "white out" the error bar at the points
      geom_point(aes(shape=clade, fill = clade), alpha=0.7, size = 3.5, stroke = 0) +
      # geom_point(aes(shape=clade), size = 3.25, stroke = 0.25) + ## uncomment if want the points to have outlines
      geom_segment(
        data = data.frame(x1=0, x2=1, y1=0, y2=1, key=lvl1_a),
        inherit.aes = F,
        aes(x = x1, y = y1, xend = x2, yend = y2), col="black", lty=2, alpha=0.2 ## adding again (w/ low alpha to give effect behind points)
        ) +
      scale_shape_manual(values = 21:24) +
      scale_colour_manual(values = bycatch_cols) +
      scale_fill_manual(values = bycatch_cols) +
      scale_x_continuous(
        #expand = c(0, 0), 
        limits=c(0,1),
        labels=percent,
        oob = rescale_none
        ) +
      scale_y_continuous(
        #limits=c(NA,1),
        limits=c(0,1),
        labels=percent,
        oob = rescale_none
        ) +
      coord_fixed() +
      labs(
        x = "Reduction in mortality \nneeded to halt decline",
        y = NULL
        ) +
      facet_wrap(~key, nrow=2, scales = "free_y", strip.position="left") +
      theme(
        axis.title.x = element_text(size=14),
        strip.text = element_text(size=14),
        strip.placement = "outside"
      )
  }

## If want to show all data (incl. start of inset plot), put the following
## code into the above function
# df %>%
#   ggplot(aes(x = pctredbpt, y = q50)) +
#   geom_segment(
#     data = data.frame(x1=0, x2=1, y1=0, y2=1, key=lvl1_a),
#     aes(x = x1, y = y1, xend = x2, yend = y2), col="black", lty=2
#   ) +
#   geom_point(aes(col = clade, fill = clade), alpha = 0.7, stroke = 0.25, size = 3) +
#   geom_point(aes(col = clade, fill = clade), shape = 1, stroke = 0.25, size = 3) +
#   # scale_color_brewer(palette = "Set1") +
#   scale_colour_manual(values = bycatch_cols) +
#   geom_errorbar(aes(ymax = q975, ymin = q025, col = clade), width = 0) +
#   geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl, col = clade), height = 0) +
#   geom_rect(
#     inherit.aes = FALSE,
#     aes(xmin=-0.01, xmax=1.01, ymin=-0.01, ymax=1.01),
#     col="red", fill=NA
#   ) +
#   scale_x_continuous(#expand = c(0, 0),
#     labels=percent,
#     oob = rescale_none) +
#   scale_y_continuous(
#     labels=percent,
#     oob = rescale_none) +
#   labs(
#     x = "Reduction in mortality \nneeded to halt decline",
#     y = NULL
#   ) +
#   facet_wrap(~key, nrow=2, scales = "free_y", strip.position="left") +
#   theme(
#     axis.title.x = element_text(size=14),
#     strip.text = element_text(size=14),
#     strip.placement = "outside"
#   )
