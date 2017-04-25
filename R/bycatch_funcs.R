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
  
  dt <- dt %>%
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
  
  dt <- dt %>%
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
  dt <- df %>%
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
    pctredb > 100,
    ycost <- 100,
    ifelse(
      # Check if there is a shortfall.
      # If not, cost = 0 (because MSY reduces bycatch enough to stop decline)
      pctredb < meanpctredmsy,
      ycost <- 0,
      
      # There is a shortfall, and it is attainable through additional target species F reductions.
      # Code below calculates costs. It is assumed that df includes weights of each target stock (which sum to 1)
      
      #1. Calculate marginal yield cost (with function 'my_calc'), 
      #     given pct reduction needed for bycatch species.
      #2. Calculate total yield cost, given marginal cost calculated in 1. 
      ycost <- redncost_giv_my(df, my_calc(df, pctredb))$cost
    )
  )
  return(ycost)
}

cost_profit <- function(df, pctredb, meanpctredmey) {
  ifelse(
    # Check if stopping all fishing is enough for bycatch target.
    # If not, cost = 100 (percent of fishing profits)
    pctredb > 100,
    pcost <- 100,
    ifelse(
      # Check if there is a shortfall.
      # If not, cost = 0 (because MEY reduces bycatch enough to stop decline)
      pctredb < meanpctredmey + 1.2,
      pcost <- 0,
      
      # There is a shortfall, and it is attainable through additional target species F reductions.
      # Code below calculates costs. It is assumed that df includes weights of each target stock (which sum to 1)
      
      #1. Calculate marginal profit cost (with function 'mp_calc'), 
      #     given pct reduction needed for bycatch species.
      #2. Calculate total profit cost, given marginal cost calculated in 1. 
      pcost <- redncost_giv_mp(df, mp_calc(df, pctredb))$cost
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
      filter(fmeyvfmsy > 0,
             marginalcost > 0) %>%
      mutate(wgt = marginalcost * ((curr_f)^beta)) %>%
      select(idorig,idoriglumped,pctredfmsy,pctredfmey,wgt,speciescat,speciescatname,fmeyvfmsy, ## GRM: ADDED speciescat/name,fmeyvfmsy
             k,fvfmsy,g,beta,phi,price,marginalcost,eqfvfmey,curr_f,f_mey) %>%
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
  function(dt2, n2, pctredb, pctredbl, pctredbu, reltdf) {
    
    samp <- 
      lapply(dt2, function (x) samp_func(x, n2, reltdf)) %>%
      bind_rows() %>%
      mutate(
        wgt = wt * (1/n2),
        pctrmsywt = wgt * 100 * (1 - ((1 - (pctredfmsy/100))^alpha_exp)),
        pctrmeywt = wgt * 100 * (1 - ((1 - (pctredfmey/100))^alpha_exp))
        )
    
    mpctmsy <- sum(samp$pctrmsywt)
    mpctmey <- sum(samp$pctrmeywt)
    
    pctb <- runif(1, min = pctredbl, max = pctredbu)
    
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
  function(dt2, n1 = 10000, n2 = 100){
    
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
      pblapply(1:n1, function(i) {
        single_worldstate_outputs(dt2, n2, pctredbt, pctredbtl, pctredbtu, rel_targets)
      }) %>%
      bind_rows()
    
    return(dists)
  }

###############################################################################
## Final step: Convenience wrapper to apply over all target stocks affecting ##
## a single bycatch species. This is what we'll call in the actual analysis. ##
###############################################################################

bycatch_func <- 
  function(z){
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


##################################################################
### Plot that asks: Will rebuilding be enough to stop decline? ###
##################################################################

bycatchdistggplot <-
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
costggplot <-
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
      filter(n()==2) ## Make sure we only keep cases that can be depicted in both facets
    df$key <- factor(df$key, levels = c(lvl0_a, lvl0_b), labels = c(lvl1_a, lvl1_b))
    
    excld_species <- anti_join(distinct(summ_df, species), distinct(df, species))$species
    
    print(noquote(paste0("Note: The following (outlier) species has been excluded from the plot: ", excld_species)))
    
    # pd = position_dodge(0.5)
    
    df %>%
      ggplot(aes(x = pctredbpt, y = q50, col = grp, fill = grp, group = species)) +
      # geom_abline(
      #   data = data.frame(key=levels(df$key), b = c(1, NA)),
      #   aes(intercept=0, slope = b), lty = 2
      #   ) +
      geom_segment(
        data = data.frame(x1=0, x2=1, y1=0, y2=1, key=lvl1_a),
        inherit.aes = F,
        aes(x = x1, y = y1, xend = x2, yend = y2), col="black", lty=2
        ) +
      geom_errorbarh(aes(xmin = pctredbl_dash, xmax = pctredbu_dash, col = grp), height = 0, linetype = "21", alpha = 0.7) +
      geom_errorbarh(aes(xmin = pctredbl, xmax = pctredbu, col = grp), height = 0, alpha = 0.7) +
      geom_point(stroke = 0.25, size = 3, alpha = 0.7) +
      geom_point(stroke = 0.25, size = 3, shape = 1) +
      geom_errorbar(aes(ymin = q025, ymax = q975, col = grp), width = 0, alpha = 0.7) +
      # geom_pointrange(aes(ymin = q025, ymax = q975), fatten = 6, alpha = 0.7, stroke = 0.25, shape = 21) +
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
        labels=percent,
        oob = rescale_none
        ) +
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
#   geom_point(aes(col = grp, fill = grp), alpha = 0.7, stroke = 0.25, size = 3) +
#   geom_point(aes(col = grp, fill = grp), shape = 1, stroke = 0.25, size = 3) +
#   # scale_color_brewer(palette = "Set1") +
#   scale_colour_manual(values = bycatch_cols) +
#   geom_errorbar(aes(ymax = q975, ymin = q025, col = grp), width = 0) +
#   geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl, col = grp), height = 0) +
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


### OLD FIGURE 3 ETC CODE
# pctmsyplot <-
#   function(bdist){
#     
#     df <- left_join(bdist, bycatch_df, by = 'species') %>% group_by(species) 
#     
#     df %>%
#       rename(MSY = ycostmsy,
#              MEY = pcostmey) %>%
#       select(MSY:species) %>%
#       gather(key, pctcost, -species) %>%
#       ggplot() +
#       # geom_line(stat = "density") + ## lines only
#       geom_errorbar(limits, width=0.2) +
#       geom_density(aes(x = pctcost, y = ..scaled.., col = key, fill = key), alpha = .5) +
#       labs(x = "Cost (% of MSY or MEY)", y = "Density") +
#       # xlim(0, 100) +
#       #scale_color_brewer(name = "", palette = "Set1") +
#       #scale_fill_brewer(name = "", palette = "Set1") +
#       scale_color_manual(values=c("red", "blue")) +
#       scale_fill_manual(values=c("red", "blue")) +
#       facet_wrap(~species) +
#       theme(
#         #text = element_text(family = font_type),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         # strip.text = element_text(size = 18, colour = "black"),
#         strip.background = element_rect(fill = "white"), ## Facet strip
#         panel.spacing = unit(2, "lines") ## Increase gap between facet panels
#       ) 
#   } 
# 
# # Multiple plot function (for Fig 3 and SI figure that shows all of sensitivty analyses)
# #
# # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# # - cols:   Number of columns in layout
# # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# #
# # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# # then plot 1 will go in the upper left, 2 will go in the upper right, and
# # 3 will go all the way across the bottom.
# #
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   library(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }
# 
# 
# 
# theme_fig3a <-
#   theme(
#     text = element_text(family = font_type),
#     legend.position = "top",
#     legend.title = element_blank(),
#     axis.line.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     strip.text = element_text(size = 18, colour = "black"),
#     strip.background = element_rect(fill = "white"), ## Facet strip
#     panel.spacing = unit(2, "lines") ## Increase gap between facet panels
#   )
# 
# theme_fig3b <-
#   theme(
#     text = element_text(family = font_type),
#     legend.position = "none",
#     legend.title = element_blank(),
#     # strip.text = element_text(size = 18, colour = "black"),
#     strip.background = element_rect(fill = "white"), ## Facet strip
#     panel.spacing = unit(2, "lines") ## Increase gap between facet panels
#   )
# 
# fig3mey <-
#   function(resultssum, yminimum){
#     
#     df <- resultssum 
#     
#     f3a <- (ggplot(data = df, aes(x = pctredbpt, y = pctredmey50, col = grp, fill = grp)) +
#               geom_abline(slope = 1, lty = 2) + 
#               geom_point(alpha = 0.7, stroke = 0.25, size = 4) + 
#               geom_point(shape = 1, stroke = 0.25, size = 4) +
#               scale_color_brewer(palette = "Set1") +
#               geom_errorbar(aes(ymax = pctredmey975, ymin = pctredmey025), width = 0) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), height = 0) +
#               #xlim(0,100) +
#               #ylim(-50,100) +
#               scale_x_continuous(expand = c(0, 0), limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(expand = c(0.05, 0), limits=c(yminimum,100), oob = rescale_none) +
#               # geom_hline(yintercept = 0, lty = 2) + 
#               # geom_hline(yintercept = 51, lty = 2, colour = "grey33") + ## mean reduction across our population samples
#               labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
#                    y = "Projected % reduction \nin mortality (MEY)") +
#               theme_fig3a
#     )
#     
#     f3b <- (ggplot(data = df, aes(x = pctredbpt, y = pcostmey50, col = grp, fill = grp)) +
#               geom_point(alpha = 0.7, stroke = 0.25, size = 4) + 
#               geom_point(shape = 1, stroke = 0.25, size = 4) +
#               scale_color_brewer(palette = "Set1") +
#               geom_errorbar(aes(ymax = pcostmey975, ymin = pcostmey025), width = 0) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), height = 0) +
#               #xlim(0,100) +
#               #ylim(0,100) +
#               scale_x_continuous(expand = c(0, 0), limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(expand = c(0.05, 0), limits=c(0,100), oob = rescale_none) +
#               # geom_vline(xintercept = 51, lty = 2, colour = "grey33") + ## mean reduction across our population samples
#               labs(x = "Required % reduction in \nmortality to halt decline", 
#                    y = "Projected profit cost \n (% of MEY)") +
#               theme_fig3b
#     )
#     # figure3 <- multiplot(f3a, f3b)
#     figure3 <- plot_grid(f3a, f3b, labels = c("A", "B"), align = "v")
#     
#     return(figure3)
#     
#   } 
# 
# fig3msy <-
#   function(resultssum, yminimum){
#     
#     df <- resultssum 
#     
#     f3a <- (ggplot(data = df, aes(pctredbpt, pctredmsy50, col = grp, fill = grp)) +
#               geom_abline(slope = 1, lty = 2) + 
#               geom_point(alpha = 0.7, stroke = 0.25, size = 4) + 
#               geom_point(shape = 1, stroke = 0.25, size = 4) +
#               scale_color_brewer(palette = "Set1") +
#               geom_errorbar(aes(ymax = pctredmsy975, ymin = pctredmsy025), width = 0) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), height = 0) +
#               #xlim(0,100) +
#               #ylim(-50,100) +
#               scale_x_continuous(expand = c(0, 0), limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
#               # geom_hline(yintercept = 0, lty = 2) + 
#               # geom_hline(yintercept = 23, lty = 2, colour = "grey33") + 
#               labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
#                    y = "Projected % reduction \nin mortality (MSY)") +
#               theme_fig3a
#     )
#     
#     f3b <- (ggplot(data = df, aes(pctredbpt, ycostmsy50, col = grp, fill = grp)) +
#               geom_point(alpha = 0.7, stroke = 0.25, size = 4) + 
#               geom_point(shape = 1, stroke = 0.25, size = 4) +
#               scale_color_brewer(palette = "Set1") +
#               geom_errorbar(aes(ymax = ycostmsy975, ymin = ycostmsy025), width = 0) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), height = 0) +
#               #xlim(0,100) +
#               #ylim(0,100) +
#               scale_x_continuous(expand = c(0, 0), limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(limits=c(0,100), oob = rescale_none) +
#               labs(x = "Required % reduction in \nmortality to halt decline", 
#                    y = "Projected profit cost \n (% of MSY)") +
#               theme_fig3b
#     )
#     figure3 <- multiplot(f3a, f3b)
#     return(figure3)
#     
#   } 
# 
# # Combined fig (Alt. Fig 3 if decided to show MEY *and* MSY)
# fig3all <-
#   function(resultssum, yminimum){
#     
#     df <- resultssum 
#     
#     fa <- (ggplot(data = df, aes(pctredbpt, pctredmey50)) +
#               geom_point(aes(colour = grp), size = 4) + 
#               scale_color_brewer(palette = "Paired") +
#               geom_errorbar(aes(ymax = pctredmey975, ymin = pctredmey025), width = 0.25) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
#               #xlim(0,100) +
#               #ylim(-50,100) +
#               scale_x_continuous(limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
#               geom_abline(slope = 1) + 
#               # geom_hline(yintercept = 0, lty = 2) + 
#               geom_hline(yintercept = 51, lty = 2, colour = "grey33") + 
#               labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
#                    y = "% reduction in mortality \n projected (MEY)") +
#              theme_fig3a 
#     )
#     
#     fb <- (ggplot(data = df, aes(pctredbpt, pcostmey50)) +
#               geom_point(aes(colour = grp), size = 4) + 
#               scale_color_brewer(palette = "Paired") +
#               geom_errorbar(aes(ymax = pcostmey975, ymin = pcostmey025), width = 0.25) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
#               #xlim(0,100) +
#               #ylim(0,100) +
#               scale_x_continuous(limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(limits=c(0,100), oob = rescale_none) +
#               geom_vline(xintercept = 51, lty = 2, colour = "grey33") + 
#               labs(x = "% reduction in mortality \n needed to halt decline", 
#                    y = "Projected profit cost \n (% of MEY)") +
#              theme_fig3b
#     )
#     
#     fc <- (ggplot(data = df, aes(pctredbpt, pctredmsy50)) +
#               geom_point(aes(colour = grp), size = 4) + 
#               scale_color_brewer(palette = "Paired") +
#               geom_errorbar(aes(ymax = pctredmsy975, ymin = pctredmsy025), width = 0.25) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
#               #xlim(0,100) +
#               #ylim(-50,100) +
#               scale_x_continuous(limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
#               geom_abline(slope = 1) + 
#               # geom_hline(yintercept = 0, lty = 2) + 
#               geom_hline(yintercept = 23, lty = 2, colour = "grey33") + 
#               labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
#                    y = "% reduction in mortality \n projected (MSY)") +
#              theme_fig3a 
#     )
#     
#     fd <- (ggplot(data = df, aes(pctredbpt, ycostmsy50)) +
#               geom_point(aes(colour = grp), size = 4) + 
#               scale_color_brewer(palette = "Paired") +
#               geom_errorbar(aes(ymax = ycostmsy975, ymin = ycostmsy025), width = 0.25) +
#               geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
#               #xlim(0,100) +
#               #ylim(0,100) +
#               scale_x_continuous(limits=c(0,100), oob = rescale_none) +
#               scale_y_continuous(limits=c(0,100), oob = rescale_none) +
#               geom_vline(xintercept = 23, lty = 2, colour = "grey33") + 
#               labs(x = "% reduction in mortality \n needed to halt decline", 
#                    y = "Projected yield cost \n (% of MSY)") +
#              theme_fig3b
#     )
#     figure3 <- multiplot(fa, fb, fc, fd, cols = 2)
#     return(figure3)
#     
#   } 


#########################
### Sensitivity plots ###
#########################

## Alt. Fig S5 with MSY

sensfigmsy <-
  function(resultssum_nouncert, resultssum_uncert, resultssum_alpha05, resultssum_alpha2, yminimum){
    
    df_nouncert <- resultssum_nouncert 
    df_uncert <- resultssum_uncert 
    df_alpha05 <- resultssum_alpha05 
    df_alpha2 <- resultssum_alpha2 
    
    fa <- (ggplot(data = df_nouncert, aes(pctredbpt, pctredmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = pctredmsy975, ymin = pctredmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(-50,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
              geom_abline(slope = 1) + 
              # geom_hline(yintercept = 0, lty = 2) + 
              geom_hline(yintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                   y = "% reduction in mortality \n projected (MEY)") +
             theme_fig3a 
    )
    
    fb <- (ggplot(data = df_nouncert, aes(pctredbpt, ycostmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = ycostmsy975, ymin = ycostmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(0,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(0,100), oob = rescale_none) +
              geom_vline(xintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = "% reduction in mortality \n needed to halt decline", 
                   y = "Projected yield cost \n (% of MSY)") +
             theme_fig3b
    )
    
    fc <- (ggplot(data = df_uncert, aes(pctredbpt, pctredmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = pctredmsy975, ymin = pctredmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(-50,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
              geom_abline(slope = 1) + 
              # geom_hline(yintercept = 0, lty = 2) + 
              geom_hline(yintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                   y = "% reduction in mortality \n projected (MEY)") +
             theme_fig3a 
    )
    
    fd <- (ggplot(data = df_uncert, aes(pctredbpt, ycostmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = ycostmsy975, ymin = ycostmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(0,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(0,100), oob = rescale_none) +
              geom_vline(xintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = "% reduction in mortality \n needed to halt decline", 
                   y = "Projected yield cost \n (% of MSY)") +
             theme_fig3b
    )
    
    fe <- (ggplot(data = df_alpha05, aes(pctredbpt, pctredmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = pctredmsy975, ymin = pctredmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(-50,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
              geom_abline(slope = 1) + 
              # geom_hline(yintercept = 0, lty = 2) + 
              geom_hline(yintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                   y = "% reduction in mortality \n projected (MEY)") +
             theme_fig3a 
    )
    
    ff <- (ggplot(data = df_alpha05, aes(pctredbpt, ycostmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = ycostmsy975, ymin = ycostmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(0,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(0,100), oob = rescale_none) +
              geom_vline(xintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = "% reduction in mortality \n needed to halt decline", 
                   y = "Projected yield cost \n (% of MSY)") +
             theme_fig3b
    )
    
    fg <- (ggplot(data = df_alpha2, aes(pctredbpt, pctredmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = pctredmsy975, ymin = pctredmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(-50,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
              geom_abline(slope = 1) + 
              # geom_hline(yintercept = 0, lty = 2) + 
              geom_hline(yintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                   y = "% reduction in mortality \n projected (MEY)") +
             theme_fig3a 
    )
    
    fh <- (ggplot(data = df_alpha2, aes(pctredbpt, ycostmsy50)) +
              geom_point(aes(colour = grp), size = 4) + 
              scale_color_brewer(palette = "Paired") +
              geom_errorbar(aes(ymax = ycostmsy975, ymin = ycostmsy025), width = 0.25) +
              geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
              #xlim(0,100) +
              #ylim(0,100) +
              scale_x_continuous(limits=c(0,100), oob = rescale_none) +
              scale_y_continuous(limits=c(0,100), oob = rescale_none) +
              geom_vline(xintercept = 23, lty = 2, colour = "grey33") + 
              labs(x = "% reduction in mortality \n needed to halt decline", 
                   y = "Projected yield cost \n (% of MSY)") +
             theme_fig3b
    )
    
    figure3 <- multiplot(fa, fb, fc, fd, fe, ff, fg, fh, cols = 4)
    return(figure3)
    
  } 

## Fig. S5

sensfigmey <-
  function(resultssum_nouncert, resultssum_uncert, resultssum_alpha05, resultssum_alpha2, yminimum){
    
    df_nouncert <- resultssum_nouncert 
    df_uncert <- resultssum_uncert 
    df_alpha05 <- resultssum_alpha05 
    df_alpha2 <- resultssum_alpha2 
    
    fa <- (ggplot(data = df_nouncert, aes(pctredbpt, pctredmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pctredmey975, ymin = pctredmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(-50,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
             geom_abline(slope = 1) + 
             # geom_hline(yintercept = 0, lty = 2) + 
             geom_hline(yintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                  y = "% reduction in mortality \n projected (MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "top",
               legend.title = element_blank(),
               axis.line.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    fb <- (ggplot(data = df_nouncert, aes(pctredbpt, pcostmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pcostmey975, ymin = pcostmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(0,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(0,100), oob = rescale_none) +
             geom_vline(xintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = "% reduction in mortality \n needed to halt decline", 
                  y = "Projected profit cost \n (% of MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "none",
               legend.title = element_blank(),
               # strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    fc <- (ggplot(data = df_uncert, aes(pctredbpt, pctredmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pctredmey975, ymin = pctredmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(-50,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
             geom_abline(slope = 1) + 
             # geom_hline(yintercept = 0, lty = 2) + 
             geom_hline(yintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                  y = "% reduction in mortality \n projected (MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "top",
               legend.title = element_blank(),
               axis.line.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    fd <- (ggplot(data = df_uncert, aes(pctredbpt, pcostmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pcostmey975, ymin = pcostmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(0,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(0,100), oob = rescale_none) +
             geom_vline(xintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = "% reduction in mortality \n needed to halt decline", 
                  y = "Projected profit cost \n (% of MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "none",
               legend.title = element_blank(),
               # strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    fe <- (ggplot(data = df_alpha05, aes(pctredbpt, pctredmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pctredmey975, ymin = pctredmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(-50,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
             geom_abline(slope = 1) + 
             # geom_hline(yintercept = 0, lty = 2) + 
             geom_hline(yintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                  y = "% reduction in mortality \n projected (MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "top",
               legend.title = element_blank(),
               axis.line.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    ff <- (ggplot(data = df_alpha05, aes(pctredbpt, pcostmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pcostmey975, ymin = pcostmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(0,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(0,100), oob = rescale_none) +
             geom_vline(xintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = "% reduction in mortality \n needed to halt decline", 
                  y = "Projected profit cost \n (% of MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "none",
               legend.title = element_blank(),
               # strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    fg <- (ggplot(data = df_alpha2, aes(pctredbpt, pctredmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pctredmey975, ymin = pctredmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(-50,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(yminimum,100), oob = rescale_none) +
             geom_abline(slope = 1) + 
             # geom_hline(yintercept = 0, lty = 2) + 
             geom_hline(yintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = NULL, #"% reduction in mortality \n needed to halt decline", 
                  y = "% reduction in mortality \n projected (MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "top",
               legend.title = element_blank(),
               axis.line.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    fh <- (ggplot(data = df_alpha2, aes(pctredbpt, pcostmey50)) +
             geom_point(aes(colour = grp), size = 4) + 
             scale_color_brewer(palette = "Paired") +
             geom_errorbar(aes(ymax = pcostmey975, ymin = pcostmey025), width = 0.25) +
             geom_errorbarh(aes(xmax = pctredbu, xmin = pctredbl), width = 0.25) +
             #xlim(0,100) +
             #ylim(0,100) +
             scale_x_continuous(limits=c(0,100), oob = rescale_none) +
             scale_y_continuous(limits=c(0,100), oob = rescale_none) +
             geom_vline(xintercept = 51, lty = 2, colour = "grey33") + 
             labs(x = "% reduction in mortality \n needed to halt decline", 
                  y = "Projected profit cost \n (% of MEY)") +
             theme(
               #text = element_text(family = font_type),
               legend.position = "none",
               legend.title = element_blank(),
               # strip.text = element_text(size = 18, colour = "black"),
               strip.background = element_rect(fill = "white"), ## Facet strip
               panel.spacing = unit(2, "lines") ## Increase gap between facet panels
             )
    )
    
    figure3 <- multiplot(fa, fb, fc, fd, fe, ff, fg, fh, cols = 4)
    return(figure3)
    
  } 

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
        q025 = quantile(value, .025),
        q50 = quantile(value, .5),
        q975 = quantile(value, .975)
      ) %>%
      left_join(bycatch_df, by = 'species') %>%
      select(species, grp, region, contains("pctred"), everything())
  }

# resultssummary <- function(df) {
#   dt <- df %>%
#     group_by(species) %>%
#     summarise(pctredmsy025 = quantile(pctredmsy, probs = 0.025),
#               pctredmsy25 = quantile(pctredmsy, probs = 0.25),
#               pctredmsy50 = median(pctredmsy),
#               pctredmsy75 = quantile(pctredmsy, probs = 0.75),
#               pctredmsy975 = quantile(pctredmsy, probs = 0.975),
#               pctredmey025 = quantile(pctredmey, probs = 0.025),
#               pctredmey25 = quantile(pctredmey, probs = 0.25),
#               pctredmey50 = median(pctredmey),
#               pctredmey75 = quantile(pctredmey, probs = 0.75),
#               pctredmey975 = quantile(pctredmey, probs = 0.975),
#               ycostmsy025 = quantile(ycostmsy, probs = 0.025),
#               ycostmsy25 = quantile(ycostmsy, probs = 0.25),
#               ycostmsy50 = median(ycostmsy),
#               ycostmsy75 = quantile(ycostmsy, probs = 0.75),
#               ycostmsy975 = quantile(ycostmsy, probs = 0.975),
#               pcostmey025 = quantile(pcostmey, probs = 0.025),
#               pcostmey25 = quantile(pcostmey, probs = 0.25),
#               pcostmey50 = median(pcostmey),
#               pcostmey75 = quantile(pcostmey, probs = 0.75),
#               pcostmey975 = quantile(pcostmey, probs = 0.975))
#   rs <- left_join(dt, bycatch_df, by = 'species')
#   return(rs)
# }
