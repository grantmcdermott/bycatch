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
    mutate(pctred_b = wgt * 100 * (1 - (f_mp/curr_f))) %>%
    
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
    mutate(pctred_b = wgt * 100 * (1 - (f_my/curr_f))) %>%
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
  imp <- inverse(function (mp) redncost_giv_mp(dt, mp)$pctred, 0.0001, mxmp)
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
      pctredb < meanpctredmey,
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


############################################################################
############################################################################
### Function 1: Obtain distribution of sampled reductions in MSY and MEY ###
############################################################################
############################################################################

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
    stocks_df <- 
      if(dt$type == 1){
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter(speciescat %in% dt$spcat) 
      } else if (dt$type == 2) {
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
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
      }
    
    stocks_df <-
      stocks_df %>%
      filter(fmeyvfmsy > 0,
             marginalcost > 0) %>%
      mutate(wgt = marginalcost * ((curr_f)^beta)) %>%
      select(idorig,pctredfmsy,pctredfmey,wgt,k,fvfmsy,g,beta,phi,price,marginalcost,eqfvfmey,curr_f,f_mey) %>%
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

# input is list of 'dt's' -> output from 'extract_func'
single_worldstate_outputs <- 
  function(dt2, n2, pctredb, reltdf) {
    samp <- lapply(dt2,
                   function (x) samp_func(x, n2, reltdf)) %>%
      bind_rows() %>%
      mutate(wgt = wt * (1/n2),
             pctrmsywt = wgt * pctredfmsy,
             pctrmeywt = wgt * pctredfmey)
    
    mpctmsy <- sum(samp$pctrmsywt)
    mpctmey <- sum(samp$pctrmeywt)
    
    stwld <- data_frame(pctredmsy = mpctmsy, 
                        pctredmey = mpctmey,
                        ycostmsy = cost_yield(samp, pctredb, mpctmsy),
                        pcostmey = cost_profit(samp, pctredb, mpctmey)) # need pctredpt from somewhere
    return(stwld)
  }

# input is list of 'dt's'
disb_func <-
  function(dt2, n1 = 10000, n2 = 100){
    
    ##########################################################
    ## Step 3.1: Create dataframe of relevant target stocks ##
    ##########################################################
    rel_targets <- lapply(dt2,
                          upsides_subset_func) %>%
      bind_rows()
    
    pctredbt <- (bycatch_df %>% 
                   filter(species == rel_targets$bycsp[1]))$pctredbpt[1]

    #########################################################################################################
    ## Step 3.2: Repeatedly sample from stocks data frame to create distribution of both pctreds and costs ##
    #########################################################################################################

    dists <-
      pblapply(1:n1, function(i) {
        single_worldstate_outputs(dt2, n2, pctredbt, rel_targets)
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

##################################################################
### Plot that asks: Will rebuilding be enough to stop decline? ###
##################################################################

bycatchdistggplot <-
  function(bdist){
    
    df <- left_join(bdist, bycatch_df) %>% group_by(species) %>%
      select(-ycostmsy, -pcostmey)
    
    df %>%
      rename(MSY = pctredmsy,
             MEY = pctredmey) %>%
      select(MSY:species) %>%
      gather(key, pctred, -species) %>%
      ggplot() +
      geom_rect(data = filter(df, row_number() == 1),
                aes(ymin = -Inf, ymax = Inf, xmin = pctredbl, xmax = pctredbu),
                alpha = .25) +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = pctred, col = key, fill = key), alpha = .5) +
      geom_vline(data = filter(df, row_number() == 1),
                 aes(xintercept = pctredbpt), lty = 2) +
      labs(x = "Reduction in mortality (%)", y = "Density") +
      # xlim(0, 100) +
      #scale_color_brewer(name = "", palette = "Set1") +
      #scale_fill_brewer(name = "", palette = "Set1") +
      scale_color_manual(values=c("red", "blue")) +
      scale_fill_manual(values=c("red", "blue")) +
      facet_wrap(~species) +
      theme(
        #text = element_text(family = font_type),
        legend.position = "bottom",
        legend.title = element_blank(),
        # strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "white"), ## Facet strip
        panel.margin = unit(2, "lines") ## Increase gap between facet panels
      ) 
  } 

##############################################################
### Plot that asks: How much will it cost to stop decline? ###
##############################################################
costggplot <-
  function(bdist){
    
    df <- left_join(bdist, bycatch_df) %>% group_by(species) %>%
      select(-pctredmsy, -pctredmey)
    
    df %>%
      rename(MSY = ycostmsy,
             MEY = pcostmey) %>%
      select(MSY:species) %>%
      gather(key, pctred, -species) %>%
      ggplot() +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = pctred, col = key, fill = key), alpha = .5) +
      labs(x = "Cost (% of MSY or MEY)", y = "Density") +
      # xlim(0, 100) +
      #scale_color_brewer(name = "", palette = "Set1") +
      #scale_fill_brewer(name = "", palette = "Set1") +
      scale_color_manual(values=c("red", "blue")) +
      scale_fill_manual(values=c("red", "blue")) +
      facet_wrap(~species) +
      theme(
        #text = element_text(family = font_type),
        legend.position = "bottom",
        legend.title = element_blank(),
        # strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "white"), ## Facet strip
        panel.margin = unit(2, "lines") ## Increase gap between facet panels
      ) 
  } 

