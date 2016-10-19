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
      mprofitf(0,price,marginalcost,g,k,phi,beta) < mp,
      output <- 0,
      output <- imp(mp)
    ) 
    return(output)
  }

inv_marg_yield <- 
  function(my,g,k,phi) {
    imy <- 
      inverse(function (f) myieldf(f,g,k,phi), 0, g)
    ifelse(
      myieldf(0,g,k,phi) < my,
      output <- 0,
      output <- imy(my)
    ) 
    return(output)
  }

#1.3 Function that takes the sample dataframe of target stocks, and a marginal yield or profit, 
#      and computes the pct reduction for the bycatch species, and the total cost 
#      (as a percentage of total MSY or MEY), given that marginal cost.

redncost_giv_mp <- function(df, mp) {
  df$gen_id <- 1:nrow(df)
  ids <- c()
  fmp <- c()
  for (i in 1:length(df$idorig)) {
    ids <- append(ids, df$gen_id[i])
    fmp <- append(fmp, inv_marg_profit(mp,df$f_mey[i],df$price[i],
                                df$marginalcost[i],df$g[i],
                                df$k[i],df$phi[i],df$beta[i]))
  }
  dt <- data_frame(gen_id = ids, f_mp = fmp)
  dt <- left_join(df, dt, by = 'gen_id')
  
  dt <- dt %>%
    mutate(pctred_b = wgt * 100 * (1 - (f_mp/curr_f))) %>%
    mutate(mey = eqprofitf(f_mey,price,marginalcost,g,k,phi,beta)) %>%
    mutate(pcost = profitcostf(f_mp,price,marginalcost,f_mey,g,k,phi,beta))
  
  output <- data_frame(pctred = 0, cost = 0)
  output$pctred <- sum(dt$pctred_b)
  output$cost <- 100 * (sum(dt$pcost)/sum(dt$mey))
  return(output)
}

redncost_giv_my <- function(df, my) {
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

##2. Given dataframe and pctredpt (the bycatch reduction target, called 'pctredb' in the function), 
#     calculate marginal costs (y and p).

my_calc <- function(df, pctredb) {
  dt <- df %>%
    mutate(maxmy = myieldf(0,g,k,phi))
  mxmy <- max(dt$maxmy)
  imy <- inverse(function (my) redncost_giv_my(dt, my)$pctred, 0, mxmy)
  ifelse(
    redncost_giv_my(dt, mxmy)$pctred < pctredb,
    output <- mxmy,
    output <- imy(pctredb)
  ) 
  return(output)
}

mp_calc <- function(df, pctredb) {
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
      
      #1. Calculate marginal yield cost, given pct reduction needed for bycatch species (my_calc function)
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
      
      #1. Calculate marginal profit cost, given pct reduction needed for bycatch species (mp_calc function)
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

#######################################################################
## Step 2: Produce a weighted data frame of the sample distributions ##
#######################################################################

disb_func <-
  function(dt, n1 = 10000, n2 = 100){
    
    ##############################################################
    ## Step 2.1: Select the target stocks per relevant criteria ##
    ##############################################################
    
    stocks_df <- 
      if(dt$type == 1){
        upsides %>%
          filter(regionfao %in% dt$faoreg) %>%
          filter(speciescat %in% dt$spcat) 
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
      }
    
    stocks_df <-
      stocks_df %>%
      filter(fmeyvfmsy > 0) %>%
      mutate(wgt = marginalcost * ((curr_f)^beta)) %>%
      select(idorig,pctredfmsy,pctredfmey,wgt,fvfmsy,g,beta,phi,price,marginalcost,eqfvfmey,curr_f,f_mey) %>%
      mutate(wgt = wgt/sum(wgt, na.rm = T))
    
    #########################################################################################################
    ## Step 2.2: Repeatedly sample from stocks data frame to create distribution of both pctreds and costs ##
    #########################################################################################################
    
    samp <-
      pblapply(1:n1, function(i) {
        stocks_df %>%
          sample_n(n2, replace = T, weight = wgt) %>%
          summarise(pctredmsy = mean(pctredfmsy),
                    pctredmey = mean(pctredfmey))
      }) %>%
      bind_rows()
    
    ##################################################################
    ## Step 2.3: Multiply sampled data frame by target stock weight ##
    ##################################################################
    
    samp_wt <- dt$wt * samp
    
    return(samp_wt)
  }

###############################################################################
## Final step: Convenience wrapper to apply over all target stocks affecting ##
## a single bycatch species. This is what we'll call in the actual analysis. ##
###############################################################################

bycatch_func <- 
  function(z){
    lapply(
      extract_func(z),
      function(x) disb_func(x, n1, n2)
    ) %>%
      Reduce("+", .) %>%
      as_data_frame() %>%
      mutate(species = z)
  }

## E.g. bycatch_func("Loggerhead_turtle")

##############################################################################
### Function 2: Plot that asks: Will rebuilding be enough to stop decline? ###
##############################################################################

bycatchdistggplot <-
  function(bdist){
    
    df <- left_join(bdist, bycatch_df) %>% group_by(species)
    
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
      scale_color_brewer(name = "", palette = "Set1") +
      scale_fill_brewer(name = "", palette = "Set1") +
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

