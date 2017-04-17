##################################################################
##################################################################
### Equilibrium profit and yield functions (for cost analysis) ###
##################################################################
##################################################################

eqprofit <- 
  function(pctred,price,marginalcost,fvfmsy,g,k,phi,beta) {
    eqp <-
      (price*(1-(0.01*pctred))*fvfmsy*g*k*
         ((1 - (((1-(0.01*pctred))*fvfmsy)/(phi + 1)))^(1/phi))) - 
      (marginalcost*((fvfmsy*g*(1-(0.01*pctred)))^beta))
      return(eqp)
  }

eqyield <- 
  function(pctred,fvfmsy,g,k,phi) {
    eqy <- 
      (1-(0.01*pctred))*fvfmsy*g*k*
      ((1 - (((1-(0.01*pctred))*fvfmsy)/(phi + 1)))^(1/phi))
    return(eqy)
  }

###################################
###################################
### Test code for cost analysis ###
###################################
###################################

# Pseudocode:
# Summary: For each bycatch draw (each of n1 draws of n2 stocks from each target sp. category),
#    input is list of drawn stocks, n2, weights for different target species, and pctredpt-
#    -the percent reduction in bycatch mortality estimated to be needed for the bycatch species of interest.
#    The objective is to calculate the minimum (yield or profit) cost needed to get bycatch to this level,
#    if MSY or MEY management does not provide enough of a reduction. This is repeated for each of the n1 draws,
#    which yields the presented distribution. The pseudocode describes the cost-minimization procedure for each draw (of size n2).
# 1. For each selected stock, calculate cost (in terms of yield or profit) of equal-sized increments of 
#    bycatch mortality reduction, from fmey (for profit) or fmsy (yield) to pctred = 100 
#    (i.e. no more fishing on that target stock). Make data frame of these increments by cost.
#    Columns of new dataframe: idorig, byc_incr (i.e. additional pctred implied for bycatch),
#    tar_incr (i.e. additional pctred implied for target), cost (yield_cost for yield, profit_cost for profit)
# 2. Combine these tables for all stocks into single dataframe. (Potential issue: data table size)
# 3. Sort by cost (from smallest to largest).
# 4. Add down the rows until the necessary reduction for the bycatch species is achieved.
# 5. Add down rows to calculate cost.
# 6. Compare cost to total profit/yield at MEY/MSY.



##### Pseudocode
#1. Subset full table into small table with weights, stock ids and pctred values (stockselect functions below)
#2. For each table, create a sampled version with n rows--which samples from pctred 
#     (i.e. a vector of pctred samples of length n)--combine MSY and MEY vectors into data frame
#3. Take weighted averages of resulting vectors for different species groups

##### end pseudocode

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
      mutate(wt = marginalcost * ((g * fvfmsy)^beta)) %>%
      mutate(pctredfmsy = 100 * (1 - (1/fvfmsy))) %>%
      mutate(pctredfmey = 100 * (1 - (fmeyvfmsy/fvfmsy))) %>%
      select(idorig,pctredfmsy,pctredfmey,wt,fvfmsy,g,beta,phi,price,marginalcost) %>%
      mutate(wt = wt/sum(wt, na.rm = T))
    
    ###############################################################################
    ## Step 2.2: Repeatedly sample from stocks data frame to create distribution ##
    ###############################################################################
    
    samp <-
      pblapply(1:n1, function(i) {
        stocks_df %>%
          sample_n(n2, replace = T, weight = wt) %>%
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



#################################################################################
### Function 3: How much does it cost to reduce mortality enough? (minimized) ###
#################################################################################

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