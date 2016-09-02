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
  function(dt, #n1, 
           n2){
    
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
      select(idorig,pctredfmsy,pctredfmey,wt,fvfmsy,g,beta,phi,k,price,marginalcost) %>% ## added k
      mutate(wt = wt/sum(wt, na.rm = T))
    
    ## Profit and yield
    stocks_df <-
      stocks_df %>%
      mutate(eqm_yield = eqyield(pctredfmsy, fvfmsy, g, k, phi),
             eqm_profit = eqprofit(pctredfmey, price, marginalcost, fvfmsy, g, k, phi, beta)
             ) %>%
      ## NEW
      mutate(
        pctredfmsy_incrmt = pctredfmsy+i*incrmt/dt$wt, ## Assume we reduce FMSY by an additional (scaled) amount c.f. upsides optimum
        pctredfmey_incrmt = pctredfmey+i*incrmt/dt$wt, ## Assume we reduce FMEY by an additional (scaled) amount c.f. upsides optimum
        eqm_yield_incrmt = eqyield(pctredfmsy+i*incrmt/dt$wt, fvfmsy, g, k, phi), ## How would this affect yield?
        eqm_profit_incrmt = eqprofit(pctredfmey+i*incrmt/dt$wt, price, marginalcost, fvfmsy, g, k, phi, beta), ## How would this affect profit?
        eqm_yield_prev_incrmt = eqyield(pctredfmsy+(i-1)*incrmt/dt$wt, fvfmsy, g, k, phi), ## Yield if we went back one (scaled) increment
        eqm_profit_prev_incrmt = eqprofit(pctredfmey+(i-1)*incrmt/dt$wt, price, marginalcost, fvfmsy, g, k, phi, beta) ## Profit if we went back one (scaled) increment
        ) %>%
      mutate(
        marg_yield = eqm_yield_incrmt - eqm_yield_prev_incrmt,
        marg_profit = eqm_profit_incrmt - eqm_profit_prev_incrmt
      ) 

    
    
    ###############################################################################
    ## Step 2.2: Repeatedly sample from stocks data frame to create distribution ##
    ###############################################################################
    
    ## NEW
    set.seed(j)
    
    samp <-
      # pblapply(1:n1, function(i) {
        stocks_df %>%
          sample_n(n2, replace = T, weight = wt) %>%
          summarise(pctredmsy = mean(pctredfmsy),
                    pctredmey = mean(pctredfmey),
                    eqm_yield = mean(eqm_yield, na.rm = T),
                    eqm_profit = mean(eqm_profit, na.rm = T),
                    ## NEW
                    pctredmsy_incrmt = mean(pctredfmsy_incrmt), 
                    pctredmey_incrmt = mean(pctredfmey_incrmt), 
                    eqm_yield_incrmt = mean(eqm_yield_incrmt, na.rm = T), 
                    eqm_profit_incrmt = mean(eqm_profit_incrmt, na.rm = T), 
                    marg_yield = mean(marg_yield, na.rm = T),
                    marg_profit = mean(marg_profit, na.rm = T),
                    pctred = i,
                    marg_pctred = incrmt/dt$wt
                    ) 
      # }) %>%
      # bind_rows()
  
    
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
## NEW
incrmt <- 1 ## one percent increment in bycatch reduction  #0.1 ## i.e. one tenth of a percent

bycatch_func <- 
  function(z){
    
    final_df <- 
      pblapply(1:n1, function(j) { ## big loop to sample over distributions
        
        ## NEW
        j <<- j
        i <<- 0
        shortfall_msy <<- 1 
        shortfall_mey <<- 1
        pctredmsy_incrmt <<- 1
        pctredmey_incrmt <<- 1
        
        ## NEW
        weighted_df <- list()
        while((shortfall_msy > 0 | shortfall_mey > 0) & (pctredmsy_incrmt < 100 & pctredmey_incrmt < 100)){
        # while(i < 25){ ## test
          weighted_df[[i+1]] <-
            
            ## Produce a data frame with one row
            lapply(
              extract_func(z), ## get lists of relevant target stocks based on bycatch species
              function(x) disb_func(x, n2)#n1, n2) 
              ) %>%
              Reduce("+", .) %>% ## collapse into single weighted row(s)
              as_data_frame() %>%
              mutate(species = z)
          
          ## shortfall = required mean reduction - estimated reduction from sampling function
          shortfall_msy <<- filter(bycatch_df, species == z)$pctredbpt - weighted_df[[i+1]]$pctredmsy_incrmt
          shortfall_mey <<- filter(bycatch_df, species == z)$pctredbpt - weighted_df[[i+1]]$pctredmey_incrmt
          ## Should not be able to reduce beyond 100%
          pctredmsy_incrmt <<- weighted_df[[i+1]]$pctredmsy_incrmt
          pctredmey_incrmt <<- weighted_df[[i+1]]$pctredmey_incrmt
          
          weighted_df[[i+1]]$req_red <- filter(bycatch_df, species == z)$pctredbpt
          weighted_df[[i+1]]$shortfall_msy <- shortfall_msy #ifelse(shortfall_msy>=0, shortfall_msy, NA)
          weighted_df[[i+1]]$shortfall_mey <- shortfall_mey #ifelse(shortfall_mey>=0, shortfall_mey, NA)
          
          ## NEW
          i <<- i + 1
          
        } ## end of while loop NEW
        
        # Cleanup
        rm(i, shortfall_msy, shortfall_mey, pctredmsy_incrmt, pctredmey_incrmt,
           envir = .GlobalEnv)
        
        weighted_df <- 
          weighted_df %>% 
          bind_rows() %>% ## combine list of increments into single data frame
          select(req_red, everything()) 
        
        yield_temp <-
          (weighted_df %>%
            arrange(marg_yield) %>%
            mutate(total_pctred = cumsum(marg_pctred),
                   total_yield = cumsum(marg_yield)
                   ) %>%
             mutate(total_yield = ifelse(total_pctred<=req_red, total_yield, NA)) %>%
             summarise(yield = max(total_yield, na.rm = T)))$yield
        
        profit_temp <-
          (weighted_df %>%
             arrange(marg_profit) %>%
             mutate(total_pctred = cumsum(marg_pctred),
                    total_profit = cumsum(marg_profit)) %>%
             mutate(total_profit = ifelse(total_pctred<=req_red, total_profit, NA)) %>%
             summarise(profit = max(total_profit, na.rm = T)))$profit
        

        ## Finally, collapse into a one-row data frame with the following columns:
        ## species, req_red, pctredmsy, pctredmey, yield, profit
        weighted_df <-
          weighted_df %>%
          group_by(species) %>% ## just to keep species column when summarising next
          summarise(
            req_red = mean(req_red, na.rm = T),
            pctredmsy = mean(pctredmsy, na.rm = T),
            pctredmey = mean(pctredmey, na.rm = T)
            ) %>%
            mutate(
              yield = ifelse(req_red <= pctredmsy, 0, yield_temp),
              profit = ifelse(req_red <= pctredmey, 0, profit_temp)
            )
        
        return(weighted_df)
        
      }) %>% ## end of big loop over sample distribution
      bind_rows() 
    
    ## Clean up
    rm(j,
       envir = .GlobalEnv)
    
    return(final_df)
  }

## E.g. bycatch_func("Loggerhead_turtle")

##############################################################################
### Function 2: Plot that asks: Will rebuilding be enough to stop decline? ###
##############################################################################

bycatchdistggplot <-
  function(bdist){
    
    df <- 
      left_join(
        bdist %>% select(pctredmsy,pctredmey, species), 
        bycatch_df
        ) %>% 
      group_by(species)
    
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

##############################################################################################
### Function 3: Plot that maps the expected yield and profit costs to reach bycatch target ###
##############################################################################################

costggplot <-
  function(df){
    
    df %>%
      select(species, yield, profit) %>%
      rename(Yield = yield, Profit = profit) %>%
      gather(key, cost, -species) %>%
      ggplot() +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = cost, col = key, fill = key), alpha = .5) +
      geom_vline(aes(xintercept = 0), lty = 2) +
      labs(x = "", y = "Density") +
      # xlim(0, 100) +
      scale_color_brewer(name = "", palette = "Set1") +
      scale_fill_brewer(name = "", palette = "Set1") +
      facet_wrap(~key, scales= "free") +
      theme(
        #text = element_text(family = font_type),
        legend.position = "none",
        legend.title = element_blank(),
        # strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "white"), ## Facet strip
        panel.margin = unit(2, "lines") ## Increase gap between facet panels
      ) 
  } 