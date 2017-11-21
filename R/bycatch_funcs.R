#########################################################################
### CONVENIENCE FUNCTION THAT SETS PARAMETERS FOR SPECIFIC MODEL RUNS ###
### AND LOADS CORRECT VERSION OF THE UPSIDES DATA.                    ###  
#########################################################################

choose_run <- 
  function(run){
    
    r_types <- 
      c("main","fcorrected","conservation","alpha=05","alpha=2",
        "nonei","sensrange95","weights","2012only")
    
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  {
        ifelse(is.double(x), abs(x - round(x) < tol), F) 
      }
    
    if(is.wholenumber(run)) {run <- r_types[run]}
    
    run <- tolower(run)
    
    ## 1. Main run (default)
    upsides_type <<- "main" ## C.f. Runs 6 and 9
    corr_factor <<- 1 ## C.f. Run 2
    scenario <<- "All stocks" ## C.f. Run 3
    alpha_exp <<- 1  ## C.f. Runs 4 and 5
    sensrange95 <<- 0 ## C.f. Run 7
    weights_sens <<- 0 ## C.f. Run 8
    
    
    ## 2. Correction factor for possible upward bias in catch-MSY projections
    ##    Here we assume that current F is twice as big as estimated by Costello 
    ##    et al. (2016).
    if(run=="fcorrected") {corr_factor <<- 2}
    
    ## 3. Conservation concern scenario: Only stocks going to MEY/MSY are those currently 
    ##    with F > Fmsy or B < Bmsy.
    if(run=="conservation") {scenario <<- "Con. Concern"}
    
    ## 4. and 5. Different "alpha" exponents, i.e. elasticity of (changes in) bycatch to 
    ##           (changes in) target stocks. See eqn (S14) in the paper.
    if(run=="alpha=05") {alpha_exp <<- 0.5}
    if(run=="alpha=2") {alpha_exp <<- 2}
    
    ## 6. Remove "NEI" (not elsewhere included) stocks from the Costello et al. 2016
    if(run=="nonei") {upsides_type <<- "nonei"}
    
    ## 7. Simulate over a 95% uncertainty range in Fe and delta
    if(run=="sensrange95") {sensrange95 <<- 1}
    
    ## 8. Simulate over a (uniform) 25% uncertainty range in bycatch weights when doing
    ##    cost analysis.
    if(run=="weights") {weights_sens <<- 1}
    
    ## 9. Use F from 2012 only. 
    if(run=="2012only") {upsides_type <<- "2012only"}
    
    ## Convenience strings for import/export, file-reading and naming conventions
    corr_str <- ""
    scenario_str <- ifelse(scenario=="All stocks", "", "_conservation")
    alpha_str <- ifelse(alpha_exp==1, "", gsub("\\.","",paste0("_alpha=",alpha_exp))) 
    upsides_str <- ifelse(upsides_type=="main", "", paste0("_", upsides_type))
    sensrange_str <- ifelse(sensrange95==0, "", "_sensrange95")
    weights_str <- ifelse(weights_sens==0, "", "_weights")
    
    ## Read in the relevant target stock data, derived from the "upsides" model of 
    ## Costello et al. (PNAS, 2016).
    upsides <<- 
      fread(paste0("Data/upsides_", upsides_type, ".csv")) %>% 
      as_data_frame()
    # upsides <<- read_csv(paste0("Data/upsides_", upsides_str, ".csv"), col_types = cols(regionfao = "c"))
    
    ## Final adjustment to upsides data in case of the "fcorrected" run
    if(run=="fcorrected") {
      corr_str <- "_fcorrected"
      upsides <<-
        upsides %>%
        mutate(curr_f = ifelse(dbase=="FAO", curr_f/corr_factor, curr_f)) %>%
        ## Adjust additionally affected variables in sequence
        mutate(
          fvfmsy = curr_f/g,
          eqfvfmey = curr_f/f_mey
        ) %>%
        mutate(
          pctredfmsy = 100 * (1 - (1/fvfmsy)),
          pctredfmey = 100 * (1 - (1/eqfvfmey)),
          pctredfmsycon = 100 * (1-(fconmsy/curr_f))
        )
    }
    
    ## Combined convenience variable for output file name suffix
    suff_str <<- paste0(upsides_str, alpha_str, corr_str, scenario_str, sensrange_str, weights_str)
    
    ## Descriptive message for user
    if(match(run,r_types) == 1) {
      message(paste0("Run ", match(run,r_types), " (", run, ")"))
    } else {
      message(paste0("Run ", match(run,r_types), " (", run, ", sensitivity analysis)"))
    } 
    
    message(
      paste0(
        "Parameters: \n",
        "   upsides_type = ", upsides_type, "\n",
        "   corr_factor = ", corr_factor, "\n",
        "   scenario = ", scenario, "\n",
        "   alpha_exp = ", alpha_exp, "\n",
        "   sensrange95 = ", sensrange95, "\n",
        "   weights_sens = ", weights_sens
        )
      )
  }


##################################################################
##################################################################
### Equilibrium profit and yield functions (for cost analysis) ###
##################################################################
##################################################################

# Equilibrium profit, as a function of f
# Based on eqn 8 in Costello et al. 2016 (PNAS; ref. 2) SI
eqprofitf <- 
  function(f,price,marginalcost,g,k,phi,beta) {
    ymsy <- ((g * k)/((phi + 1)^(1/phi)))
    bvbm <- ((1 + phi - ((f/g) * phi))^(1/phi))
    eqp <- 
      (price * ymsy * (f/g) * bvbm) -
      (marginalcost * ((f)^beta))
    return(eqp)
  }

# Equilibrium yield, as a function of f
eqyieldf <- 
  function(f,g,k,phi) {
    eqy <- 
      f *
      (k * ((1 - ((f * phi)/(g * (phi + 1))))^(1/phi)))
    return(eqy)
  }


# Equilibrium profit cost of fishing at f, relative to MEY:
profitcostf <- 
  function(f,price,marginalcost,fmey,g,k,phi,beta) {
    pc <-
      eqprofitf(fmey,price,marginalcost,g,k,phi,beta) -
      eqprofitf(f,price,marginalcost,g,k,phi,beta)
    return(pc)
  }

# Equilibrium yield cost of fishing at f, relative to MSY:
yieldcostf <- 
  function(f,fmsy,g,k,phi) {
    yc <- 
      eqyieldf(fmsy,g,k,phi) -
      eqyieldf(f,g,k,phi)
    return(yc)
  }

# Marginal equilibrium yield, as a function of f:
myieldf <- 
  function(f,g,k,phi) {
    my <- 
      (((g - f) * k *
          ((1 - ((f * phi)/(g * (phi + 1))))^((1/phi) - 1)))
       /g)
    return(my)
  }

# Marginal equilibrium profit, as a function of f:
mprofitf <- 
  function(f,price,marginalcost,g,k,phi,beta) {
    mp <-
      ((price * ((1 + phi - ((f * phi)/g))^(1/phi)) * (f - g) * k * ((1 + phi)^((phi - 1)/phi)))/
         ((f * phi) - (g * (1 + phi)))) -
      (beta * marginalcost * (f^(beta - 1)))
    return(mp)
  }


###############################
###############################
### Cost analysis functions ###
###############################
###############################

##1. Given dataframe of sampled target stocks, and marginal cost of bycatch reduction 
#          (denoted 'mpb' for profit, 'myb' for yield), return total cost.

#1.1 Define numerical inverse function
inverse = function (f, lower = 0, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]$root
}
# inverse <- possibly(inverse, NA_real)

#1.2 Inverse marginal cost (in eq. yield or profit) of reducing f for one target stock:
#     Input is marginal cost of reducing f (denoted 'mp', equal to marginal profit from increasing f, 
#           denoted pi'(Fij) in SM text ('my' is analog for yield)).
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

#1.3 Function that takes the sample dataframe of target stocks (df), 
#      and a marginal yield or profit of bycatch reduction (mpb or myb), 
#      and computes the pct reduction for the bycatch species, and the total cost 
#      (as a percentage of total MSY or MEY), given that marginal cost.

redncost_giv_mpb <- function(df, mpb) { # here, mpb is the marginal profit cost of bycatch reduction.
    
  if (scenario == "All stocks") { # compare to MEY for all stocks for reduction and cost
    dt <- 
      df %>% 
      group_by(1:n()) %>%
      # Calculate fishing mortality for each stock (f_mpb) that corresponds to the marginal profit cost, 'mpb'
      # wgt * mpb = mp (from eq. S9 in the SM, where wgt = wj/n2)
      mutate(f_mpb = inv_marg_profit(wgt * mpb, f_mey, price, marginalcost, g, k, phi, beta)) %>%
      ungroup() %>%
      # Create new column 'pctred_b' which specifies the % reduction in fishing
      # mortality for the bycatch stock resulting from the % reduction
      # for each target stock (relative to base year) at f_mpb.
      mutate(pctred_b = wgt * 100 * (1 - ((f_mpb/curr_f)^alpha_exp))) %>%
      # Create new column 'mey' which specifies the mey for
      # each target stock. Because scenario is "All stocks", mey assumes Fmey for all target stocks. 
      mutate(mey = eqprofitf(f_mey,price,marginalcost,g,k,phi,beta)) %>%
      # Create new column 'pcost' which specifies the lost profit relative to mey for
      # each target stock.
      mutate(pcost = profitcostf(f_mpb,price,marginalcost,f_mey,g,k,phi,beta))
  } else { # compare to Con. Concern (i.e. only Con. Concern stocks at MEY, no other stocks can fish harder)
    dt <- 
      df %>% 
      group_by(1:n()) %>%
      # Calculate fishing mortality for each stock that corresponds to the marginal profit cost, 'mpb'. 
      mutate(f_mpb1 = inv_marg_profit(wgt * mpb, f_mey, price, marginalcost, g, k, phi, beta)) %>%
      # If f_mpb1 exceeds fconmey, set f_mpb = fconmey, because non-con. concern stocks cannot be fished harder
      # in this scenario.
      mutate(f_mpb = if_else(f_mpb1 < fconmey, f_mpb1, fconmey)) %>%
      select(-f_mpb1) %>%
      ungroup() %>%
      mutate(pctred_b = wgt * 100 * (1 - ((f_mpb/curr_f)^alpha_exp))) %>%
      # Create new column 'mey' which specifies the mey for
      # each target stock. Because scenario is "Con. Concern", mey assumes Fmey for all overfished stocks, 
      # and constant biomass for other stocks, following Costello et al. 2016.
      mutate(mey = eqprofitf(fconmey,price,marginalcost,g,k,phi,beta)) %>%
      mutate(pcost = profitcostf(f_mpb,price,marginalcost,fconmey,g,k,phi,beta))
    }
  
  # Create empty data frame to contain the results.
  output <- data_frame(pctred = 0, cost = 0)
  
  # % reduction in mortality for bycatch species (pctred), given marginal profit cost, (mpb),
  # is the sum of pctred_b across all target stocks in df.
  output$pctred <- sum(dt$pctred_b)
  
  # 'cost' is the difference between the total profit at the f's having marginal profit cost, mpb,
  # (and therefore required for % bycatch reduction 'pctred') and mey, as a fraction of mey,
  # cumulatively across all target stocks in df.
  output$cost <- 100 * (sum(dt$pcost)/sum(dt$mey))
  return(output)
}

redncost_giv_myb <- function(df, myb) { # function is analogous to 'redncost_giv_mp' above
  
  if (scenario == "All stocks") {
    dt <- 
      df %>% 
      group_by(1:n()) %>%
      mutate(f_myb = inv_marg_yield(wgt * myb, g, k, phi)) %>%
      ungroup() %>%
      mutate(pctred_b = wgt * 100 * (1 - ((f_myb/curr_f)^alpha_exp))) %>%
      mutate(msy = eqyieldf(g,g,k,phi)) %>%
      mutate(ycost = yieldcostf(f_myb,g,g,k,phi))
  } else {
    dt <- 
      df %>% 
      group_by(1:n()) %>%
      mutate(f_myb1 = inv_marg_yield(wgt * myb, g, k, phi)) %>%
      mutate(f_myb = if_else(f_myb1 < fconmsy, f_myb1, fconmsy)) %>%
      select(-f_myb1) %>%
      ungroup() %>%
      mutate(pctred_b = wgt * 100 * (1 - ((f_myb/curr_f)^alpha_exp))) %>%
      mutate(msy = eqyieldf(fconmsy,g,k,phi)) %>%
      mutate(ycost = yieldcostf(f_myb,fconmsy,g,k,phi))
  }
  
  output <- data_frame(pctred = 0, cost = 0)
  output$pctred <- sum(dt$pctred_b)
  output$cost <- 100 * (sum(dt$ycost)/sum(dt$msy))
  return(output)
}

##2. Given dataframe of sampled target stocks ('df') 
#     and pctredpt (the bycatch reduction target, called 'pctredb' in the function), 
#     calculate marginal profit costs (y and p).

myb_calc <- function(df, pctredb) {
  dt <- 
    df %>%
    # Generate new column in df of maximum marginal yield cost (i.e. myb at 0 fishing)
    mutate(maxmyb = myieldf(0,g,k,phi)/wgt)
  
  # store largest possible 'maxmyb' among stocks in df as 'mxmyb'
  mxmyb <- max(dt$maxmyb)
  
  # Create 'imy' function which calculates the marginal yield cost, given pctredb
  imyb <- inverse(function (myb) redncost_giv_myb(dt, myb)$pctred, 0, mxmyb)
  # Check if the maximum marginal yield cost ('mxmy') corresponds 
  # to a still-insufficient reduction in bycatch.
  ifelse(
    redncost_giv_myb(dt, mxmyb)$pctred < pctredb,
    # If yes, then require the max marginal yield.
    output <- mxmyb,
    # If no, use 'imyb' function to calculate required marginal yield
    output <- imyb(pctredb)
  ) 
  return(output)
}

mpb_calc <- function(df, pctredb) { # function is analogous to 'myb_calc' above
  dt <- df %>%
    mutate(maxmpb = mprofitf(0,price,marginalcost,g,k,phi,beta)/wgt)
  mxmpb <- max(dt$maxmpb)
  impb <- inverse(function (mpb) redncost_giv_mpb(dt, mpb)$pctred, 1, mxmpb)
  ifelse(
    redncost_giv_mpb(dt, mxmpb)$pctred < pctredb,
    output <- mxmpb,
    output <- impb(pctredb)
  ) 
  return(output)
}

##3 Main function: Takes as input a dataframe of sampled stocks 
#                  and relevent columns from these, as well as total reduction needed
#                  for bycatch species, and percent reduction from MSY/MEY
#                  Returns 0 if there is no shortfall. 
#                  Returns 100 if pctredpt > 100 (i.e. if stopping all fishing isn't enough to stop decline). 
#                  Returns profit (yield) cost as fraction of total MEY profit (MSY yield) otherwise.

cost_yield <- function(df, pctredb, meanpctredmsy) {
  
  ifelse(
    # Check if stopping all fishing is enough for bycatch target.
    # If not, cost = 100 (percent of fishing yield)
    round(pctredb) >= 100, 
    ycost <- 100,
    ifelse(
      # Check if there is a shortfall.
      round(pctredb) < round(meanpctredmsy), 
      # If not, cost = 0 (because MSY reduces bycatch enough to stop decline)
      ycost <- 0,
      # There is a shortfall, and it is attainable through additional target species F reductions.
      # Code below calculates costs. It is assumed that df includes weights of each target stock (which sum to 1)
      #1. Calculate marginal yield cost (with function 'my_calc'), 
      #     given pct reduction needed for bycatch species.
      #2. Calculate total yield cost, given marginal cost calculated in 1.
      ycost <-
        tryCatch(
          redncost_giv_myb(df, myb_calc(df, pctredb))$cost,
          error = function(e) as.double(NA)
          )
    )
  )
  return(ycost)
}

cost_profit <- function(df, pctredb, meanpctredmey) {
  ifelse(
    # Check if stopping all fishing is enough for bycatch target.
    # If not, cost = 100 (percent of fishing profits)
    round(pctredb) >= 100, 
    pcost <- 100,
    ifelse(
      # Check if there is a shortfall.
      round(pctredb) < (meanpctredmey), 
      # If not, cost = 0 (because MEY reduces bycatch enough to stop decline)
      pcost <- 0,
      # Else, there is a shortfall, and it is attainable through additional target species F reductions.
      # Code below calculates costs. It is assumed that df includes weights of each target stock (which sum to 1)
      #1. Calculate marginal profit cost (with function 'mp_calc'), 
      #     given pct reduction needed for bycatch species.
      #2. Calculate total profit cost, given marginal cost calculated in 1. 
      pcost <-
        tryCatch(
          redncost_giv_mpb(df, mpb_calc(df, pctredb))$cost, 
          error = function(e) as.double(NA)
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
        fconmsy > 0,
        fconmey > 0
        ) %>%
      mutate(sampling_wgt = marginalcost * ((curr_f)^beta)) %>%
      select(
        dbase,idorig,idoriglumped,pctredfmsy,pctredfmey,pctredfmey,pctredfmsycon,pctredfmeycon,
        sampling_wgt,speciescat,speciescatname,fmeyvfmsy,k,fvfmsy,g,beta,phi,price,marginalcost,eqfvfmey,
        curr_f,f_mey,fconmsy,fconmey
        ) %>%
      mutate(trgcat = dt$target) %>%
      mutate(wt = dt$wt) %>%
      mutate(bycsp = dt$species) %>%
      mutate(sampling_wgt = sampling_wgt/sum(sampling_wgt, na.rm = T))
    
    return(stocks_df)
  }

########################################################################
## Step 3: Produce a data frame of the sample distributions and costs ##
########################################################################

## Function for sampling within a particular target stock category (e.g. demersal).
## Sampling is weighted according to the marginal costs of that stock (defined in 
## 'stocks_df' above).
samp_func <- 
  function(dt, n2, reltdf) {
    smple <- reltdf %>%
      filter(trgcat == dt$target) %>%
      sample_n(n2, replace = T, weight = sampling_wgt)
    return(smple)
  }

## Function to randomly adjust bycatch weights within +/- 25% uniform range.
## Only relevant to "weights" sensitivity run.
wt_func <-
  function(wts) {
    prbs <- runif(length(wts), 0.75, 1.25)
    new_wts <- wts * prbs / sum(wts)
    new_wts_normalised <- new_wts / sum(new_wts)
    return(new_wts_normalised)
  }

## Inputs are a list of "dt's" (i.e. the output from `extract_func`) and the set 
## of parameters for determining %T. This latter metric will either be a point
## estimate as in the main run, or drawn from underlying parameter distributions 
## as in the "sensrange95" run.
single_worldstate_outputs <- 
  # function(dt2, n2, pctredb, pctredbl, pctredbu, reltdf, sensrangept95, sensrange_type) { 
  ### REVISION
  function(
    dt2, n2, reltdf, 
    delta_mean, delta_q025, delta_q975, deltaN_mean, deltaN_q025, deltaN_q975, fe_mean, fe_q025, fe_q975, 
    sensrange95, sensrange_type
    ) { ### END REVISION
    
    ## Sample within each of the relevant target categories (demersal, shrimp, etc.) 
    ## and the combine into a common data frame.
    samp <- 
      lapply(dt2, function (x) samp_func(x, n2, reltdf)) %>%
      bind_rows()
    
    ## weights adjustment
    if(weights_sens==1){
      samp <-
        samp %>%
        distinct(trgcat, wt) %>%
        mutate(wt = wt_func(wt)) %>%
        right_join(samp %>% select(-wt))
    }
    ## end weights adjustment
    
    ## The impact of each stocks is then weighted according to bycatch estimates 
    ## for its overall target category (taken from the literature)
    samp <- 
      samp %>%
      mutate(wgt = wt * (1/n2))
    
   ## Conservation concern scenario adjustments
    if (scenario == "All stocks") { ## i.e. scenario != "Con. concern"
    samp <- 
      samp %>%
      mutate(
        pctrmsywt = wgt * 100 * (1 - ((1 - (pctredfmsy/100))^alpha_exp)),
        pctrmeywt = wgt * 100 * (1 - ((1 - (pctredfmey/100))^alpha_exp))
        )
    } else { ## i.e. scenario == "Con. concern"
      samp <- 
        samp %>%
        mutate(
          pctrmsywt = wgt * 100 * (1 - ((1 - (pctredfmsycon/100))^alpha_exp)),
          pctrmeywt = wgt * 100 * (1 - ((1 - (pctredfmeycon/100))^alpha_exp))
        )
    }
    ## End conservation concern scenario adjustments
    
    mpctmsy <- sum(samp$pctrmsywt, na.rm = T) ## GRM: Added na.rm = T
    mpctmey <- sum(samp$pctrmeywt, na.rm = T) ## GRM: Added na.rm = T
    
    # if (sensrangept95 == 1) {
    #   pctb <- runif(1, min = pctredbl, max = pctredbu) # if 95% uncertainty in Fe and delta is on
    #   } else {
    #   pctb <- pctredb
    #   }
    ### REVISION
    ## Parameter uncertainty / 95% CI sensitivity scenario adjustments
    ## Calculate "%T" according to whether we are just taking the mean values
    ## of delta, deltan and Fe as given point estimates... Or sampling from the 
    ## full underlying parameter distributions ("sensrange95" run only).
    if (sensrange95 == 0) { ## Normal run. Don't consider uncertainty in %T
      # pctb <- pctredb
      if(is.na(deltaN_mean)) {
        deltaN_mean <- delta_mean + fe_mean ## Should only be relevent for type 4 cases
      }
      pctT <- 100 * ((fe_mean-deltaN_mean) / fe_mean)
    } else { ## The "sensrange95" run. Sample %T according to distribution of underlying parameters
      # pctb <- runif(1, min = pctredbl, max = pctredbu) # if 95% uncertainty in Fe and delta is on
      delta_sd <- abs((delta_mean - delta_q025)) / 1.96
      deltaN_sd <- abs((deltaN_mean - deltaN_q025)) / 1.96
      fe_sd <- abs((fe_mean - fe_q025)) / 1.96
      ## Take draws from assumed normal distributions. Will override with 
      ## uniform draws for sensrange_type(s) 1 and 5.
      delta_draw <- suppressWarnings(rnorm(1, mean = delta_mean, sd = delta_sd))
      deltaN_draw <- suppressWarnings(rnorm(1, mean = deltaN_mean, sd = deltaN_sd))
      fe_draw <- suppressWarnings(rnorm(1, mean = fe_mean, sd = fe_sd))
      ## Type 1: delta~N(.) & deltaN~U(.) 
      if(sensrange_type==1) {
        deltaN_draw <- runif(1, min = max(delta_draw , deltaN_q025), max = deltaN_q975) ## Make sure deltaN_draw>=delta_draw
        pctT <- 100 * (delta_draw / (delta_draw-deltaN_draw))
      }
      ## Type 2: Fe~N(.) & deltaN~N(.)
      if(sensrange_type==2) {
        pctT <- 100 * ((fe_draw-deltaN_draw) / fe_draw)
      }
      ## Type 3: delta~N(.) & deltaN~N(.)
      if(sensrange_type==3) {
        if(deltaN_draw<=delta_draw){
          deltaN_draw <- truncnorm::rtruncnorm(1, a=delta_draw, b=Inf, mean = deltaN_mean, sd = deltaN_sd) ## Make sure deltaN_draw>=delta_draw
        }
        pctT <- 100 * (delta_draw / (delta_draw-deltaN_draw))
      }
      ## Type 4: delta~N(.) & Fe~N(.)
      if(sensrange_type==4) {
        pctT <- 100 * (-delta_draw / fe_draw)
      }
      ## Type 5: delta~U(.) & deltaN~U(.)
      if(sensrange_type==5) {
        delta_draw <- runif(1, min = delta_q025, max = delta_q975)
        deltaN_draw <- runif(1, min = max(delta_draw , deltaN_q025), max = deltaN_q975) ## Make sure deltaN_draw>=delta_draw
        pctT <- 100 * (delta_draw / (delta_draw-deltaN_draw))
      }
    }
    ## Manual correction if randomly ended up dividing by zero (v. low probability of occuring)
    if(pctT==Inf) {pctT <- 0}
    ## End parameter uncertainty / 95% CI sensitivity scenario adjustments
    ### END REVISION
    
    stwld <- 
      data_frame(
        pctredmsy = mpctmsy, 
        pctredmey = mpctmey,
        # ycostmsy = cost_yield(samp, pctb, mpctmsy),
        # pcostmey = cost_profit(samp, pctb, mpctmey)
        ycostmsy = cost_yield(samp, pctT, mpctmsy), ### REVISION changed from pctb to pctT
        pcostmey = cost_profit(samp, pctT, mpctmey), ### REVISION changed from pctb to pctT
        pctT = pctT ### REVISION added
        ) # need pctredpt from somewhere
    return(stwld)
  }

## Input is list of 'dt's' -> output from 'extract_func'
disb_func <-
  function(dt2, n1, n2){
    
    ##########################################################
    ## Step 3.1: Create dataframe of relevant target stocks ##
    ##########################################################
    rel_targets <- 
      lapply(dt2, upsides_subset_func) %>%
      bind_rows()
    
    ## Required reduction parameters (i.e. What is the rate of population decline for this bycatch species?)
    # pctredbt <- (filter(bycatch_df, species==rel_targets$bycsp[1]))$pctredbpt[1] ## Point estimate
    # pctredbtl <- (filter(bycatch_df, species==rel_targets$bycsp[1]))$pctredbl[1] ## Lower bound
    # pctredbtu <- (filter(bycatch_df, species==rel_targets$bycsp[1]))$pctredbu[1] ## Upper bound
    ### REVISION
    ## Collapse bycatch_df to single species of interest/relevance
    rel_bycsp <- filter(bycatch_df, species==rel_targets$bycsp[1])
    ## Now extract the relevant parameters for determining (the distribution 
    ## around "%T". See Eq. (2b) in the paper. Only the mean values will be used 
    ## in the main run. The full distributions (i.e. 95% CIs) will be used in  
    ## the "sensrange95" run.
    delta_mean <- rel_bycsp$delta_mean
    delta_q025 <- rel_bycsp$delta_q025
    delta_q975 <- rel_bycsp$delta_q975
    deltaN_mean <- rel_bycsp$deltaN_mean
    deltaN_q025 <- rel_bycsp$deltaN_q025
    deltaN_q975 <- rel_bycsp$deltaN_q975
    fe_mean <- rel_bycsp$fe_mean
    fe_q025 <- rel_bycsp$fe_q025
    fe_q975 <- rel_bycsp$fe_q975
    ## Similarly, what combination of underlying parameter distributions govern the 
    ## overall uncertainty in "%T"? (Only relevant to the "sensrange95" run.)
    sensrange_type <- rel_bycsp$sensrange_type
    ## END REVISION
    
    
    #########################################################################################################
    ## Step 3.2: Repeatedly sample from stocks data frame to create distribution of both pctreds and costs ##
    #########################################################################################################
    
    dists <-
      pblapply(1:n1, 
               possibly(function(i) {
                 evalWithTimeout(
                   # single_worldstate_outputs(dt2, n2, pctredbt, pctredbtl, pctredbtu, rel_targets, sensrange95),
                   ### REVISION
                   single_worldstate_outputs(
                     dt2, n2, rel_targets, 
                     delta_mean, delta_q025, delta_q975, deltaN_mean, deltaN_q025, deltaN_q975, fe_mean, fe_q025, fe_q975, 
                     sensrange95, sensrange_type
                     ), 
                   ### END REVISION
                   timeout = 20, ## i.e. Time out after 20 seconds if can't resolve 
                   TimeoutException = function(ex) "TimedOut"
                   )
                 }, 
                 # otherwise = NULL
                 otherwise = data_frame(pctredmsy=as.double(NA), pctredmey=as.double(NA), ycostmsy=as.double(NA), pcostmey=as.double(NA)) ## To catch failed uniroot cases
                 ),
               cl=num_cores
               )
    return(dists)
    }

###############################################################################
## Final step: Convenience wrapper to apply over all target stocks affecting ##
## a single bycatch species. This is what we'll call in the actual analysis. ##
###############################################################################

bycatch_func <-
  function(z){
    message(paste0(paste0((bycatch_df %>% mutate(n=row_number()) %>% filter(species==z))$n), ". ", z))
    disb_func(extract_func(z), n1, n2) %>%
      bind_rows() %>%
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
  function(bdist, series=NULL){
    
    df1 <- 
      left_join(bdist, bycatch_df) %>% 
      group_by(species) %>%
      # select(-ycostmsy, -pcostmey) %>%
      select(-ycostmsy, -pcostmey, -pctT) %>% ### REVISION
      mutate_if(is.double, funs(. / 100)) 
    
    df2 <-
      df1 %>%
      rename(
        MSY = pctredmsy,
        MEY = pctredmey
        ) %>%
      select(MSY:species) %>%
      gather(key, pctred, -species) %>%
      mutate(key = factor(key, levels = c("MSY", "MEY"))) 
    
    ### REVISION
    pctT_df <- df1 %>% distinct(species, .keep_all = T) %>% ungroup()
    pctT_df <-
      lapply(seq_len(nrow(pctT_df)), function(i) {
        rel_bycsp <- slice(pctT_df, i)
        bycsp <- rel_bycsp$species
        # N <- df1 %>% filter(species==bycsp) %>% count(species) %>% pull(n)
        N <- 10000
        ##
        delta_mean <- rel_bycsp$delta_mean
        delta_q025 <- rel_bycsp$delta_q025
        delta_q975 <- rel_bycsp$delta_q975
        deltaN_mean <- rel_bycsp$deltaN_mean
        deltaN_q025 <- rel_bycsp$deltaN_q025
        deltaN_q975 <- rel_bycsp$deltaN_q975
        fe_mean <- rel_bycsp$fe_mean
        fe_q025 <- rel_bycsp$fe_q025
        fe_q975 <- rel_bycsp$fe_q975
        ## 
        sensrange_type <- rel_bycsp$sensrange_type
        ##
        if(is.na(deltaN_mean)) {
          deltaN_mean <- delta_mean + fe_mean ## Should only be relevent for type 4 cases
        }
        pctT_mean <- 100 * ((fe_mean-deltaN_mean) / fe_mean)
        ##
        delta_sd <- abs((delta_mean - delta_q025)) / 1.96
        deltaN_sd <- abs((deltaN_mean - deltaN_q025)) / 1.96
        fe_sd <- abs((fe_mean - fe_q025)) / 1.96
        ## Take draws from assumed normal distributions. Will override with 
        ## uniform draws for sensrange_type(s) 1 and 5.
        delta_draw <- suppressWarnings(rnorm(N, mean = delta_mean, sd = delta_sd))
        deltaN_draw <- suppressWarnings(rnorm(N, mean = deltaN_mean, sd = deltaN_sd))
        fe_draw <- suppressWarnings(rnorm(N, mean = fe_mean, sd = fe_sd))
        ## Type 1: delta~N(.) & deltaN~U(.) 
        if(sensrange_type==1) {
          # deltaN_draw <- runif(N, min = max(delta_draw , deltaN_q025, na.rm=T), max = deltaN_q975) ## Make sure deltaN_draw>=delta_draw
          deltaN_draw <- runif(N, min = deltaN_q025, max = deltaN_q975)
          # pctT <- 100 * (delta_draw / (delta_draw - deltaN_draw))
        }
        ## Type 2: Fe~N(.) & deltaN~N(.)
        if(sensrange_type==2) {
          # pctT <- 100 * ((fe_draw - deltaN_draw) / fe_draw)
        }
        ## Type 3: delta~N(.) & deltaN~N(.)
        if(sensrange_type==3) {
          # if(deltaN_draw<=delta_draw){
          #   deltaN_draw <- truncnorm::rtruncnorm(N, a=delta_draw, b=Inf, mean = deltaN_mean, sd = deltaN_sd) ## Make sure deltaN_draw>=delta_draw
          # }
          deltaN_draw <- rnorm(N, mean = deltaN_mean, sd = deltaN_sd)
          # pctT <- 100 * (delta_draw / (delta_draw - deltaN_draw))
        }
        ## Type 4: delta~N(.) & Fe~N(.)
        if(sensrange_type==4) {
          # pctT <- 100 * (-delta_draw / fe_draw)
        }
        ## Type 5: delta~U(.) & deltaN~U(.)
        if(sensrange_type==5) {
          delta_draw <- runif(N, min = delta_q025, max = delta_q975)
          # deltaN_draw <- runif(N, min = max(delta_draw , deltaN_q025, na.rm=T), max = deltaN_q975) ## Make sure deltaN_draw>=delta_draw
          deltaN_draw <- runif(N, min = deltaN_q025, max = deltaN_q975)
          # pctT <- 100 * (delta_draw / (delta_draw - deltaN_draw))
        }
        
        T_df <-
          # data.frame(pctT = pctT/100) %>%
          data.frame(
            delta_draw = delta_draw,
            deltaN_draw = deltaN_draw,
            fe_draw = fe_draw
            ) %>%
          as_data_frame() %>%
          mutate(
            species = bycsp,
            sensrange_type = sensrange_type,
            pctT_mean = pctT_mean/100,
            delta_mean = delta_mean,
            delta_sd = delta_sd,
            delta_q025 = delta_q025,
            delta_q975 = delta_q975,
            deltaN_mean = deltaN_mean,
            deltaN_sd = deltaN_sd,
            deltaN_q025 = deltaN_q025,
            deltaN_q975 = deltaN_q975,
            fe_mean = fe_mean,
            fe_sd = fe_sd,
            fe_q025 = fe_q025,
            fe_q975 = fe_q975
            ) %>% 
          ## Correction 
          mutate(
            deltaN_draw = ifelse(
              deltaN_draw < delta_draw,
              ifelse(
                sensrange_type==3,
                truncnorm::rtruncnorm(1, a=delta_draw, b=Inf, mean = deltaN_mean, sd = deltaN_sd),
                ifelse(
                  sensrange_type %in% c(1, 5),
                  runif(1, min = max(delta_draw , deltaN_q025, na.rm=T), max = deltaN_q975),
                  deltaN_draw
                  )
                ),
              deltaN_draw
              )
            ) %>%
          ## Now calculate pctTs (NOTE: Don't multiply by 100 because plot uses raw percentages)
          mutate(
            pctT = ifelse(
              ## Type 1
              sensrange_type==1,
              delta_draw / (delta_draw - deltaN_draw),
              ifelse(
                ## Type 2
                sensrange_type==2,
                (fe_draw - deltaN_draw) / fe_draw,
                ifelse(
                  ## Type 3
                  sensrange_type==3,
                  delta_draw / (delta_draw - deltaN_draw),
                  ifelse(
                    ## Type 4
                    sensrange_type==4,
                    -delta_draw / fe_draw,
                    ## Type 5
                    delta_draw / (delta_draw - deltaN_draw)
                    )
                  )
                )
              )
            ) %>%
          ## Manual correction if randomly ended up dividing by zero (v. low probability of occuring)
          mutate(pctT = if_else(pctT==Inf, 0, pctT)) %>%
          ## Truncate disbs to 95% range to aid visual inspection
          mutate(
            q025 = quantile(pctT, 0.025, na.rm=T),
            q975 = quantile(pctT, 0.975, na.rm=T)
            ) %>%
          mutate(pctT = ifelse(pctT<q025, NA, pctT)) %>%
          mutate(pctT = ifelse(pctT>q975, NA, pctT)) %>%
          select(species, pctT, pctT_mean)
        
        return(T_df)
        
      }) %>%
      bind_rows()
    ### END REVISION
    
    if(!is.null(series)) df2 <- filter(df2, key==series)
    
    df2 %>% 
      ggplot() +
      # geom_rect(
      #   data = filter(df1, row_number() == 1),
      #   aes(ymin = -Inf, ymax = Inf, xmin = pctredbl, xmax = pctredbu),
      #   alpha = .25
      #   ) +
      geom_density(
        data = pctT_df, 
        aes(x = pctT, y = ..scaled..), col = "gray", fill = "gray", alpha = 0.75,
        adjust = 2#, ## use more smoothing 
        # show.legend = F
        ) +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = pctred, y = ..scaled.., col = key, fill = key), alpha = .5) +
      # geom_vline(
      #   data = filter(df1, row_number() == 1),
      #   aes(xintercept = pctredbpt), lty = 2
      #   ) +
      geom_vline(
        data = pctT_df,
        aes(xintercept = pctT_mean), lty = 2
        ) +
      labs(x = "Reduction in mortality", y = "Density") + 
      scale_x_continuous(labels = percent) + 
      # scale_color_brewer(palette = "Set1") + ## Replaced with below to match MEY filter above
      # scale_fill_brewer(palette = "Set1") + ## Ditto
      scale_color_manual(name = "", values = c("MSY"="#E41A1C", "MEY"="#377EB8")) + ## show_col(brewer_pal(palette = "Set1")(2))
      scale_fill_manual(name = "", values = c("MSY"="#E41A1C", "MEY"="#377EB8")) + ## show_col(brewer_pal(palette = "Set1")(2))
      facet_wrap(~species) +
      theme(legend.position = "bottom") 
  } 

##############################################################
### Plot that asks: How much will it cost to stop decline? ###
##############################################################
cost_plot <-
  function(bdist, series=NULL){
    
    df1 <- 
      left_join(bdist, bycatch_df) %>% 
      group_by(species) %>%
      # select(-pctredmsy, -pctredmey) %>%
      select(-pctredmsy, -pctredmey, -pctT) %>% ### REVISION
      mutate_if(is.double, funs(. / 100)) 
    
    df2 <-
      df1 %>%
      rename(MSY = ycostmsy,MEY = pcostmey) %>%
      select(MSY:species) %>%
      gather(key, pctcost, -species) %>%
      mutate(key = factor(key, levels = c("MSY", "MEY"))) 
    
    if(!is.null(series)) df2 <- filter(df2, key==series)
    
    if(!is.null(series)){
      x_lab <-  paste0("Cost (percent of ", series,")")
    }else{
      x_lab <- "Cost (percent of MSY or MEY)"
    }
    
    df2 %>% 
      ggplot() +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = pctcost, y = ..scaled.., col = key, fill = key), alpha = .5, adjust = 0.01) +
      labs(x = x_lab, y = "Density") + 
      # xlim(0, 100) +
      scale_x_continuous(limits=c(0,1), oob = rescale_none, labels = percent) + 
      # scale_color_brewer(palette = "Set1") + ## Replaced with below to match MEY filter above
      # scale_fill_brewer(palette = "Set1") + ## Ditto
      scale_color_manual(values = c("MSY"="#E41A1C", "MEY"="#377EB8")) + ## show_col(brewer_pal(palette = "Set1")(2))
      scale_fill_manual(values = c("MSY"="#E41A1C", "MEY"="#377EB8")) + ## show_col(brewer_pal(palette = "Set1")(2))
      facet_wrap(~species) +
      theme(legend.position = "bottom") 
  } 

###################################################################################
### Plot that asks: How much would targeting need to improve to stop decline? ###
###################################################################################
targeting_plot <-
  function(bdist, series=NULL){
    
    df1 <- 
      left_join(bdist, bycatch_df) %>% 
      group_by(species) %>%
      mutate(delta_mey = delta + (fe * (pctredmey/100)),
             delta_msy = delta + (fe * (pctredmsy/100))) %>% # calculate growth rate after rebuilding
      mutate(targeting_pct_mey1 = if_else(delta_mey >= 0, 0, 100 * (1-((delta + fe)/(delta + fe - delta_mey)))),
             targeting_pct_mey = if_else(targeting_pct_mey1 > 100, 100, targeting_pct_mey1),
             targeting_pct_msy1 = if_else(delta_msy >= 0, 0, 100 * (1-((delta + fe)/(delta + fe - delta_msy)))),
             targeting_pct_msy = if_else(targeting_pct_msy1 > 100, 100, targeting_pct_msy1)) %>% # calculate selectivity change needed (% reduction in Fe at MEY)
      select(-pctredmsy, -pctredmey, -pcostmey, -ycostmsy, -delta_mey, -delta_msy, -targeting_pct_mey1, -targeting_pct_msy1) %>%
      mutate_if(is.double, funs(. / 100)) 
    
    df2 <-
      df1 %>%
      rename(MSY = targeting_pct_msy, MEY = targeting_pct_mey) %>%
      select(MSY, MEY, species) %>%
      gather(key, targeting_pct, -species) %>%
      mutate(key = factor(key, levels = c("MSY", "MEY"))) 
    
    if(!is.null(series)) df2 <- filter(df2, key==series)
    
    x_lab <- "Req. targeting improvement"
    
    df2 %>% 
      ggplot() +
      # geom_line(stat = "density") + ## lines only
      geom_density(aes(x = targeting_pct, y = ..scaled.., col = key, fill = key), alpha = .5, adjust = 0.01) +
      labs(x = x_lab, y = "Density") + 
      # xlim(0, 100) +
      scale_x_continuous(limits=c(0,1), oob = rescale_none, labels = percent) + 
      # scale_color_brewer(palette = "Set1") + ## Replaced with below to match MEY filter above
      # scale_fill_brewer(palette = "Set1") + ## Ditto
      scale_color_manual(values = c("MSY"="#E41A1C", "MEY"="#377EB8")) + ## show_col(brewer_pal(palette = "Set1")(2))
      scale_fill_manual(values = c("MSY"="#E41A1C", "MEY"="#377EB8")) + ## show_col(brewer_pal(palette = "Set1")(2))
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
      filter(key=="MEY") %>% ## MEY filter added
      ggplot(aes(x = pctred, y = wt_usd, col = key)) +
      geom_point(alpha=0.5, size = 3) +
      # scale_color_brewer(palette = "Set1") + ## Replaced with below to match MEY filter above
      scale_color_manual(values = c("MSY"="#E41A1C", "MEY"="#377EB8")) + ## show_col(brewer_pal(palette = "Set1")(2))
      scale_x_continuous(limits = c(-1, 1), labels = percent) +
      scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
      labs(x = "Reduction in mortality", y = "Effort (USD)") +
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
      gather(key, value, -c(species, pctT)) %>%
      group_by(species, key) %>%
      # do(q=(quantile(.$value))) %>%
      summarise(
        q025 = quantile(value, .025, na.rm = T),
        q50 = quantile(value, .5, na.rm = T),
        q975 = quantile(value, .975, na.rm = T)
      ) %>%
      left_join(bycatch_df, by = 'species') %>%
      select(species, grp, clade, region, contains("pctred"), everything()) %>%
      ungroup()
  }


tradeoffs_plot <-
  function(summ_df, goal) {
    
    if(toupper(goal)=="MEY"){
      df1_filter <- "pctredmey"
      df2_filter <- "pcostmey"
    }else{
      df1_filter <- "pctredmsy"
      df2_filter <- "ycostmsy"
    }
    
    # Make dataframe and calculate new growth rate at mey (median), and selectivity change needed
    df1 <- 
      summ_df %>%
      ## Get pctredmey/pctredmsy rows from results_summary
      filter(key == df1_filter) %>% 
      ## Calculate growth rate post rebuilding
      # mutate(delta_post = delta + (fe * (q50/100))) %>% 
      mutate(delta_post = delta_mean + (fe_mean * (q50/100))) %>% ### REVISION (Should this be deltaN? Same q for below...)
      ## Calculate required selectivity change (% reduction in Fe at MEY or MSY)
      mutate(
        # targeting_req1 = if_else(delta_post>=0, 0, 1-((delta+fe)/(delta+fe-delta_post))),
        targeting_req1 = if_else(delta_post>=0, 0, 1-((delta_mean+fe_mean)/(delta_mean+fe_mean-delta_post))), ### REVISION
        targeting_req = if_else(targeting_req1>1, 1, targeting_req1)
        ) %>% 
      mutate(clade=paste0(stringr::str_to_title(clade), "s")) %>%
      # select(species, grp, clade, delta, delta_post, targeting_req, contains("sens")) 
      select(species, grp, clade, delta_mean, delta_post, targeting_req, contains("sens")) ### REVISION
    
    # Add median cost estimate
    df2 <- 
      summ_df %>%
      filter(key == df2_filter) %>% 
      mutate(cost = q50/100) %>%
      select(species, cost, contains("sens")) 
    
    df <- 
      left_join(df1, df2) %>% 
      ungroup() %>%
      mutate(species = factor(species)) %>%
      # mutate(species = fct_reorder(species, delta, .desc=TRUE)) 
      mutate(species = fct_reorder(species, delta_mean, .desc=TRUE)) ### REVISION
    
    df %>%
      ggplot(aes(y=species)) +
      geom_point(
        aes(x=delta_post, size=cost, col=targeting_req), 
        alpha=0.7
        ) +
      geom_point(
        data = filter(df, cost>=1),
        aes(x=delta_post, size=cost), 
        col="grey", 
        show.legend = F
        ) + 
      geom_point(
        # aes(x=delta), 
        aes(x=delta_mean),  ### REVISION
        size=3, stroke=1, shape=21, col="red"
        ) +
      geom_vline(xintercept = 0, lty=2) +
      geom_segment(
        # aes(yend=species, x=delta, xend=delta_post),
        aes(yend=species, x=delta_mean, xend=delta_post), ### REVISION
        arrow = arrow(length = unit(.25, "lines"))
        ) +
      scale_size_continuous(
        name=paste0("Cost (percent\nof ", toupper(goal), ")"),
        labels=percent, range=c(2,9)
        ) +
      scale_color_viridis(
        name=paste0("Req. targeting\nimprovement"),
        trans="reverse", direction=-1, 
        option="plasma", labels=percent
        ) +
      guides(
        size = guide_legend(order = 1),
        col = guide_colourbar(order = 2)
        ) +
      # coord_fixed(ratio = .025) +
      labs(x = expression(Rate~of~population~change~(Delta))) +
      facet_grid(clade~., scales = "free", space = "free", switch = "both") +
      theme(
        legend.title = element_text(),
        axis.title.y=element_blank(), 
        strip.placement = "outside" ## Alongside `switch="both"` in facet_grid() call above
        ) 
  }
