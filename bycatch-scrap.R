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
