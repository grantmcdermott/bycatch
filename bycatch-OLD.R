rm(list = ls())

library(readr)
library(ggplot2)
library(cowplot)
library(tibble)
library(tidyr)
library(dplyr)

##### 2 functions #####

### Function 1: Plot that asks: Will rebuilding be enough? 
### (i.e. distribution plot/frequency weight one) 

## Inputs: 
# 1: A dataframe with 4 columns: idorig, pctredfmsy, pctredfmey, weight (weight pre-normalized to add up to 1)
#    Sample of this dataframe attached as "sample-table-f1.csv"
# 2: 3 estimates of percent change needed for the bycatch species of interest: point/median, lower, upper

## Outputs: Plot with % Change in F on x-axis, and freq. distributions 
##  for each of msy and mey (i.e. pctredfmsy and pctredfmey), and vertical gridline at point/median estimate, 
##  with shaded region between lower and upper values given

df1 <- read_csv("sample-table-f1.csv")
## GM: Matt, when exporting CSV files via write.csv, remember to set row.names = F
## I've adjusted the sample datasets in this folder already, but just a tip for the future

df1 %>%
  rename(MSY = pctredfmsy,
         MEY = pctredfmey) %>%
  gather(key, reduction, MSY:MEY) %>%
  group_by(key, reduction) %>%
  summarise(culm_weight = sum(weight, na.rm = T)) %>%
  ggplot(aes(x = reduction, y = culm_weight, col = key)) +
  geom_point() +
  geom_line() +
  ## GM: We can test the geom_smooth function with bigger data...
  # geom_smooth(span = 0.2, se = F) +
  ## GM: You have to plug these next values in maually at the moment. 
  ## Ideally it would call an another dataframe automatically. This is
  ## easy to do, but I'd need to see the data to set up the function.
  annotate("rect", xmin = 50, xmax = 60, ymin = -Inf, ymax = Inf,
           alpha = .2) +
  geom_vline(xintercept = 55, lty = 2) +
  labs(x = "% Reduction in Effort", y = "Frequency") +
  theme(legend.title = element_blank())


### Function 2: How much does it cost to reduce mortality enough? (minimized) 

## Inputs: 
# 1: A dataframe with columns: idorig, pctredfmsy, pctredfmey, weight, fvfmsy, g, beta, phi, price, marginalcost
#    Sample of this dataframe attached as "sample-table-f2.csv"
# 2: 3 estimates of percent change needed for the bycatch species of interest: point/median, lower, upper
# 3: A desired percent chance that the bycatch reduction threshold is met
# 4: A function (below) describing the relationship for a target stock between the equilibrium profit 
#     and the percent reduction in F (pctredf; a variable) in terms of columns in the input table.
#     The function will be the same for all stocks.
eqprofit <- function(pctred,price,marginalcost,fvfmsy,g,k,phi,beta) {
  eqp <-  (0.01 * price * fvfmsy * g * k * (100 - pctred) * (((1 + phi - (fvfmsy * phi) + (0.01 * fvfmsy * pctred * phi))/(1 + phi))^(1/phi))) - 
    (marginalcost * ((fvfmsy * g * (1 - (pctred/100)))^beta)) + 
    return(eqp)
}

# 5: A function (below) describing the relationship for a target stock between the equilibrium yield (i.e. catch) 
#     and the percent reduction in F (pctredf; a variable) in terms of columns in the input table.
#     The function will be the same for all stocks.
eqyield <- function(pctred,fvfmsy,g,k,phi) {
  eqy <- 0.01 * fvfmsy * g * k * (100 - pctred) * (((1 + phi - (fvfmsy * phi) + (0.01 * fvfmsy * pctred * phi))/(1 + phi))^(1/phi))
  return(eqy)
}

# 6: A function (below) describing the relationship between bycatch benefit (in terms of percent reduction in bycatch mortality)
#     and the percent reduction in F for a target stock (pctredf; a variable)
bycatchbenefit <- function(pctred,weight) {
  benefit <- pctred * weight
  return(benefit)
}

## Outputs: 
# 1. 2 values (not a data frame): (i) minimum profit cost of achieving the point/median percent change needed for bycatch species
#                                 (ii) minimum yield cost of achieving the point/median percent change needed for bycatch species
# 2,3. Same as 1. but for each of the lower and upper percent change needed.
