library(dplyr)
library(readr)

bycatch_df <-
  data_frame(
    species = c("Loggerhead turtle", "Leatherback turtle", "Olive ridley turtle (NEI)", "Olive ridley turtle (WI)", 
                "Australian sea lion", "NZ sea lion", "Finless porpoise", "Humpback dolphin", "Vaquita",
                "Amsterdam albatross", "Sooty shearwater", "Tristan albatross", "White-chinned petrel"),
    region = c("NW Atlantic", "East Pacific", "NE Indian Ocean", "W Indian Ocean", 
               "Australia", "New Zealand", "NW Pacific", "SW Pacific", "Gulf of California - Mexico",
               "S Indian Ocean", "Atlantic, Pacific, S Indian Ocean, S Ocean", "S Atlantic", "S Pacific, S Atlantic, S Indian, S Ocean"),
    grp = c("turtle", "turtle", "turtle", "turtle", 
            "mammal", "mammal", "mammal", "mammal", "mammal",
            "bird", "bird", "bird", "bird"), 
    req_red = c(NA, 0.83, 0.23, 1.29, 
                0.065, 0.59, 0.57, 0.38, 0.76, 
                0.53, 6.00, 0.91, 0.28),
    lambda = c(0.98, 0.846, 0.97, 0.82, 
               0.9923, 0.96, 0.947, 0.9754, 0.87, 
               0.953, 0.982, 0.971, 0.964),
    fe = c(0.067, 0.334, 0.0871, 0.14, 
           0.1177, 0.06734, 0.093, 0.0646, 0.17, 
           0.088, 0.003, 0.032, 0.13),
    lambda2 = c(NA, 0.97, NA, NA, 
                NA, NA, NA, NA, NA, 
                NA, NA, NA, NA),
    fe2 = c(NA, 0.036, NA, NA, 
            NA, NA, NA, NA, NA, 
            NA, NA, NA, NA),
    ref = c(NA, "Kaplan et al. (2005)", NA, NA, 
            NA, NA, NA, NA, NA,
            NA, NA, NA, NA)
  ) 

bycatch_df <-
  bycatch_df %>%
  mutate(pctredbpt = 100 * (1 - lambda)/fe, # Point estimate
         pctredbl = 100 * ((0.75 * (1 - lambda))/(fe * 1.25)), # Lower bound: Assumes decline 25% smaller, fe 25% larger
         pctredbu = 100 * ((1.25 * (1 - lambda))/(fe * 0.75)) # Upper bound: Assumes decline 25% larger, fe 25% smaller
  )

write_csv(bycatch_df, "Data/bycatch_species.csv")
