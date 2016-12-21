### Update the marginal cost parameters for the input files

# Libraries
library(readr)
library(dplyr)
library(tidyr)

# Set working directory to where input data is (files too large for Github)
setwd("\\\\babylon/gaineslab/bowashi/SFG/Threat Potential in Fisheries/Model_Data")

# Read in input data
data <- read_csv("bycatch-upuncert-input.csv")




##################################################################################
# MARGINAL COST CODE FROM TYLER
BOA <- pmin(1.99,RecentStockData$BvBmsyOpenAccess[1] *rlnorm(1,0,ErrorSize))

c_num <-  Price*(2-BOA)*BOA*MSY*2^beta

c_den = ((2-BOA)*r)^beta

cost = c_num/c_den
##################################################################################