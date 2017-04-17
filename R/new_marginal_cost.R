### Update the marginal cost parameters for the input files

# Libraries
library(readr)
library(dplyr)
library(tidyr)

# Set working directory to where input data is (files too large for Github)
setwd("\\\\babylon/gaineslab/bowashi/SFG/Threat Potential in Fisheries/Model_Data")

# Read in input data
data <- read_csv("bycatch-upuncert-input.csv")

# Attempt at calculating marginal cost
# pmin is parallel version of min
BOA <- pmin(1.99, data$bvbmsy * rlnorm(1, 0, ErrorSize)) # What is 1.99 from and what should ErrorSize be? 
                                                          # bvbmsy from data isn't open access either is it?
c_num <- data$price * (2 - BOA) * BOA * data$msy * 2 ^ data$beta
c_den <- ((2 - BOA) * r) ^ data$beta  # Don't have an r value (use g?)
cost <- c_num / c_den


##################################################################################
# MARGINAL COST CODE FROM TYLER
BOA <- pmin(1.99, RecentStockData$BvBmsyOpenAccess[1] * rlnorm(1, 0, ErrorSize))

c_num <-  Price * (2 - BOA) * BOA * MSY * 2 ^ beta

c_den <- ((2 - BOA) * r) ^ beta

cost <- c_num / c_den
##################################################################################