# DCAC and DB-SRA
# Version 4d

# time the whole enchilada
start.enchilada <- Sys.time()

# number of CPUs available for optional parallel processing
NumCPUs <- 2

# turn on parallel processing if #CPUs > 1
do.parallel <- ifelse(NumCPUs > 1, TRUE, FALSE)

# toggle M correction term in production model
M.correction <- TRUE

# specify number of random draws for DCAC and delay-difference model
Niter <- 10000

############################
# LOAD DATA AND R PACKAGES #
############################

# load snowfall package (also requires snow)
require(snowfall)

# get life history data
parms.df <- read.csv(paste(getwd(), "/input_files/parameter_file.csv", sep = ""), header = T)

# remove spp. with "exclude" flag
parms.df <- parms.df[parms.df$exclude < 1, ]

# vector of species codes
sp.vec <- as.character(parms.df$species.code)
N.spp <- length(sp.vec)

# import landings data
lands.df <- read.csv(paste(getwd(), "/input_files/landings.csv", sep = ""), header = T)

# get first and last year of landings for each species and append to parms.df
first.yr <- numeric(N.spp)
last.yr <- numeric(N.spp)
for (i in 1:N.spp)
{
  first.yr[i] <- min(lands.df[lands.df$sp.code == sp.vec[i], "year"])
  last.yr[i]  <- max(lands.df[lands.df$sp.code == sp.vec[i], "year"])
}

parms.df <- cbind.data.frame(parms.df, first.yr, last.yr)
dimnames(parms.df)[[1]] <- sp.vec
rm(i, first.yr, last.yr)

# print warning if DCAC start year or end year is beyond landings 
if (any(parms.df$DCAC.start.yr < parms.df$first.yr)) print("DCAC start year is earlier than landings")
if (any(parms.df$DCAC.end.yr > parms.df$last.yr)) print("DCAC end year is later than landings")

# print warning if any species has fewer catch records than years in time series
for (i in sp.vec)
{
  if ( nrow(lands.df[lands.df$sp.code == i, ]) != length(parms.df[i, "first.yr"]:parms.df[i, "last.yr"]))
  {
    print(paste("Species", i, "has missing catch records."))
  }
}

# sort annual landings data frame by species, year
lands.df <- lands.df[with(lands.df, order(sp.code, year)), ]

# convert time series of catch to list object (convenient format for DB-SRA)
Catch.list <- as.list(1:N.spp)
names(Catch.list) <- sp.vec
for (i in 1:N.spp)
{
  Catch.list[[i]] <- lands.df[lands.df$sp.code == sp.vec[i], "catch.mt"]
  names(Catch.list[[i]]) <- lands.df[lands.df$sp.code == sp.vec[i], "year"]
}

# sum catch for species i between start and target years for DCAC
DCAC.SumC <- numeric(N.spp)
names(DCAC.SumC) <- sp.vec
for (i in 1:N.spp)
{
  catch.i <- lands.df[lands.df$sp.code == sp.vec[i] &
                        lands.df$year >= parms.df[sp.vec[i], "DCAC.start.yr"] &
                        lands.df$year <= parms.df[sp.vec[i], "DCAC.end.yr"], "catch.mt"]
  DCAC.SumC[i] <- sum(catch.i)
}

parms.df <- cbind.data.frame(parms.df, DCAC.SumC)
DCAC.Nyrs <- parms.df$DCAC.end.yr - parms.df$DCAC.start.yr + 1
parms.df <- cbind.data.frame(parms.df, DCAC.Nyrs, AvgCatch = (DCAC.SumC / DCAC.Nyrs))
rm(DCAC.SumC, DCAC.Nyrs, i, catch.i)

###############################################################################
# function to simulate from bounded beta distribution with
# m=mean and s=standard deviation;
# NOTE: the mean is the mean of the truncated dist, but the SD is from a standard (0,1) beta
# the user specifies the upper and lower values [a,b]

rbeta.ab <- function(n, m, s, a, b)
{
  # calculate mean of corresponding standard beta dist
  mu.std <- (m - a) / (b - a)
  
  # calculate parameters of std. beta with mean=mu.std and sd=s
  alpha <- (mu.std ^ 2 - mu.std ^ 3 - mu.std * s ^ 2) / s ^ 2
  beta  <- (mu.std - 2 * mu.std ^ 2 + mu.std ^ 3 - s ^ 2 + mu.std * s ^ 2) / s ^ 2
  
  # generate n draws from standard beta
  b.std <- rbeta(n, alpha, beta)
  
  # linear transformation from beta(0,1) to beta(a,b)
  b.out <- (b - a) * b.std + a
  
  return(b.out)
}
###############################################################################

# initialize list object for random draws
rand.list <- as.list(1:N.spp)
names(rand.list) <- sp.vec

# generate random draws from input distributions
for (i in 1:N.spp)
{
  M <- parms.df[sp.vec[i], "M.est"]
  SD.lnM <- parms.df[sp.vec[i], "SD.lnM"]
  FMSYtoMratio <- parms.df[sp.vec[i], "FMSYtoMratio"]
  SD.FMSYtoMratio <- parms.df[sp.vec[i], "SD.FMSYtoMratio"]
  Delta <- parms.df[sp.vec[i], "Delta"]
  SD.Delta <- parms.df[sp.vec[i], "SD.Delta"]
  DeltaLowerBound <- parms.df[sp.vec[i], "DeltaLowerBound"]
  DeltaUpperBound <- parms.df[sp.vec[i], "DeltaUpperBound"]
  BMSYtoB0ratio <- parms.df[sp.vec[i], "BMSYtoB0ratio"]
  SD.BMSYtoB0ratio <- parms.df[sp.vec[i], "SD.BMSYtoB0ratio"]
  BMSYtoB0LowerBound <- parms.df[sp.vec[i], "BMSYtoB0LowerBound"]
  BMSYtoB0UpperBound <- parms.df[sp.vec[i], "BMSYtoB0UpperBound"]
  myseed <- parms.df[sp.vec[i], "random.seed"] + 9167 * (0:3)
  
  # lognormal distribution for natural mortality (M)
  set.seed(myseed[1])
  M.vec <- rlnorm(Niter, meanlog = (log(M) - 0.5 * SD.lnM ^ 2), sdlog = SD.lnM)
  
  # lognormal distribution for ratio of ratio FMSY/M
  set.seed(myseed[2])
  FMSYtoM.vec  <- rlnorm(Niter, meanlog = (log(FMSYtoMratio) - 0.5 * SD.FMSYtoMratio ^ 2), sdlog = SD.FMSYtoMratio)
  
  # simulate bounded beta draws for delta
  set.seed(myseed[3])
  Delta.vec <- rbeta.ab(Niter, Delta, SD.Delta,
                        DeltaLowerBound, DeltaUpperBound)
  
  # simulate draws for ratio of BMSY/B0
  set.seed(myseed[4])
  BMSYtoB0.vec  <- rbeta.ab(Niter, BMSYtoB0ratio, SD.BMSYtoB0ratio,
                            BMSYtoB0LowerBound, BMSYtoB0UpperBound)
  
  rand.list[[i]] <- cbind(M.vec, FMSYtoM.vec, Delta.vec, BMSYtoB0.vec)
  
}

# stop if RNGs generarate NAs (happens when delta SD is large relative to mean)
if(any(unlist(lapply(rand.list, FUN = function(x) any(is.na(x))))))
{
  stop("NAs in random draws -- check distribution parameters")
}

###############################################################################
# DCAC function that accepts vector of parameters from input distributions
# parm.vec is a vector of length 4 with values for M, Fmsy/M, Depletion Delta, and Bmsy/B0
# call this function from apply(); faster than using for loop
DCAC.fun <- function(parm.vec, SumOfCatch, NumberOfYears)
{
  
  # assign names to elements of parameter vector (parm.vec) to clarify code
  M        <- parm.vec[1]
  FMSYtoM  <- parm.vec[2]
  Delta    <- parm.vec[3] # proportion that biomass is reduced relative to biomass in start year (NOT "depletion")
  BMSYtoB0 <- parm.vec[4]
  
  # calculate DCAC
  DCAC <- SumOfCatch / (NumberOfYears + (Delta/(BMSYtoB0 * FMSYtoM * M)))
  
  return(DCAC)
}
###############################################################################

# PARALLEL EXECUTION OF DCAC

# initialize "cluster" (use of multiple cores)
sfInit(parallel = do.parallel, cpus = NumCPUs)

# export variables & functions needed by slaves for parallel operation
sfExport(list = c("DCAC.fun", "rand.list", "parms.df", "sp.vec"))

DCAC.start.time <- Sys.time() # record time to see how fast this works

# apply function DCAC.fun to elements of rand.list, parms.df$DCAC.SumC, and parms.df$DCAC.Nyrs
DCAC.list <- sfLapply(1:N.spp, function(i) apply(rand.list[[i]], MARGIN = 1, FUN = DCAC.fun,
                                                 SumOfCatch = parms.df[sp.vec[i], "DCAC.SumC"],
                                                 NumberOfYears = parms.df[sp.vec[i], "DCAC.Nyrs"])
)

# how long did DCAC calculation take?
print(Sys.time() - DCAC.start.time)

names(DCAC.list) <- sp.vec
DCAC.df <- as.data.frame(DCAC.list)
rm(DCAC.list)

# print summary (mean and quantiles) of DCAC distributions
DCAC.summary <- t(apply(DCAC.df, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
DCAC.summary <- cbind.data.frame(species = parms.df[, "common.name"],
                                 mean.DCAC = apply(DCAC.df, 2, mean),
                                 DCAC.summary)

# stop parallel cluster
sfStop()

# write DCAC output
write.csv(DCAC.df, "DCAC_results.csv", row.names = FALSE)
write.csv(DCAC.summary, "DCAC_summary.csv", row.names = TRUE)

####################################################################
# functions to get exponent in Fletcher parameterization of P-T model
# based on input value of BMSY/B0
get.n <- function(BMSYtoB0ratio)
{
  # define half-width of region around 1/e (Fox model) to randomly
  # assign a value that won't blow up but is reasonable approximation
  eps <- 1e-3
  
  # search for n when BMSY/B0 < B0/e
  if(BMSYtoB0ratio < (1 / exp(1) - eps))
  {
    # lower bound of 0.05 for n sets minimum BMSY/B0 just below 5% (about 4.3%); plenty low
    n.out <- optimize(phi.fun, interval = c(0.05, 1 - eps ^ 2), phi.true = BMSYtoB0ratio)[[1]]
  }
  
  # approximate Fox model for values of BMSY/B0 near B0/e
  if( all( BMSYtoB0ratio >= (1 / exp(1) - eps),
           BMSYtoB0ratio < (1 / exp(1) + eps)))
  {
    n.out <- sample(c(1 - eps, 1 + eps), size = 1)
  }
  
  # search for n when BMSY/B0 > B0/e
  if(BMSYtoB0ratio >= (1 / exp(1) + eps))
  {
    # upper bound of 90 for n allows BMSY to occur at up to approx. 95% of B0
    n.out <- optimize(phi.fun, interval = c(1 + eps ^ 2, 90), phi.true = BMSYtoB0ratio)[[1]]
  }
  
  return(n.out)
}

# objective function to minimize error between input BMSY/B0 (phi.true) and estimate based on current n value
phi.fun <- function(n, phi.true)
{
  phi.out <- (1 / n) ^ (1 / (n - 1))
  out <- (phi.true - phi.out) ^ 2
  return(out)
}

# estimate lower and upper bounds for B0 for each species
# lower bound = average catch used in DCAC
# upper bound = 110% of sumCatch/min(Delta) (or twice the hake B0, whichever is smaller)
B0.low.vec <- numeric(N.spp)
B0.high.vec <- numeric(N.spp)
for (i in 1:N.spp)
{
  B0.low.vec[i] <- parms.df[sp.vec[i], "AvgCatch"]
  B0.high.vec[i] <- min(3000000, 1.1 * sum(Catch.list[[i]]) / min(rand.list[[i]][, "Delta.vec"]))
}

# append B0 bounds to parms.df, so you can check if the optimization step hit the boundaries
parms.df <- cbind.data.frame(parms.df, B0.low = B0.low.vec, B0.high = B0.high.vec)
rm(B0.low.vec, B0.high.vec)

###############################################################################
# Pella-Tomlinson-Fletcher-Schaefer-MacCall model

Delay.Diff.fun <- function(i)
{
  first.yr  <- parms.df[sp.vec[i], "first.yr"]
  last.yr   <- parms.df[sp.vec[i], "last.yr"]
  
  # "delta.yr" IS THE YEAR THAT WILL BE FIT TO THE CURRENT DEPLETION VALUE ###
  delta.yr <- parms.df[sp.vec[i], "delta.yr"]
  
  catch.vec <- Catch.list[[i]]
  Amat <- parms.df[sp.vec[i], "age.mat"]
  N.yrs <- last.yr - first.yr + 1
  
  # assume B0 is greater than average catch
  B0.low  <- parms.df[sp.vec[i], "B0.low"]
  
  # assume B0 is less than it would be if surplus production were zero,
  # and depletion equaled (B0 - sumCatch)/B0; i.e. some catch is from production
  B0.high <- parms.df[sp.vec[i], "B0.high"]
  
  # create dataframe to hold results for this species:
  # 4 input distributions, FMSY, EMSY, BMSY, MSY,
  # OFL vector, biomass (B) vector, production (P) vector,
  # Bjoin, Flag.NegBiomass, Flag.MissTarget, Flag.OptimWarn, Flag.NAs
  temp.mat <- matrix(NA, ncol = (N.yrs + 1) * 3 + 8, nrow = Niter) # project time series one year beyond last landing
  results.df <- cbind.data.frame(rand.list[[i]], temp.mat)
  rm(temp.mat)
  
  # give names to columns in the results data frame
  names(results.df) <- c(names(results.df)[1:4],
                         c("FMSY", "EMSY", "BMSY", "MSY"),
                         paste("OFL", first.yr:(last.yr + 1), sep = ""),		# note extra year (forecast 1 yr past landings)
                         paste("B", first.yr:(last.yr + 1), sep = ""),		# note extra year (forecast 1 yr past landings)
                         paste("P", first.yr:(last.yr + 1), sep = ""),		# note extra year (forecast 1 yr past landings)
                         "Bjoin", "Flag.NegBiomass", "Flag.NAs","Flag.OptimWarn")
  
  results.df[, "Flag.OptimWarn"] <- 0
  results.df[, "Flag.NAs"] <- 0
  
  # calculate distribution of FMSY (instantaneous rate) from input values
  results.df[, "FMSY"] <- results.df[, "M.vec"] * results.df[, "FMSYtoM.vec"]
  
  # calculate distribution of EMSY (exploitation rate) from input values
  results.df[, "EMSY"] <- (1 - exp(-(results.df[, "M.vec"] + results.df[, "FMSY"]))) *
    ( results.df[, "FMSY"] / (results.df[, "M.vec"] + results.df[, "FMSY"]))
  
  for (j in 1:Niter)
  {  print(j)
    M        <- results.df[j, "M.vec"]
    FMSYtoM  <- results.df[j, "FMSYtoM.vec"]
    Delta    <- results.df[j, "Delta.vec"] # fraction biomass is reduced relative to biomass in start year (Delta is NOT "depletion")
    BMSYtoB0 <- results.df[j, "BMSYtoB0.vec"]
    FMSY     <- results.df[j, "FMSY"]
    EMSY     <- results.df[j, "EMSY"]
    
    if (any(Delta >= 1, Delta <= 0))
    {
      stop("Delta must be between 0 and 1 given assumptions of current production model")
    }
    Depletion <- 1 - Delta
    
    n <- get.n(BMSYtoB0)
    g <- (n ^ (n / (n - 1))) / (n - 1)
    
    # create vector to hold biomass time series for this set of parameters
    B.vec <- numeric(N.yrs + 1)
    
    # create a vector to hold production time series for this set of parameters
    P.vec <- numeric(N.yrs + 1)
    
    # production functions that minimize difference between obs/pred depletion in target year
    # over a range of B0 values; parm.vec is current draw of c(M, FMSYtoM, Delta, BMSYtoB0, plus FMSY, EMSY)
    prod.fun <- function(B0, parm.vec, C.vec, Amat, detailed.output = TRUE)
    {
      M        <- parm.vec[1]
      FMSYtoM  <- parm.vec[2]
      Delta    <- parm.vec[3]
      BMSYtoB0 <- parm.vec[4]
      FMSY     <- parm.vec[5]
      EMSY     <- parm.vec[6]
      
      # calculate BMSY for current guess at B0 and current value of BMSYtoB0
      BMSY <- B0 * BMSYtoB0
      MSY <- BMSY * EMSY
      
      # calculate Bjoin and "s" (slope of production/biomass ratio) if BMSYtoB0 <= 0.5
      if (BMSYtoB0 < 0.3)
      {
        Bjoin <- B0 * (0.5 * BMSYtoB0)
        s <- (1 - n) * g * MSY * (Bjoin ^ (n - 2)) * B0 ^ (-n)
      }
      if (all(BMSYtoB0 >= 0.3, BMSYtoB0 <= 0.5))
      {
        Bjoin <- B0 * (0.75 * BMSYtoB0 - 0.075)
        s <- (1 - n) * g * MSY * (Bjoin ^ (n - 2)) * B0 ^ (-n)
      }
      
      if (BMSYtoB0 > 0.5)
      {
        Bjoin <- NA
      }
      
      # production function for Pella-Tomlinson-Fletcher model
      PTF.fun <- function(n, g, MSY, Blag, B0)
      {
        if (Blag > 0)
        {
          P <- g * MSY * (Blag/B0) - g * MSY * (Blag / B0) ^ n
        } else {
          P <- 0
        }
        return(P)
      }
      
      # production function for hybrid PTF-Schaefer model
      # that approximates BHSRR productivity when BMSY/B0 < 0.5
      hybrid.fun <- function(n, g, MSY, Blag, B0, Bjoin, s)
      {
        if (Blag > 0)
        {
          if (Blag >= Bjoin)
          {
            P <- g * MSY * (Blag / B0) - g * MSY * (Blag / B0) ^ n
          } else {
            P <-  Blag * ( ((g * MSY * (Bjoin / B0) - g * MSY * (Bjoin / B0) ^ n) / Bjoin) + s * (Blag - Bjoin))
          }
        } else {
          P <- 0
        }
        return(P)
      }
      
      # assume B0 = B(t=1)
      B.vec[1] <- B0
      
      # production in first year is always zero (starting from B0)
      P.vec[1] <- 0
      
      # PROJECT FORWARD FROM CURRENT B0 USING PRODUCTION MODEL
      # two options: with and without M correction
      if (M.correction)
      {
        for (k in 2:(N.yrs + 1))
        {
          # production is zero (at B0) until source of spawning output has been harvested
          if (k <= Amat)
          {
            # remove production term (P) because production is based on unfished biomass until k>Amat,
            # but include M correction and replace B.vec[k-Amat] with B0
            B.vec[k] <- B.vec[k - 1] + (1 - exp(-M)) * (B0 - B.vec[k - 1]) - C.vec[k - 1]
          } else {
            
            # if BMSY/B0 > 0.5, then use P-T-F model
            if (BMSYtoB0 > 0.5)
            {
              P.vec[k] <- PTF.fun(n, g, MSY, B.vec[k - Amat], B0)
            } else {
              
              # if BMSY/B0 <= 0.5, use hybrid model
              P.vec[k] <- hybrid.fun(n, g, MSY, B.vec[k - Amat], B0, Bjoin, s)
            }
            # production model with M correction term
            B.vec[k] <- B.vec[k - 1] + P.vec[k] + (1 - exp(-M)) * (B.vec[k - Amat] - B.vec[k - 1]) - C.vec[k - 1]
          }
        }
      } else { # loop without M correction term
        for (k in 2:(N.yrs + 1))
        {
          # production is zero (at B0) until source of spawning output has been harvested
          if (k <= Amat)
          {
            # remove production term (P) because production is based on unfished biomass until k>Amat
            B.vec[k] <- B.vec[k - 1] - C.vec[k - 1]
          } else {
            
            # if BMSY/B0 > 0.5, then use P-T-F model
            if (BMSYtoB0 > 0.5)
            {
              P.vec[k] <- PTF.fun(n, g, MSY, B.vec[k - Amat], B0)
            } else {
              
              # if BMSY/B0 <= 0.5, use hybrid model
              P.vec[k] <- hybrid.fun(n, g, MSY, B.vec[k - Amat], B0, Bjoin, s)
            }
            # production model without M correction term
            B.vec[k] <- B.vec[k - 1] + P.vec[k] - C.vec[k - 1]
          }
        }
      }
      
      # add year labels to the time series vectors
      names(B.vec) <- names(P.vec) <- first.yr:(last.yr + 1)
      names(C.vec) <- first.yr:last.yr
      
      Btgt.to.B0.ratio <- as.numeric(B.vec[as.character(delta.yr)] / B0)
      obj.fun <- as.numeric(abs(B.vec[as.character(delta.yr)] - Depletion * B0))
      
      if (detailed.output == FALSE)
      {
        return(as.numeric(obj.fun))
      }
      if (detailed.output == TRUE)
      {
        return(list(Species = sp.vec[i],
                    ObjectiveFunction = as.numeric(obj.fun),
                    B0.Bounds = c(B0.low, B0.high),
                    # set catch in last.yr+1 to zero (no effect)
                    TimeSeries = cbind.data.frame(B.vec, P.vec, C.vec = c(C.vec, 0)), 
                    BMSY = BMSY,
                    MSY = MSY,
                    Bjoin = Bjoin))
      }
      
    } # end of prod.fun
    
    # global assignment (used in tryCatch)
    tag <<- 0
    
    tryCatch(B0.opt <- optimize(prod.fun,
                                interval = c(B0.low, B0.high),
                                ### the following are arguments passed to prod.fun()
                                parm.vec = as.numeric(results.df[j, 1:6]),
                                C.vec = catch.vec,
                                Amat = Amat,
                                detailed.output = FALSE,
                                tol = 0.0001),
             # if optimize() generates warning,
             # then tag this iteration and store result in "Flag.OptimWarn" column
             warning = function(warn) { tag <<- 1 },
             finally = B0.opt <- optimize(prod.fun,
                                          interval = c(B0.low, B0.high),
                                          ### the following are arguments passed to prod.fun()
                                          parm.vec = as.numeric(results.df[j, 1:6]),
                                          C.vec = catch.vec,
                                          Amat = Amat,
                                          detailed.output = FALSE,
                                          tol = 0.0001)
    )
    
    # generate model results as a list
    prod.out <- prod.fun(B0.opt[[1]],
                         parm.vec = as.numeric(results.df[j, 1:6]),
                         catch.vec,
                         Amat,
                         detailed.output = TRUE)
    
    # store model results for this iteration in this species
    Flag.NegBiomass <- any(prod.out[["TimeSeries"]]$B.vec <= 0)
    Flag.NAs <- any(is.na(prod.out[["TimeSeries"]]))
    result.vec <- as.numeric(c(prod.out[["BMSY"]],
                               prod.out[["MSY"]],
                               results.df[j,"EMSY"] * prod.out[["TimeSeries"]][, "B.vec"], # time series of OFL (EMSY*B(t))
                               prod.out[["TimeSeries"]]$B.vec, # time series of biomass
                               prod.out[["TimeSeries"]]$P.vec, # time series of production
                               prod.out[["Bjoin"]],
                               Flag.NegBiomass,
                               Flag.NAs,
                               tag)
    )
    results.df[j, 7:ncol(results.df)] <- result.vec
    
  } # end of iteration loop
  
  # add results.df to list containing results for all species
  return(results.df)
  
} # end of species loop (function)

### PARALLEL EXECUTION OF DB-SRA

# initialize "cluster" (use of multiple cores)
sfInit(parallel = do.parallel, cpus = NumCPUs)

# export variables & functions needed by slaves for parallel operation
sfExport(list = c("Catch.list",
                  "M.correction",
                  "rand.list",
                  "parms.df",
                  "sp.vec",
                  "Niter",
                  "get.n",
                  "phi.fun")
)

options(warn=1) # print any warnings as they occur -- only works in sequential mode (1 cpu)

DelayDiff.start.time <- Sys.time()

# Load-balanced parallel computation of all species (requires snowfall package)
results.list <- sfClusterApplyLB(1:N.spp, Delay.Diff.fun)

# how long did it take?
print(Sys.time() - DelayDiff.start.time)

# stop parallel cluster
sfStop()

names(results.list) <- sp.vec

N.yrs.vec <- parms.df$last.yr - parms.df$first.yr + 1		# need to redefine outside of delay diff fxn
names(N.yrs.vec) <- sp.vec

#### FLAG RUNS WITH OBVIOUS ERRORS ######

# create flag for runs in which final depletion doesn't match current draw
for (i in 1:N.spp)
{
  B.target <- results.list[[i]][, paste("B", parms.df[i, "delta.yr"], sep = "")]
  B0 <- results.list[[i]][, paste("B", parms.df[i, "first.yr"], sep = "")]
  Depletion <- 1 - results.list[[i]][, "Delta.vec"]
  Flag.Miss <- as.numeric(abs(B.target / B0 - Depletion) > 0.005) # depletion must be within 0.5% of target
  # if biomass in target year is NA, call it a miss (usually due to oscillations)
  Flag.Miss[is.na(Flag.Miss)] <- 1
  results.list[[i]] <- cbind(results.list[[i]], Flag.MissTarget = Flag.Miss)
  rm(B.target, B0, Depletion, Flag.Miss)
}

# create flag that identifies any runs that hit the bounds for B0
# DEFINITION: "hitting" the bound means coming within 1 mt of upper OR 1% of lower bound
for (i in 1:N.spp)
{
  B0.low.tol <- 0.01 * parms.df[i, "B0.low"] # B0 "hits" lower bound if it is within 1% of avg. catch
  B0.low.logical  <- abs(parms.df[i, "B0.low"] - results.list[[i]][, N.yrs.vec[i] + 10]) < B0.low.tol
  # upper bound is "hit" if B0 is within 1 ton
  B0.high.logical <- abs(parms.df[i, "B0.high"] - results.list[[i]][, N.yrs.vec[i] + 10]) < 1
  test.mat <- cbind(B0.low.logical, B0.high.logical)
  results.list[[i]] <- cbind(results.list[[i]], Flag.HitBounds = as.numeric(apply(test.mat, 1, any)))
}
rm(B0.low.logical, B0.high.logical, test.mat)

# create flag that identifies any runs with an error from optimize(), but no other errors,
# and also flag "good runs" for easy extraction into "results.good"
results.good <- as.list(1:N.spp)
for (i in 1:N.spp)
{
  # identify flagged runs from optimize()
  Optim.flagged  <- results.list[[i]][, "Flag.OptimWarn"]
  # identify runs with no other errors
  No.other.flags <- apply(results.list[[i]][, c("Flag.NegBiomass", "Flag.MissTarget", "Flag.NAs", "Flag.HitBounds")],
                          MARGIN = 1, sum) < 1
  Good.run.vec <- as.numeric(No.other.flags)
  ## replace NAs with 1
  Good.run.vec[is.na(Good.run.vec)] <- 1
  results.list[[i]] <- cbind(results.list[[i]],
                             Good.Run = Good.run.vec,
                             Flag.OptimWarnOnly = as.numeric(Optim.flagged > 0 & No.other.flags > 0))
  results.good[[i]] <- results.list[[i]][Good.run.vec > 0, ]
}
names(results.good) <- sp.vec
rm(Optim.flagged, No.other.flags, Good.run.vec)

# RECORD RESULTS OF ALL ITERATIONS FROM PTFS MODEL AS .CSV FILES
dir.create(paste(getwd(), "/model_output",sep = ""))
# write separate .csv files for each species' PTFS results
for (i in 1:N.spp)
{
  write.table(results.list[[i]], paste(getwd(), "/model_output/", sp.vec[i], "_all_results.csv", sep = ""),
              sep = ",", row.names = FALSE)
}

### save all trajectories that don't go negative, hit bounds, miss the target depletion level, or have NAs in timeseries
dir.create(paste(getwd(), "/retained_model_output", sep = ""))
for (i in 1:N.spp)
{
  write.table(results.good[[i]], paste(getwd(), "/retained_model_output/", sp.vec[i], "_good_results.csv", sep = ""),
              sep = ",", row.names = FALSE)
}

# summarize error messages (negative biomass, hitting bounds of B0, missing target, etc.)
errors.list <- list()
errors.df   <- data.frame(t(rep(NA, 9)))
names(errors.df) <- c("Species", "Niter", "GoodRuns", "Flag.NegBiomass", "Flag.MissTarget",
                      "Flag.HitBounds", "Flag.NAs", "Flag.OptimWarn", "Flag.OptimWarnOnly")
for(i in 1:N.spp)
{
  errors.list[[i]] <- results.list[[i]][, c("Good.Run", "Flag.NegBiomass", "Flag.MissTarget", "Flag.HitBounds",
                                            "Flag.NAs", "Flag.OptimWarn", "Flag.OptimWarnOnly")]
  errors.df[i, 1]   <- sp.vec[i]
  errors.df[i, 2]   <- Niter
  errors.df[i, 3:9] <- apply(errors.list[[i]], MARGIN = 2, sum)
}
rm(errors.list)
errors.df
write.table(errors.df, paste(getwd(), "/error_summary.csv", sep = ""), sep = ",", row.names = F, quote = F)

# save parms.df as a comma-delimited text file
write.table(parms.df, paste(getwd(), "/parms.csv", sep = ""), sep = ",", row.names = F, quote = F)



# CALC SUMMARY STATS FOR UNFISHED BIOMASS
# create list with B0 draws from each species
B0.yrs <- paste("B", parms.df[, "first.yr"], sep = "")
B0.list <- as.list(1:N.spp)
for (i in 1:N.spp)
{
  B0.list[[i]] <- results.good[[sp.vec[i]]][, B0.yrs[i]]
}
names(B0.list) <- sp.vec

x.tmp <- matrix(NA, nrow = N.spp, ncol = 6)
for (i in 1:N.spp)
{
  x.tmp[i, 1] <- mean(B0.list[[i]])
  x.tmp[i, 2:6] <- quantile(B0.list[[i]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
}

B0.stats <- cbind.data.frame(Common.Name = parms.df[, "common.name"], x.tmp)
names(B0.stats ) <- c("Common.Name", "B0.Mean", "2.5%", "25%", "50%", "75%", "97.5%")
print(B0.stats)
write.table(B0.stats, file = paste(getwd(), "/DB-SRA_B0_summary_stats.csv", sep = ""), sep = ",", row.names = F, quote = F)
rm(x.tmp)



# CALC SUMMARY STATS FOR OFL IN "DBSRA.OFL.yr" (doesn't have to be the same as "delta.yr")
# create list with OFL timeseries
OFL.list <- as.list(1:N.spp)
for (i in 1:N.spp)
{
  # use "grep" to identify column names that begin (\\b) with "OFL" (period after OFL is wildcard)
  OFL.list[[i]] <- results.good[[sp.vec[i]]][, grep("\\bOFL.", names(results.good[[sp.vec[i]]]))]
}
names(OFL.list) <- sp.vec

OFL.yrs <- paste("OFL", parms.df[, "DBSRA.OFL.yr"], sep = "")
x.tmp <- matrix(NA, nrow = N.spp, ncol = 6)
for (i in 1:N.spp)
{
  x.tmp[i, 1] <- mean(OFL.list[[i]][, OFL.yrs[i]])
  x.tmp[i, 2:6] <- quantile(OFL.list[[i]][, OFL.yrs[i]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
}

OFL.stats <- cbind.data.frame(Common.Name = parms.df[, "common.name"],
                              Delta.Year = parms.df[, "delta.yr"],
                              OFL.Year = parms.df[, "DBSRA.OFL.yr"],
                              x.tmp)
names(OFL.stats) <- c("Common.Name", "Delta.Year", "OFL.Year", "OFL.Mean", "2.5%", "25%", "50%", "75%", "97.5%")
print(OFL.stats)
write.table(OFL.stats, file = paste(getwd(), "/DB-SRA_OFL_summary_stats.csv", sep = ""), sep = ",", row.names = F, quote = F)
rm(x.tmp)

# time the whole enchilada
print("DB-SRA Total Time:")
total.run.time <- Sys.time() - start.enchilada
print(total.run.time)
