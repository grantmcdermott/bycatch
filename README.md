# Protecting marine mammals, turtles and birds by rebuilding global fisheries

This repo contains data and *R* code for reproducing [Burgess *et al*. (2017)](http://dx.doi.org/10.1126/science.aao4248), "Protecting marine mammals, turtles and birds by rebuilding global fisheries". 

> **Abstract:** Reductions in global fishing pressure are needed to end overfishing of target species and maximize the value of fisheries. We ask whether such reductions would also be sufficient to protect non–target species threatened as bycatch. We compare changes in fishing pressure needed to maximize profits from 4713 target fish stocks—accounting for >75% of global catch—to changes in fishing pressure needed to reverse ongoing declines of 20 marine mammal, sea turtle, and seabird populations threatened as bycatch. We project that maximizing fishery profits would halt or reverse declines of approximately half of these threatened populations. Recovering the other populations would require substantially greater effort reductions or targeting improvements. Improving commercial fishery management could thus yield important collateral benefits for threatened bycatch species globally.

## Reproducibility

Assuming that you met the [requirements](#requirements), the entire analysis can executed by running the **`R/bycatch.R`** master script. For example:

```
## Run these commands in the shell
git clone git@github.com:grantmcdermott/bycatch.git
cd bycatch
Rscript R/bycatch.R
```

Of course, you may run the analysis interactively in your R console too. The `R/bycatch.R` file contains self-explanatory code for easily reproducing the results from the ten different model runs described in the paper, as well as all of the figures.

## Requirements

#### Step 1: Install R and R libaries

**Note:** The code was most recently tested and updated against *R* 3.4.2. *R* is free, open-source and available for download [here](https://www.r-project.org/).

We use [**renv**](https://rstudio.github.io/renv/) to snapshot the project environment. Run the following command(s) from your *R* console to pull in all of the necessary libraries.

```r
# renv::init()   ## Only necessary if you didn't clone/open the repo as an RStudio project
renv::restore()  ## Enter "y" when prompted
```

#### Step 2: Install system dependencies (only if applicable)

While the `renv::restore()` command above should install [package binaries](https://packagemanager.rstudio.com/) on most operating systems (OS), it will not necessarily import *system* dependencies on some Linux builds. In our case, the most obvious system dependencies are related to the underlying geospatial libaries that power the [**sf**](https://r-spatial.github.io/sf/#installing) package. You can double check that you have met the requirements for your OS by clicking the previous link. As an alternative, you can also try running the [`remotes::system_requirements()`](https://remotes.r-lib.org/reference/system_requirements.html) command. For example, if we wanted to check what requirements were required for Ubuntu 20.04, we could run:

```r
remotes::system_requirements(os = 'ubuntu', os_release = '20.04', 
                             path = 'renv/library/R-4.0/x86_64-pc-linux-gnu/sf/')
```

#### Optional: Fonts

We use the `extrafont` package embed [Open Sans](https://fonts.google.com/specimen/Open+Sans) fonts in some of the figures. Please note that the Open Sans fonts would have to be installed separately on your system and thus requires some minor setup upon first use. See [here](https://github.com/wch/extrafont/blob/master/README.md) for instructions. Feel free to skip this step if that all sounds like too much work. The code will automatically use one *R*'s default Arial fonts if others are not available.

## Performance

The core analysis in this paper involves a series of computationally intensive Monte Carlo simulations. The code is optimized to run in parallel and will automatically exploit any multi-core capability on your machine. That being said, even on a 24 core Linux server, a single model run can take around 40 minutes to complete. (See [here](http://grantmcdermott.com/2017/05/30/rstudio-server-compute-engine/) for a tutorial on how to set up your own personal server using Google Compute Engine, complete with *R* and RStudio Server.) We strongly suggest reducing the "n1" and "n2" parameters (+/- line 75 in `R/bycatch.R`) before running the analysis on a local machine with few CPUs. Failing that, the saved results from our own model runs can be found in the `Results` sub-directory.

## Problems

If you have any trouble running the code, or find any errors, please file an issue on this repo and we'll look into it.
