# Protecting marine mammals, turtles, and birds by rebuilding global fisheries

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1188538.svg)](https://doi.org/10.5281/zenodo.1188538)

This repo contains data and *R* code for reproducing Burgess *et al*. (2018), "Protecting marine mammals, turtles, and birds by rebuilding global fisheries". Click on the "fork" button at the very top right of the page to create an independent copy of the repo within your own GitHub account. Alternately, click on the green "clone or download" button just below that to download the repo to your local computer.

The main file for conducting the analysis is `R/bycatch.R`. This file contains self-explanatory code for easily reproducing the results from the ten different model runs described in the paper, as well as all of the figures.

## Requirements

The entire analysis is conducted in the *R* programming environment. *R* is free, open-source and available for download [here](https://www.r-project.org/). We highly recommend running *R* in the RStudio IDE, which you can also download for free [here](https://www.rstudio.com/products/rstudio/download/).

You will need to install a number of external *R* packages to run the code successfully. These are listed at the top of the main `R/bycatch.R` file. An easy way to ensure that you have the correct versions of all the packages is to run the following code chunk in your *R* console:

```
if (!require("pacman")) install.packages("pacman")
pacman::p_install(data.table, pbapply, parallel, R.utils, truncnorm, scales, grid, rworldmap, sf, rgeos, tidyverse, forcats, cowplot, ggthemes, RColorBrewer, viridis, extrafont, here)
devtools::install_github("tidyverse/ggplot2") ## Using dev. version of ggplot2 for geom_sf()
pacman::p_update()
```

Please note that the `extrafont` package requires some minor setup upon first use. The main analysis will still work without this initial setup, but it is perhaps useful for reproducing some of the figures in full. See [here](https://github.com/wch/extrafont/blob/master/README.md) for instructions.

## Performance

The core analysis in this paper involves a series of computationally intensive Monte Carlo simulations. The code is optimized to run in parallel and will automatically exploit any multi-core capability on your machine. That being said, even on a 24 core Linux server, a single model run can take around 40 minutes to complete. (See [here](http://grantmcdermott.com/2017/05/30/rstudio-server-compute-engine/) for a tutorial on how to set up your own personal server using Google Compute Engine, complete with *R* and RStudio Server.) We strongly suggest reducing the "n1" and "n2" parameters (+/- line 75 in `R/bycatch.R`) before running the analysis on a local machine with few CPUs. Failing that, the saved results from our own model runs can be found in the `Results` sub-directory.

## Problems

If you have any trouble running the code, or find any errors, please file an issue on this repo and we'll look into it.

## License

The software code contained within this repository is made available under the [MIT license](http://opensource.org/licenses/mit-license.php). The data and figures are made available under the [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/) license.
