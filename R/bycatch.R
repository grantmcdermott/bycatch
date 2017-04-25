## Clear environment
rm(list = ls())

#####################################
########### Load packages ###########
#####################################
library(data.table) ## Mostly for super fast reading/writing of large csv files
library(pbapply)
# library(pbmcapply) ## Now assign parallel computation through pbapply (with cl option)
library(parallel)
library(scales)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(extrafont) ## See https://github.com/wch/extrafont for first-time use instructions

#######################################################
########## Load functions and global elements #########
#######################################################

### Assign global elements for figures. 

## Assign font. Register to extrafont package DB first. If below font is not
## available, then extrafont package will use Arial default. Again, see: https://github.com/wch/extrafont
font_type <- choose_font(c("Open Sans", "sans")) ## Download here: https://fonts.google.com/specimen/Open+Sans
## Assign color scheme
bycatch_cols <- c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506",
                  "#a6cee3","#fb9a99","#984ea3","#ffff33")
## Make some adjustments to (now default) cowplot ggplot2 theme for figures
theme_update(
  text = element_text(family = font_type),
  legend.title = element_blank(),
  strip.background = element_rect(fill = "white"), ## Facet strip
  panel.spacing = unit(2, "lines") ## Increase gap between facet panels
)

## Load functions
source("R/bycatch_funcs.R")


##############################
########## Load data #########
##############################

### Load bycatch data
bycatch_df <- read_csv("Data/bycatch_species.csv")
target_df <- read_csv("Data/target_species.csv")

### Load target stock data, derived from the "upsides" model of Costello et al. 
### (PNAS, 2016).

## First choose which version of the upsides data to: 1) With uncertainty (have 
## to exclude NEI stocks), or 2) no uncertainty (can include NEI stocks). The 
## main results of the paper use the latter. The former are used for the 
## sensitivty analysis in the SI.
uncert_type <- c("uncert", "nouncert")[2] ## Change as needed.

## Now read in the data
upsides <- 
  fread(paste0("Data/upsides_", uncert_type, ".csv")) %>% 
  as_data_frame()


################################
########### Analysis ###########
################################

### MCMC sampling parameters

## How many draws (states of the world) are we simulating for each species?
n1 <- 10000
## How many times do we sample (with replacement) over target stocks to resolve 
## uncertainty for a single draw?
n2 <- 100 
## Set "alpha" parameter, i.e. elasticity of (changes in) bycatch to (changes 
## in) target stocks. Default is 1. Other options used in sensitivity analysis.
## See equation (S14). 
alpha_exp <- c(1, 0.5, 2)[1] 

### Results for all species

## First get a vector of bycatch species
all_species <- bycatch_df$species 

## Now apply the MCMC simulation function over all species (and bind into a 
## common data frame). A progress bar (PB) will give you an indication of how  
## long you have to wait. For the non-parallel version, you'll see one PB per
## species, i.e. 20 in total. For the parallel version, you'll see a single PB
## updated in clumps, i.e. corresponding to how many CPUs you have.
# all_dt <- pblapply(all_species, bycatch_func) %>% bind_rows() ## Non-parallel version (slower)
all_dt <- pblapply(all_species, bycatch_func, cl = detectCores()) %>% bind_rows() ## Parallel version (faster)

## Write results for convenient later use
write_csv(all_dt, paste0("Results/bycatch_results_", uncert_type, ".csv"))

# rm(upsides)

################################################
### READ IN THE PREVIOUSLY GENERATED RESULTS ###
################################################

all_dt <- read_csv(paste0("Results/bycatch_results_", uncert_type, ".csv"))

alldistplots <- 
  bycatchdistggplot(all_dt) +
  facet_wrap(~species, ncol = 3, scales = "free_x") 
alldistplots + ggsave(paste0("Figures/dist-", uncert_type,".png"), width = 10, height = 13)
alldistplots + ggsave(paste0("Figures/PDFs/dist-", uncert_type,".pdf"), width = 10, height = 13)

allcostplots <- 
  costggplot(all_dt) +
  facet_wrap(~species, ncol = 3, scales = "free_x")
allcostplots + ggsave(paste0("Figures/cost-", uncert_type,".png"), width = 10, height = 13)
allcostplots + ggsave(paste0("Figures/PDFs/cost-", uncert_type,".pdf"), width = 10, height = 13)

rm(alldistplots, allcostplots)
dev.off()

####################################################
########### Summary Figs and sensitivity ###########
####################################################

### Get median, 2.5th, 25th, 75th, 97.5th percentiles from main run for Figs. 3a,b

## GRM: ASSIGN NAME 
# results_summary_nouncert <- resultssummary(all_dt) ## GRM: DELETE OLD FUNCTION
results_summary_nouncert <- summ_func(all_dt)
# rm(all_dt)


#########################################
###########  Sensitivity runs ########### 
#########################################

## Uncertainty, alpha = 1
all_dt_uncert <- read_csv("Results/bycatch_results_uncert.csv")
# results_summary_uncert <- resultssummary(all_dt_uncert)
results_summary_uncert <- summ_func(all_dt_uncert)
rm(all_dt_uncert)

## No uncertainty, alpha = 0.5
all_dt_alpha05 <- read_csv("Results/bycatch_results_nouncert_alpha=05.csv")
# results_summary_alpha05 <- resultssummary(all_dt_alpha05)
results_summary_alpha05 <- summ_func(all_dt_alpha05)
rm(all_dt_alpha05)

## No uncertainty, alpha = 2
all_dt_alpha2 <- read_csv("Results/bycatch_results_nouncert_alpha=05.csv")
# results_summary_alpha2 <- resultssummary(all_dt_alpha2)
results_summary_alpha2 <- summ_func(all_dt_alpha2)
rm(all_dt_alpha2)


### Sensitivity for Supplement, where FbMEY/Fbcurrent = (FtargetMSY/Ftargetcurrent) ^ b

### Part 1: show how this assumption drives a result for a hypothetical species
### (i.e. in general)

## Fig. S4
fig_s4 <-
  ggplot(data_frame(x = c(0, 5)), aes(x = x)) +
  stat_function(fun = function(x) (1 - (0.5^x))) +
  geom_vline(aes(xintercept = 1), lty = 2) +
  geom_hline(aes(yintercept = 0.5), lty = 2) +
  scale_y_continuous(label = percent) +
  annotate("text", label = paste(expression(alpha==1)), x = 1.5, y = 0.03, parse=T, family = font_type) +
  annotate("text", label = "        Reduction in target \nspecies mortality (50%)", x = 3.5, y = 0.75, family = font_type) +
  labs(
    x = expression(alpha), 
    y = "Reduction in bycatch mortality"
    ) 
fig_s4 + ggsave("Figures/Fig-S4.png", width = 4, height = 4)
fig_s4 + ggsave("Figures/PDFs/Fig-S4.pdf", width = 4, height = 4)
rm(fig_s4)

### Part 2: Show how it affects our results

## Fig. 3
fig_3mey <- tradeoffs_plot(results_summary_nouncert, "MEY")
fig_3mey + ggsave("Figures/Fig-3-mey.png", width=10*.6, height=13*.6)
fig_3mey + ggsave("Figures/PDFs/Fig-3-mey.pdf", width=10*.6, height=13*.6)
dev.off()

fig_3msy <- tradeoffs_plot(results_summary_nouncert, "MSY")
fig_3msy + ggsave("Figures/Fig-3-msy.png", width=10*.6, height=13*.6)
fig_3msy + ggsave("Figures/PDFs/Fig-3-msy.pdf", width=10*.6, height=13*.6)
dev.off()

# ## MEY
# save_plot("Figures/Fig-3-mey.pdf", fig3mey(results_summary_nouncert, -50),
#           ncol = 1, # we're saving a grid plot of 2 columns
#           nrow = 2, # and 2 rows
#           # each individual subplot should have an aspect ratio of 1.3
#           base_aspect_ratio = 1.3
#           )
# ## MSY
# save_plot("Figures/Fig-3-msy.pdf", fig3msy(results_summary_nouncert, -75),
#           ncol = 1, # we're saving a grid plot of 2 columns
#           nrow = 2, # and 2 rows
#           # each individual subplot should have an aspect ratio of 1.3
#           base_aspect_ratio = 1.3
#           )


# Joint figure
fig3all(results_summary_no_uncert,-75)

# Sensitivity figure
sensfigmsy(results_summary_no_uncert, 
           results_summary_uncert, 
           results_summary_alpha05, 
           results_summary_alpha2, 
           -75)

sensfigmey(results_summary_no_uncert, 
           results_summary_uncert, 
           results_summary_alpha05, 
           results_summary_alpha2, 
           -75)


###################################################
########### Fig. 2 - Loggerhead Example ###########
###################################################

sp_type <- "Loggerhead turtle (NW Atlantic)"

### Fig 2.A

library(png)
library(gridGraphics)
# library(viridis)

img1 <- readPNG("Figures/AnimalSilhouettes/turtle-silhouette.png")
g1 <- rasterGrob(img1, interpolate=FALSE)

lh_delta <- (filter(bycatch_df, species == sp_type))$delta
lh_fe <- (filter(bycatch_df, species == sp_type))$fe

fig2a <-
  crossing(delta = seq(-5,0, length.out = 100), fe = seq(0,10, length.out = 100)) %>%
  mutate(z = abs(delta/fe)*100) %>%
  mutate(z = ifelse(z>100, 100, z)) %>%
  mutate_all(funs(./100)) %>%
  ggplot(aes(delta, fe, fill = z)) + 
  geom_raster(interpolate = T) +
  scale_fill_gradientn(
    name = expression(Delta/~italic(F)[e]),#bquote(atop("Reduction in"~italic(F)[e], "to halt decline ")), 
    colours = c("#F2F2F2FF", brewer_pal(palette = "Spectral")(11)), 
    trans = "reverse",
    labels = percent
    ) +
  annotation_custom(g1, xmin=lh_delta-0.005, xmax=lh_delta+0.005, ymin=lh_fe-0.005, ymax=lh_fe+0.005) +
  # geom_vline(xintercept = lh_delta*c(.75, 1.25), lty = 2) +
  geom_segment(aes(y=-Inf, yend=.1, x=lh_delta*.75, xend=lh_delta*.75), lty=2) +
  geom_segment(aes(y=-Inf, yend=.1, x=lh_delta*1.25, xend=lh_delta*1.25), lty=2) +
  # geom_hline(yintercept = lh_fe*c(.75, 1.25), lty = 2) +
  geom_segment(aes(x=-Inf, xend=0, y=lh_fe*.75, yend=lh_fe*.75), lty=2) +
  geom_segment(aes(x=-Inf, xend=0, y=lh_fe*1.25, yend=lh_fe*1.25), lty=2) +
  labs(
    x = expression(Rate~of~population~decline~(Delta)),
    y = expression(Bycatch~mortality~rate~(italic(F)[e]))
  ) +
  theme(legend.title = element_text())
# detach(package:gridGraphics)
# detach(package:png)


### Fig 2.B

### Fig 2.C
## First select the relevant stocks
fig2c <- 
  stockselect_func(sp_type) %>%
  samples_plot()
# ## Then plot the results figure (Note: Log scale)
# fig2c +
#   ggsave("Figures/Fig-2c.png", height = 3, width = 9)
# fig2c +
#   ggsave("Figures/PDFs/Fig-2c.pdf", height = 3, width = 9)

### Fig 2.D
fig2d <- bycatchdistggplot(all_dt %>% filter(species==sp_type)) 

### Fig 2. E
fig2e <- costggplot(all_dt %>% filter(species==sp_type))

### Extract Legend (to serve as common legend at bottom of composite figure) 
g_legend <- 
  function(a_ggplot){ 
    tmp <- ggplot_gtable(ggplot_build(a_ggplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)
    } 
legend <- g_legend(fig2d) 
# grid::grid.draw(legend)

### Tweak plots before putting theme together in composite figure
fig2d <- fig2d + theme(strip.text.x = element_blank(), legend.position = "none")
fig2e <- fig2e + theme(strip.text.x = element_blank(), legend.position = "none")

### Finally, draw composite figure
# fig2 <-
#   ggdraw() +
#   # draw_plot(figureName, xpos, ypos, width, height) +
#   draw_plot(fig2a, 0, 0.7, 0.5, 0.3) +
#   draw_plot(fig2a, 0.5, 0.7, 0.5, 0.3) +
#   draw_plot(fig2c, 0, 0.375, 1, 0.3) +
#   draw_plot(fig2d, 0, 0.05, 0.5, 0.3) +
#   draw_plot(fig2e, 0.5, 0.05, 0.5, 0.3) +
#   draw_plot(legend, 0, 0, 1, 0.05) +
#   draw_plot_label(c("A", "B", "C", "D", "E", ""), c(0, 0.5, 0, 0, 0.5, 0), c(1, 1, 0.675, 0.35, 0.35, 0), size = 15) 
fig2 <-
  ggdraw() +
  # draw_plot(figureName, xpos, ypos, width, height) +
  draw_plot(fig2a, 0, 0.7, 0.475, 0.3) +
  draw_plot(fig2a, 0.525, 0.7, 0.475, 0.3) +
  draw_plot(fig2c, 0, 0.375, 1, 0.3) +
  draw_plot(fig2d, 0, 0.05, 0.475, 0.3) +
  draw_plot(fig2e, 0.525, 0.05, 0.475, 0.3) +
  draw_plot(legend, 0, 0, 1, 0.05) +
  draw_plot_label(c("A", "B", "C", "D", "E", ""), c(0, 0.525, 0, 0, 0.525, 0), c(1, 1, 0.675, 0.35, 0.35, 0), size = 15)

fig2

save_plot("Figures/fig2.png", fig2,
          base_height = 13,
          base_aspect_ratio = 1/1.3
          )
save_plot("Figures/PDFs/fig2.pdf", fig2,
          base_height = 13,
          base_aspect_ratio = 1/1.3
          )


### WORKING HERE ###

# plot.iris <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) + 
#   geom_point() + facet_grid(. ~ Species) + 
#   panel_border() # and a border around each panel

ggdraw() +
  draw_plot(fig2c, 0, .5, 1, .5) +
  draw_plot(fig2d, 0, 0, .5, .5) +
  draw_plot(fig2e, .5, 0, .5, .5) +
  draw_plot_label(c("C", "D", "E"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)

################################################################
########### Fig 1b. Overall Reductions (MEY and MSY) ###########
################################################################

ovrred <- 
  upsides %>%
  group_by(idoriglumped) %>%
  summarise(margc = mean(marginalcost),
            bet = mean(beta),
            g = mean(g),
            fvfmey = mean(eqfvfmey),
            fvfmsy = mean(fvfmsy),
            pctmey = mean(pctredfmey),
            pctmsy = mean(pctredfmsy)) %>%
  mutate(wt = margc * ((g * fvfmsy)^bet),
         cstcurr = wt,
         cstmey = margc * (((g * fvfmsy)/fvfmey)^bet),
         cstmsy = margc * ((g)^bet)) %>%
  ungroup() %>%
  mutate(wt = wt/sum(wt, na.rm = T)) %>%
  mutate(wtpctmey = wt * pctmey,
         wtpctmsy = wt * pctmsy) 

avpctmey <- sum(ovrred$wtpctmey, na.rm = T)
avpctmsy <- sum(ovrred$wtpctmsy, na.rm = T)

avpctmey <- 100 * (1 - (sum(ovrred$cstmey, na.rm = T)/sum(ovrred$cstcurr, na.rm = T)))
avpctmsy <- 100 * (1 - (sum(ovrred$cstmsy, na.rm = T)/sum(ovrred$cstcurr, na.rm = T)))

## List of species categories
# list("Shads" = 24, "Flounders, halibuts, soles" = 31, 
#   "Cods, hakes, haddocks" = 32,"Miscellaneous coastal fishes" = 33,
#  "Miscellaneous demersal fishes" = 34,"Herrings, sardines,anchovies" = 35,
# "Tunas,bonitos,billfishes" = 36,"Miscellaneous pelagic fishes" = 37,
#"Sharks, rays, chimeras" = 38,"Shrimps, prawns" = 45,
#"Carps, barbels and other cyprinids" = 11,"Sturgeons, paddlefishes" = 21,
#"Salmons, trouts, smelts" = 23,"Miscellaneous diadromous fishes" = 25,
#"Crabs, sea-spiders" = 42,"Lobsters, spiny rock lobsers" = 43,
#"King crabs, squat lobsters" = 44,"Miscellaneous marine crustaceans" = 47,
#"Abalones, winkles, conchs" = 52,"Oysters" = 53,"Mussels" = 54,
#"Scallops, pectens" = 55,"Clams, cockles, arkshells" = 56,
#"Squids, cuttlefishes, octopuses" = 57,"Horseshoe crabs and other arachnoids" = 75,
#"Sea-urchins and other echinoderms" = 76,"Miscellaneous aquatic invertebrates" = 56)
