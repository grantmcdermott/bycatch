## Clear environment
rm(list = ls())

#####################################
########### Load packages ###########
#####################################
library(data.table) ## Mostly for super fast reading/writing of large csv files
library(pbapply)
library(parallel)
library(R.utils)
library(truncnorm)
library(scales)
library(grid)
# library(maps)
library(rworldmap) ## Better shape files
library(sf)
library(rgeos)
library(tidyverse) ## NOTE: Using dev. version of ggplot2 for geom_sf() devtools::install_github("tidyverse/ggplot2")
library(forcats)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(extrafont) ## See https://github.com/wch/extrafont for first-time use instructions


#######################################################
########## LOAD FUNCTIONS AND GLOBAL ELEMENTS #########
#######################################################

set.seed(123) 

### Assign global elements for figures. 

## Assign font. Register to extrafont package DB first. If below font is not
## available, then extrafont package will use Arial default. Again, see: https://github.com/wch/extrafont
font_type <- choose_font(c("Open Sans", "sans")) ## Download here: https://fonts.google.com/specimen/Open+Sans

## Make some adjustments to the (now default) cowplot ggplot2 theme for figures
## to match publication requirements. See: 
## http://www.sciencemag.org/site/feature/contribinfo/prep/prep_revfigs.xhtml
theme_science <-
  theme_cowplot(font_size = 7, line_size = 0.15) +
  theme(
    text = element_text(family = font_type),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "lines"),
    legend.margin = margin(t=3, r=0, b=3, l=0, unit="pt"), ## Change last digit to -ve to reduce space between right legend and plot
    legend.spacing = unit(0.2, "cm"),
    legend.justification = "center", 
    plot.margin = margin(t=3, r=3, b=3, l=3, unit="pt"),#unit(c(3, 3, 3, 3), "points"), 
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.spacing = unit(1, "lines"), ## Increase gap between facet panels
    strip.text = element_text(size = 7)
    )
theme_set(theme_science)
fig_width <- 2.3 ## Width of one column in Science (in inches)
## Also adjust geoms (https://stackoverflow.com/a/21175042)
params <- ls(pattern = '^geom_', env = as.environment('package:ggplot2'))
params <- setdiff(params, paste0("geom_", c("bin2d", "count", "freqpoly", "histogram", "jitter", 
                                            "qq", "qq_line")))
geoms <- gsub("geom_", "", params)
lapply(geoms, update_geom_defaults, list(size = 1.3, stroke = 0.225))
lapply(c("line", "vline", "hline", "segment", "density"), update_geom_defaults, list(size = 0.2))
rm(params, geoms)

## Decide on number of cores for parallel computation
num_cores <- detectCores() ## i.e. Use all available CPUs. Subtract 1 or 2 if you are running additional processes on your computer

## Load functions
source("R/bycatch_funcs.R")


##############################
########## LOAD DATA #########
##############################

### Load bycatch data
bycatch_df <- read_csv("Data/bycatch_species.csv")
target_df <- read_csv("Data/target_species.csv")

## Get a vector of bycatch species
all_species <- bycatch_df$species 


################################
########### ANALYSIS ###########
################################

#### WARNING: FULL ANALYSIS CAN TAKE A LONG TIME TO RUN (+/- 40 MIN ON A
#### 24 CORE LINUX SERVER). SKIP DIRECTLY TO FIGURES SECTION (LINE 130)    
#### TO PLOT PREVIOUSLY RUN (AND SAVED) RESULTS. 

### Monte Carlo (MC) sampling parameters

## How many draws (states of the world) are we simulating for each species?
n1 <- 1000
## How many times do we sample (with replacement) over target stocks to resolve 
## uncertainty for a single draw?
n2 <- 100 

### Choose model run. The `choose_run` function below sets the correct parameters  
### for each run and also loads the correct version of the upsides data into the 
### global environment. It takes as argument one of 10 shorthand model descriptions
### listed in the `run` vector. Alternatively, you may simply enter the integer between 
### 1 and 10 that corresponds to the run of your choice.
run <- 
  c("main", "fcorrected", "conservation", "alpha=05", "alpha=2", 
    "nonei", "weights", "2012only", "doubleuncert", "kitchen")[1] ## Change as needed

choose_run(run) ## choose_run(1) works equally as well 

### Results for all species

## Apply the MC simulation function over all species (and bind into a common
## data frame). A series of progress bars will give you an indication of how
## long you have to wait, with one progress shown per bycatch species. In other 
## words, you'll see 20 progress bars in total if you run the full sample.
results <- lapply(all_species, bycatch_func) %>% bind_rows() 

## Get summary results
results_summary <- summ_func(results)

## Write results to disk for convenient later use
write_csv(results, paste0("Results/bycatch_results", suff_str, ".csv"))
write_csv(results_summary, paste0("Results/bycatch_summary_results", suff_str, ".csv"))

####################################
########### MAIN FIGURES ########### 
####################################

## First, read the main results back in (no uncertainty, alpha = 1)
results <- read_csv("Results/bycatch_results.csv")
results_summary <- read_csv("Results/bycatch_summary_results.csv")

## Choose map projection (See http://spatialreference.org)
proj_string <- 
  c("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", ## Default (Mercator?)
    "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", ## Robinson World
    "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs" ## Google projection
  )[2]

## Load basic world/countries spatial data for maps
countries <-
  st_as_sf(countriesLow) %>% ## data(countriesLow) is from the rworldmap package
  st_transform(proj_string)

#####################################
### Fig. 1 (Threats and upsides) ####
#####################################

sp_type <- all_species ## All species

#### Fig 1.A (Heatmap: Bycatch mortality VS. Population decline) #####

fig1a <-
  crossing(
    delta_mean = seq(-45, 0, length.out = 100), 
    fe_mean = seq(0, 45, length.out = 100)
    ) %>%
  mutate(pctT = abs(delta_mean/fe_mean)*100) %>%
  mutate(pctT = ifelse(pctT>100, 100, pctT)) %>%
  mutate_all(funs(./100)) %>%
  ggplot(aes(x = delta_mean, y = fe_mean)) + 
  geom_raster(aes(fill = pctT), interpolate = T) +
  scale_fill_gradientn(
    name = expression('%'~italic(T)),
    # colours = brewer_pal(palette = "Spectral")(11), 
    colours = rev(brewer_pal(palette = "YlOrRd")(9)), 
    trans = "reverse",
    labels = percent
    ) +
  geom_polygon(
    data=data_frame(delta_mean=c(-45,-45,0)/100, fe_mean=c(0,45,0)/100), 
    fill="#F2F2F2FF", col="#F2F2F2FF", 
    lwd = 0.6
    ) +
  labs(
    x = expression(Rate~of~population~change~(Delta)),
    y = expression(Bycatch~mortality~rate~(italic(F)[e]))
  ) +
  theme(
    legend.title = element_text(),
    legend.margin = margin(l=-8, unit="pt")
    ) 

fig1a <-
  fig1a +
    geom_point(
      data = bycatch_df %>% mutate(clade = stringr::str_to_title(clade)),
      aes(shape=clade), col="black"
      ) +
  scale_shape_manual(values = 21:24) +
  guides(
    fill = guide_colourbar(order = 1),
    shape = guide_legend(order = 2, title = NULL)
    ) + 
  coord_fixed()

#### Fig 1.B (Upsides FAO summary map) ####

overall_red <- 
  upsides %>%
  group_by(idoriglumped, regionfao, speciescat, speciescatname) %>%
  summarise(
    margc = mean(marginalcost, na.rm=T), 
    bet = mean(beta, na.rm=T),
    g = mean(g, na.rm=T),
    fvfmey = mean(eqfvfmey, na.rm=T),
    fvfmsy = mean(fvfmsy, na.rm=T),
    pctmey = mean(pctredfmey, na.rm=T),
    pctmsy = mean(pctredfmsy, na.rm=T)
    ) %>%
  mutate(
    fvfmey = ifelse(fvfmey==-Inf, NA, fvfmey),
    fvfmsy = ifelse(fvfmsy==-Inf, NA, fvfmsy),
    pctmey = ifelse(pctmey==-Inf, NA, pctmey),
    pctmsy = ifelse(pctmey==-Inf, NA, pctmsy)
    ) %>%
  mutate(
    wt = margc * ((g * fvfmsy)^bet),
    cstcurr = wt,
    cstmey = margc * (((g * fvfmsy)/fvfmey)^bet),
    cstmsy = margc * ((g)^bet)
    ) %>%
  ungroup() %>%
  mutate(wt = wt/sum(wt, na.rm = T)) %>%
  mutate(
    wtpctmey = wt * pctmey,
    wtpctmsy = wt * pctmsy
    ) 

fao_red <-
  overall_red %>% 
  select(-c(idoriglumped, speciescat, speciescatname)) %>%
  separate_rows(regionfao) %>%
  group_by(regionfao) %>%
  summarise_all(funs(sum(., na.rm=T))) %>% 
  group_by(regionfao) %>%
  mutate(
    avpctmey = 100 * (1 - (cstmey/cstcurr)),
    avpctmsy = 100 * (1 - (cstmsy/cstcurr))
    ) %>%
  mutate(
    fvfmey = ifelse(fvfmey==-Inf, NA, pctmey),
    fvfmsy = ifelse(fvfmsy==-Inf, NA, pctmey),
    pctmey = ifelse(pctmey==-Inf, NA, pctmey),
    pctmsy = ifelse(pctmey==-Inf, NA, pctmsy)
    ) 

## Load (and filter) FAO spatial data, before joining with the fao_red DF above
fao_sf <- 
  st_read("Data/Shapefiles/FAO_AREAS/FAO_AREAS.shp") %>%
  st_transform(proj_string) %>%
  filter(F_LEVEL=="MAJOR") %>%
  as_data_frame() %>%
  mutate(regionfao = as.character(F_AREA)) %>%
  left_join(fao_red)

## Plot the figure
fig1b <-
  ggplot() + 
  geom_sf(data = countries, fill = "white", col="white") +
  geom_sf(data = fao_sf, mapping = aes(fill = avpctmey/100), lwd = 0.08) +
  scale_fill_gradientn(
    name = "Reduction in fishing effort (MEY vs. 2010-2012)",
    colours = brewer_pal(palette = "YlOrRd")(9),
    labels = percent, limits=c(min(fao_sf$avpctmey)/100, 1)
    ) +
  guides(
    fill=guide_colourbar(barwidth=10.5, label.position="bottom", title.position="top")
    ) +
  theme(
    legend.title = element_text(), ## Turn legend text back on
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(l=-3, t=-5, unit="pt"),
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major = element_line(colour = "white")
  )

# fig1b

#### Composite Fig. 1 ####

fig1 <-
  ggdraw() +
  # draw_plot(figureName, xpos, ypos, width, height) +
  draw_plot(fig1a, 0, 0.49, 1, 0.49) +
  draw_plot(fig1b, 0, 0, 1, 0.49) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.49), size = 9)
# fig1

save_plot(
  "Figures/fig-1.png", fig1,
  base_width = fig_width,
  base_height = fig_width*1.5
  )
save_plot(
  "Figures/PDFs/fig-1.pdf", fig1,
  base_width = fig_width,
  base_height = fig_width*1.5,
  device=cairo_pdf
  )
rm(fig1, fig1a, fig1b)
dev.off()


##############################################
### Fig. 2 (NWA Loggerhead turtle example) ###
##############################################

sp_type <- "Loggerhead turtle (NW Atlantic)"

#### Fig 2.A (Range) ####

## Load shape file of NWA LH Regional Mgmt Units (based on Wallace et. al, PLoSONE 2010)
lh_rmus <- 
  read_sf("Data/Shapefiles/NW_Atl_Loggerhead/NW_Atl_Loggerhead_RMUs.shp") %>%
  st_transform(proj_string)

## Similarly for the the nesting sites
lh_nesters <- 
  read_sf("Data/Shapefiles/NW_Atl_Loggerhead/NW_Loggerhead_nesters.shp") %>%
  st_set_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(proj_string)

fig2a <-
  ggplot() +
  geom_sf(data = countries, col="black", fill="black", size = 0.1) +
  geom_sf(data = lh_rmus, col="#20b261", fill="#20b261", size = 0.1, alpha = 0.6) + 
  geom_sf(data = lh_nesters, col="#b22071", alpha = 0.2, size = 0.2) +
  theme(
    axis.line=element_blank(), 
    axis.ticks=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank(), 
    axis.title.x=element_blank(), axis.title.y=element_blank(),
    panel.grid.major = element_line(colour = "grey60", size = 0.5)
  )

# fig2a

#### Fig 2.B (Heatmap) ####

set.seed(123) ## First reset seed for disb on %T parameters

fig2b <- 
  sensrange_func(sp_type, 1000) %>% 
  mutate(pctT = pctT/100) %>%
  ggplot(aes(delta_draw, deltaN_draw, col = pctT)) +
  geom_point(alpha = 0.5) + 
  geom_point(shape = 21) +
  scale_color_distiller( 
    name = expression('%'~italic(T)), 
    palette = "YlOrRd", 
    trans = "reverse",
    labels = percent, limits = c(1, 0)
    ) +
  labs(
    x = expression(Delta),
    y = expression(Delta[n])
    ) +
  theme(
    legend.title = element_text(),
    legend.margin = margin(l=-8, unit="pt")
    ) 

#### Fig 2.C (Target species) ####
fig2c <- 
  stockselect_func(sp_type) %>%
  samples_plot() ## Note log scale

#### Fig 2.D (Bycatch reduction disb) ####
fig2d <- bycatchdist_plot(filter(results, species==sp_type), series = "MEY", truncate95 = T) 

#### Fig 2.E (Cost disb) ####
fig2e <- cost_plot(results %>% filter(species==sp_type), series = "MEY")

#### Fig 2.F (Targeting disb) ####
fig2f <- 
  targeting_plot(results %>% filter(species==sp_type), series = "MEY") +
  labs(x = "Targeting requirement")

#### Composite Fig. 2 ####

# ## Extract legend
# legend_fig2 <- g_legend(fig2d) 

### Tweak plots before putting theme together in composite figure
fig2d <- fig2d + theme(strip.text.x = element_blank(), legend.position = "none")
fig2e <- fig2e + theme(strip.text.x = element_blank(), legend.position = "none")
fig2f <- fig2f + theme(strip.text.x = element_blank(), legend.position = "none")

### Draw the figure (without legend)
fig2 <-
  ggdraw() +
  # draw_plot(fig, xpos,  ypos, width, height) +
  draw_plot(fig2a, 0,     0.67, 0.55, 0.33) +
  draw_plot(fig2b, 0.55, 0.67, 0.45, 0.33) +
  draw_plot(fig2c, 0,     0.34, 0.99, 0.33) +
  draw_plot(fig2d, 0,     0,    0.33, 0.33) +
  draw_plot(fig2e, 0.33,  0,    0.33, 0.33) +
  draw_plot(fig2f, 0.66,  0,    0.33, 0.33) +
  draw_plot_label(
    c("A", "B", "C", "D", "E", "F"), 
    c(0, 0.525, 0, 0, 0.33, 0.66), 
    c(1, 1, 0.66, 0.33, 0.33, 0.33), 
    size = 9
    )

# fig2

save_plot(
  "Figures/fig-2.png", fig2,
  base_width = 2 * fig_width,
  base_height = 2 * fig_width
  )
save_plot(
  "Figures/PDFs/fig-2.pdf", fig2,
  base_width = 2 * fig_width,
  base_height = 2 * fig_width, 
  device=cairo_pdf
  )

rm(fig2, fig2a, fig2b, fig2c, fig2d, fig2e, fig2f, lh_rmus, lh_nesters)
dev.off()


###############################
### Fig. 3 (Tradeoff plots) ###
###############################

fig_3mey <- tradeoffs_plot(results_summary, "MEY") + 
  theme(strip.text.y = element_text(angle = 90))  #https://github.com/tidyverse/ggplot2/issues/2356
fig_3mey + ggsave("Figures/fig-3-mey.png", width=2*fig_width, height=2*fig_width)
fig_3mey + ggsave("Figures/PDFs/fig-3-mey.pdf", width=2*fig_width, height=2*fig_width, device=cairo_pdf)
rm(fig_3mey)
dev.off()

fig_3msy <- tradeoffs_plot(results_summary, "MSY") + 
  theme(strip.text.y = element_text(angle = 90))  #https://github.com/tidyverse/ggplot2/issues/2356
fig_3msy + ggsave("Figures/fig-3-msy.png", width=2*fig_width, height=2*fig_width)
fig_3msy + ggsave("Figures/PDFs/fig-3-msy.pdf", width=2*fig_width, height=2*fig_width, device=cairo_pdf)
rm(fig_3msy)
dev.off()


###############################
### Fig. 4 (Recovery rates) ###
###############################

## First get a data frame of recovery rates by cost and targeting level
recovery_df <- recovery_func(results)

fig4 <- recovery_plot(recovery_df)
fig4 + ggsave("Figures/fig-4.png", width=2*fig_width, height=fig_width)
fig4 + ggsave("Figures/PDFs/fig-4.pdf", width=2*fig_width, height=fig_width, device=cairo_pdf)
rm(fig4)
dev.off()


#############################################
########### SUPPLEMENTARY FIGURES ########### 
#############################################

############################################################
##### Fig S1 (Upsides by FAO region & taxonomic group) #####
############################################################

## Read in taxonomy CSV for faceting categories
tax_df <- read_csv("Data/taxonomies.csv") 

## Start with `overall_red` DF created above (Fig. 1B)
fao_tax_red <-
  overall_red %>% 
  left_join(tax_df) %>%
  filter(taxonomy != "Not Included") %>%
  select(-c(idoriglumped, speciescatname)) %>%
  separate_rows(regionfao) %>%
  group_by(regionfao, taxonomy) %>% 
  summarise_all(funs(sum(., na.rm=T))) %>% 
  group_by(taxonomy, regionfao) %>%
  mutate(
    avpctmey = 100 * (1 - (cstmey/cstcurr)),
    avpctmsy = 100 * (1 - (cstmsy/cstcurr))
  ) %>%
  mutate(
    fvfmey = ifelse(fvfmey==-Inf, NA, pctmey),
    fvfmsy = ifelse(fvfmsy==-Inf, NA, pctmey),
    pctmey = ifelse(pctmey==-Inf, NA, pctmey),
    pctmsy = ifelse(pctmey==-Inf, NA, pctmsy)
  ) 

## Load (and filter) FAO spatial data, before joining with the `fao_tax_red` DF above
fao_tax_sf <- 
  st_read("Data/Shapefiles/FAO_AREAS/FAO_AREAS.shp") %>%
  st_transform(proj_string) %>%
  filter(F_LEVEL=="MAJOR") %>%
  as_data_frame() %>%
  mutate(regionfao = as.character(F_AREA)) 
fao_tax_sf <- 
  fao_tax_sf %>%
  ## Next step ensures all FAO regions are represented for each taxonomy (even if NA)
  right_join(
    crossing(
      taxonomy=unique(fao_tax_red$taxonomy), 
      regionfao=unique(fao_tax_sf$regionfao)
      )
    ) %>%
  left_join(fao_tax_red)
## Now join with the total reduction sf DF created in Fig 1B above
fao_tax_sf <-
  rbind(
    fao_tax_sf,
    fao_sf %>%
      mutate(
        taxonomy = "All",
        speciescat = 0
        ) %>%
      select_(.dots = colnames(fao_tax_sf))
    )

## Plot the figure
fig_s1 <-
  ggplot() + 
  geom_sf(data = countries, fill = "white", col="white") +
  geom_sf(data = fao_tax_sf, mapping = aes(fill = avpctmey/100), lwd = 0.08) +
  scale_fill_viridis(
    name = "Reduction in fishing effort (MEY vs. 2010-2012)",
    labels = percent, limits=c(min(fao_tax_sf$avpctmey)/100, 1)
    )  +
  guides(
    fill=guide_colourbar(barwidth=10.5, label.position="bottom", title.position="top")
    ) +
  facet_wrap(~taxonomy, ncol=2) +
  theme(
    legend.title = element_text(), ## Turn legend text back on
    legend.position = "bottom",
    legend.justification = "center", 
    legend.margin = margin(t=-5, unit="pt"),
    axis.line=element_blank(), axis.text.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major = element_line(colour = "white")
    )

fig_s1 + ggsave("Figures/fig-S1.png", width = 2*fig_width, height = 2*fig_width)
fig_s1 + ggsave("Figures/PDFs/fig-S1.pdf", width = 2*fig_width, height = 2*fig_width)
rm(fig_s1)
dev.off()


#############################################################
##### Fig S2 (Combined bycatch reduction distributions) #####
#############################################################
fig_s2 <- 
  bycatchdist_plot(results, combined_avg = T, truncate95 = T) +
  facet_wrap(~species, ncol = 3, scales = "free_x") +
  theme(legend.text = element_text(size = 7))
fig_s2 + ggsave("Figures/fig-S2.png", width=2.5*fig_width, height=2.5*fig_width*1.3)
fig_s2 + ggsave("Figures/PDFs/fig-S2.pdf", width=2.5*fig_width, height=2.5*fig_width*1.3, device=cairo_pdf)
rm(fig_s2)
dev.off()

################################################
##### Fig S3 (Combined cost distributions) #####
################################################
fig_s3 <- 
  cost_plot(results, combined_avg = T) +
  facet_wrap(~species, ncol = 3, scales = "free_x") +
  theme(legend.text = element_text(size = 7))
fig_s3 + ggsave("Figures/fig-S3.png", width=2.5*fig_width, height=2.5*fig_width*1.3)
fig_s3 + ggsave("Figures/PDFs/fig-S3.pdf", width=2.5*fig_width, height=2.5*fig_width*1.3, device=cairo_pdf)
rm(fig_s3)
dev.off()

############################################################
##### Fig S4 (Combined targeting change distributions) #####
############################################################
fig_s4 <- 
  targeting_plot(results, combined_avg = T) +
  facet_wrap(~species, ncol = 3, scales = "free_x") +
  theme(legend.text = element_text(size = 7))
fig_s4 + ggsave("Figures/fig-S4.png", width=2.5*fig_width, height=2.5*fig_width*1.3)
fig_s4 + ggsave("Figures/PDFs/fig-S4.pdf", width=2.5*fig_width, height=2.5*fig_width*1.3, device=cairo_pdf)

rm(fig_s4)
dev.off()


######################################
### Fig. S5 ("Pretty good" profit) ###
######################################

## If not already loaded from Fig. 4 above
# recovery_df <- recovery_func(results)

fig_s5 <- pgp_plot(recovery_df)
fig_s5 + ggsave("Figures/fig-S5.png", width=fig_width, height=fig_width)
fig_s5 + ggsave("Figures/PDFs/fig-S5.pdf", width=fig_width, height=fig_width, device=cairo_pdf)
rm(fig_s5)
dev.off()

#################################################
#### Fig. S6 (Theoretical alpha sensitivity) ####
#################################################
fig_s6 <-
  ggplot(data_frame(x = c(0, 5)), aes(x = x)) +
  stat_function(fun = function(x) (1 - (0.5^x)), size = 0.2) +
  geom_vline(aes(xintercept = 1), lty = 2) +
  geom_hline(aes(yintercept = 0.5), lty = 2) +
  scale_y_continuous(label = percent) +
  annotate(
    "text", 
    label = paste(expression(alpha==1)), 
    x = 1.5, y = 0.03, parse=T, family = font_type,
    size = 2
    ) +
  annotate(
    "text", 
    label = "        Reduction in target \nspecies mortality (50%)", 
    x = 3.5, y = 0.75, family = font_type,
    size = 2
    ) +
  labs(
    x = expression(alpha), 
    y = "Reduction in bycatch mortality"
  ) 
fig_s6 + ggsave("Figures/fig-S6.png", width = fig_width, height = fig_width)
fig_s6 + ggsave("Figures/PDFs/fig-S6.pdf", width = fig_width, height = fig_width)

rm(fig_s6)
dev.off()

############################################################
#### Figs. S7-S9 (Trade-off plots for sensitivity runs) ####
############################################################

s_runs <- c("fcorrected", "conservation", "alpha=05", "alpha=2", 
              "nonei", "weights", "2012only", "doubleuncert", "kitchen")

s_runs_list <- lapply(paste0("Results/bycatch_summary_results_", s_runs,".csv"), read_csv)

fig_s7_df <-
  bind_rows(
    s_runs_list[[1]] %>% mutate(sens = "A"),
    s_runs_list[[2]] %>% mutate(sens = "B"),
    s_runs_list[[3]] %>% mutate(sens = "C")
    )
fig_s8_df <-
  bind_rows(
    s_runs_list[[4]] %>% mutate(sens = "D"),
    s_runs_list[[5]] %>% mutate(sens = "E"),
    s_runs_list[[6]] %>% mutate(sens = "F")
    )
fig_s9_df <-
  bind_rows(
    s_runs_list[[7]] %>% mutate(sens = "G"),
    s_runs_list[[8]] %>% mutate(sens = "H"),
    s_runs_list[[9]] %>% mutate(sens = "I")
    )

## Find the low and high extremes of "delta_post" to set common bounds on the 
## x-axes across the figs.
x_lim <- 
  lapply(s_runs_list, function(df) {
    bounds <- 
      df %>% 
      filter(key=="pctredmey") %>% 
      mutate(delta_post = delta_mean + (fe_mean * (q50/100))) %>%
      summarise(low = min(delta_post), high = max(delta_post))
    }) %>%
  bind_rows() %>%
  summarise(low = min(low), high = max(high))

## Plot the figures
lapply(list(fig_s7_df, fig_s8_df, fig_s9_df), function(df) {
  p <-
    df %>%
    tradeoffs_plot("MEY") +
    xlim(c(x_lim$low, x_lim$high)) +
    facet_grid(clade~sens, scales = "free_y", space = "free", switch = "y") +
    theme(
      strip.text.y = element_text(angle = 90),  #https://github.com/tidyverse/ggplot2/issues/2356
      strip.text.x = element_text(face="bold", size = 9, hjust=0),
      panel.background = element_rect(fill = "#F2F2F2FF", colour = "#F2F2F2FF"),
      panel.spacing = unit(0.5, "lines")
    )
  fig_name <- ifelse(grepl(df$sens[1], "A"), "S7", ifelse(grepl(df$sens[1], "D"), "S8", "S9"))
  p + ggsave(paste0("Figures/fig-", fig_name, ".png"), width=3*fig_width, height=2*fig_width)
  p + ggsave(paste0("Figures/PDFs/fig-", fig_name, ".pdf"), width=3*fig_width, height=2*fig_width, device=cairo_pdf)
})

rm(s_runs, s_runs_list, x_lim, fig_s7_df, fig_s8_df, fig_s9_df)
dev.off()


###########################################################################
#### Fig. S10 (Illustrating role of uncertainty with the Māui dolphin) ####
###########################################################################

sp_type <- "Māui dolphin"
## Run 1 (main)
main_results <- read_csv("Results/bycatch_results.csv") %>% filter(species==sp_type)
## Run 9 (doubleuncert)
doubleuncert_results <- read_csv("Results/bycatch_results_doubleuncert.csv") %>% filter(species==sp_type)
## Run 10 (kitchen)
kitchen_results <- read_csv("Results/bycatch_results_kitchen.csv") %>% filter(species==sp_type)
 
#### Fig. 10 Bycatch distributions (panels A, D and G) ####
## Define common bounds on bycatch disb plots for comparison
lim_x <- quantile(doubleuncert_results$pctT, c(0.025, 0.975), na.rm = T)/100
## Fig 10.A (Bycatch distribution, Main run)
fig_s10a <- 
  bycatchdist_plot(main_results, series = "MEY", truncate95 = T) + 
  scale_x_continuous(
    name = "Reduction in mortality (MEY) vs. %T", 
    limits = lim_x, labels = percent
    ) + 
  theme(strip.text.x = element_blank(), legend.position = "none")
## Fig 10.D (Bycatch distribution, Double uncertainty run) 
fig_s10d <- 
  bycatchdist_plot(doubleuncert_results, series = "MEY", truncate95 = T) + 
  scale_x_continuous(
    name = "Reduction in mortality (MEY) vs. %T", 
    limits = lim_x, labels = percent
    ) +
  theme(strip.text.x = element_blank(), legend.position = "none")
## Fig 10.G (Bycatch distribution, Kitchen sink run) 
fig_s10g <- 
  bycatchdist_plot(kitchen_results, series = "MEY", truncate95 = T) + 
  scale_x_continuous(
    name = "Reduction in mortality (MEY) vs. %T", 
    limits = lim_x, labels = percent
    ) +
  theme(strip.text.x = element_blank(), legend.position = "none")

#### Fig. 10 Cost distributions (panels B, E and H) ####
## Fig 10.D (Costs distribution, Main run) 
fig_s10b <- 
  cost_plot(main_results, series = "MEY") + 
  theme(strip.text.x = element_blank(), legend.position = "none")
## Fig 10.E (Costs distribution, Double uncertainty run) 
fig_s10e <- 
  cost_plot(doubleuncert_results, series = "MEY") + 
  theme(strip.text.x = element_blank(), legend.position = "none")
## Fig 10.H (Costs distribution, Kitchen sink run) 
fig_s10h <- 
  cost_plot(kitchen_results, series = "MEY") + 
  theme(strip.text.x = element_blank(), legend.position = "none")

#### Fig. 10 Recovery plots (panels C, F and I) ####

## First get data frames of recovery rates by cost and targeting level
recovery_main <- recovery_func(read_csv("Results/bycatch_results.csv"))
recovery_doubleuncert <- recovery_func(read_csv("Results/bycatch_results_doubleuncert.csv"))
recovery_kitchen <- recovery_func(read_csv("Results/bycatch_results_kitchen.csv"))

## Fig 10.C (Recovery rate, Main run) 
fig_s10c <- recovery_plot(recovery_main, goal = "cost") + coord_cartesian()
## Fig 10.F (Recovery rate, Double uncertainty run) 
fig_s10f <- recovery_plot(recovery_doubleuncert, goal = "cost") + coord_cartesian()
## Fig 10.I (Recovery rate, Kitchen sink run) 
fig_s10i <- recovery_plot(recovery_kitchen, goal = "cost") + coord_cartesian()

#### Composite Fig. S 10 ####

fig_s10 <-
  ggdraw() +
  # draw_plot(fig, xpos,  ypos, width, height) +
  draw_plot(fig_s10a, 0,    0.66, 0.33, 0.33) +
  draw_plot(fig_s10b, 0.33, 0.66, 0.33, 0.33) +
  draw_plot(fig_s10c, 0.66, 0.66, 0.33, 0.33) +
  draw_plot(fig_s10d, 0,    0.33, 0.33, 0.33) +
  draw_plot(fig_s10e, 0.33, 0.33, 0.33, 0.33) +
  draw_plot(fig_s10f, 0.66, 0.33, 0.33, 0.33) +
  draw_plot(fig_s10g, 0,    0,    0.33, 0.33) +
  draw_plot(fig_s10h, 0.33, 0,    0.33, 0.33) +
  draw_plot(fig_s10i, 0.66, 0,    0.33, 0.33) +
  draw_plot_label(
    c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
    c(0, 0.33, 0.66, 0, 0.33, 0.66, 0, 0.33, 0.66), 
    c(1, 1, 1, 0.66, 0.66, 0.66, 0.33, 0.33, 0.33), 
    # size = 15
    size = 9
  )

# fig_s10

save_plot(
  "Figures/fig-S10.png", fig_s10,
  base_width = 3 * fig_width,
  base_height = 3 * fig_width
  )
save_plot(
  "Figures/PDFs/fig-S10.pdf", fig_s10,
  base_width = 3 * fig_width,
  base_height = 3 * fig_width, 
  device=cairo_pdf
  )

rm(fig_s10, fig_s10a, fig_s10b, fig_s10c, fig_s10d, fig_s10e, fig_s10f, fig_s10g, fig_s10h, fig_s10i)
dev.off()



# ##########################################################
# #### Replicas of Figs. S2 and S3 for sensitivity runs ####
# ##########################################################
# 
# ## Figs. S2 and S3 already made for main run (1) ##
# ## Remaining sensitivty runs 2-10 as described above
# sensitivity_runs <- 
#   c("fcorrected", "conservation", "alpha=05", "alpha=2", "nonei", 
#     "weights", "2012only", "doubleuncert", "kitchen")
# 
# ## Plot the figures over all sensitivity runs
# lapply(
#   sensitivity_runs, function(s){
#     s_df <- read_csv(paste0("Results/bycatch_results_", s, ".csv"), col_types=c("ddddc"))
#     
#     bycatchdist_plot(s_df, combined_avg = T, truncate95 = T) +
#       facet_wrap(~species, ncol = 3, scales = "free_x") + 
#       ggsave(paste0("Figures/SensitivityS2S3/fig-S2-", s, ".png"), width = 10, height = 13)
#     
#     cost_plot(s_df) +
#       facet_wrap(~species, ncol = 3, scales = "free_x") +
#     ggsave(paste0("Figures/SensitivityS2S3/fig-S3-", s, ".png"), width = 10, height = 13)
#     
#     Sys.sleep(4)
#   }
# )

### END ###