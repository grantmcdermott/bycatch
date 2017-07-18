## Clear environment
rm(list = ls())

#####################################
########### Load packages ###########
#####################################
library(data.table) ## Mostly for super fast reading/writing of large csv files
library(pbapply)
library(parallel)
library(R.utils)
library(scales)
library(grid)
library(png)
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
## Assign color scheme
bycatch_cols <- c("#ef3b2c","#386cb0","#fdb462","#7fc97f",
                  "#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
## Make some adjustments to the (now default) cowplot ggplot2 theme for figures
theme_update(
  text = element_text(family = font_type),
  legend.title = element_blank(),
  strip.background = element_rect(fill = "white"), ## Facet strip
  panel.spacing = unit(2, "lines") ## Increase gap between facet panels
)

## Decide on number of cores for parallel computation
num_cores <- detectCores() ## i.e. Use all available CPUs. Subtract 1 or 2 if you are running additional processes on your computer

## Load functions
source("R/bycatch_funcs.R")


#########################################################################
########## LOAD DATA AND CHOOSE ANALYSIS TYPE / SENSITIVITY RUN #########
#########################################################################

### Load bycatch data
bycatch_df <- read_csv("Data/bycatch_species.csv")
target_df <- read_csv("Data/target_species.csv")

## Get a vector of bycatch species
all_species <- bycatch_df$species 

## Set "alpha" parameter, i.e. elasticity of (changes in) bycatch to (changes 
## in) target stocks. Default is 1. Other options used in sensitivity analysis.
## See equation (S14). 
alpha_exp <- c(1, 0.5, 2)[1] 
alpha_str <- ifelse(alpha_exp==1, "", gsub("\\.","",paste0("_alpha=",alpha_exp))) 

## Select whether analysis simulates over 25% sensitivity range in Fe and delta
sensrange25 <- c(0, 1)[1] # 0 = off, 1 = on
sensrange_str <- ifelse(sensrange25==0, "", "_sensrange25")

## Select whether cost analysis simulates over 25% sensitivity range in bycatch weights
weights_sens <- c(0, 1)[1] # 0 = off, 1 = on
weights_str <- ifelse(weights_sens==0, "", "_weights")

## Select scenario (are all stocks going to MEY/MSY ["All stocks"],
##   or just those currently with F > Fmsy or B < Bmsy)
scenario <- c("All stocks", "Con. Concern")[1]
scenario_str <- ifelse(scenario=="All stocks", "", "_conservation")

### Load target stock data, derived from the "upsides" model of Costello et al. 
### (PNAS, 2016).
## First choose which version of the upsides data to use: 1) Main (includes NEI 
## stocks and uses a mean F over 2010-2012), 2) Excluding NEI stocks, 3) Use F
## from 2012 only. The main results of the paper rely on (1). The latter two are 
## used for sensitivity analysis in the SM.
upsides_type <- c("main", "nonei", "2012only")[1] ## Change as needed.
upsides_str <- ifelse(upsides_type=="main", "", paste0("_", upsides_type))
## Now read in the data
upsides <- 
  fread(paste0("Data/upsides_", upsides_type, ".csv")) %>% 
  as_data_frame()
# upsides <- read_csv(paste0("Data/upsides_", upsides_str, ".csv"), col_types = cols(regionfao = "c"))

## Decide whether to add in a "correction factor" for possible upward bias in C-MSY projections
## Default is 1 (i.e. no bias). Alternative is 2 (i.e. current F is twice as big as estimated by
## the upsides model.)
corr_factor <- c(1, 2)[1] ## Change as required 
corr_str <- ""
if(corr_factor==2){
  corr_str <- "_fcorrected"
  upsides <-
    upsides %>%
    mutate(curr_f = ifelse(dbase=="FAO", curr_f/corr_factor, curr_f)) %>%
    ## Adjust additionally affected variables in sequence
    mutate(
      fvfmsy = curr_f/g,
      eqfvfmey = curr_f/f_mey
    ) %>%
    mutate(
      pctredfmsy = 100 * (1 - (1/fvfmsy)),
      pctredfmey = 100 * (1 - (1/eqfvfmey)),
      pctredfmsycon = 100 * (1-(fconmsy/curr_f))
    )
}

## Convenience variable for output file name suffix
suff_str <- paste0(upsides_str, alpha_str, corr_str, scenario_str, sensrange_str, weights_str)

################################
########### ANALYSIS ###########
################################

#### WARNING: FULL ANALYSIS CAN TAKE A LONG TIME TO RUN (+/- 40 MIN ON A
#### 24 CORE LINUX SERVER). SKIP DIRECTLY TO FIGURES SECTION (LINE 156)    
#### TO PLOT PREVIOUSLY RUN (AND SAVED) RESULTS. 

### MCMC sampling parameters

## How many draws (states of the world) are we simulating for each species?
n1 <- 1000
## How many times do we sample (with replacement) over target stocks to resolve 
## uncertainty for a single draw?
n2 <- 100 

### Results for all species

## Apply the MCMC simulation function over all species (and bind into a common
## data frame). A series of progress bars will give you an indication of how
## long you have to wait, with one progress shown per bycatch species. In other 
## words, you'll see 20 progress bars in total if you run the full sample.
all_dt <- lapply(all_species, bycatch_func) %>% bind_rows() 

## Write results for convenient later use
write_csv(all_dt, paste0("Results/bycatch_results", suff_str, ".csv"))

## Get summary results and write those to disk too
results_summary <- summ_func(all_dt)
write_csv(results_summary, paste0("Results/bycatch_summary_results", suff_str, ".csv"))

####################################
########### MAIN FIGURES ########### 
####################################

## First, read the main results back in (no uncertainty, alpha = 1)
all_dt <- read_csv("Results/bycatch_results.csv")
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
  crossing(delta = seq(-40,0, length.out = 100), fe = seq(0,40, length.out = 100)) %>%
  mutate(z = abs(delta/fe)*100) %>%
  mutate(z = ifelse(z>100, 100, z)) %>%
  mutate_all(funs(./100)) %>%
  ggplot(aes(x = delta, y = fe)) + 
  geom_raster(aes(fill = z), interpolate = T) +
  scale_fill_gradientn(
    name = expression('%'~italic(T)),
    colours = brewer_pal(palette = "Spectral")(11), 
    trans = "reverse",
    labels = percent
  ) +
  geom_polygon(
    data=data_frame(delta=c(-40,-40,0)/100, fe=c(0,40,0)/100), 
    fill="#F2F2F2FF", col="#F2F2F2FF", lwd=1.5
    ) +
  labs(
    x = expression(Rate~of~population~change~(Delta)),
    y = expression(Bycatch~mortality~rate~(italic(F)[e]))
  ) +
  theme(legend.title = element_text())

fig1a <-
  fig1a +
  geom_point(
    data = bycatch_df %>% mutate(clade = stringr::str_to_title(clade)), 
    aes(shape=clade), fill="black", alpha=0.5, size = 3.5, stroke = 0
    ) +
    geom_point(
      data = bycatch_df %>% mutate(clade = stringr::str_to_title(clade)),
      aes(shape=clade), size = 3.5
    ) +
  scale_shape_manual(values = 21:24) +
  guides(
    fill = guide_colourbar(order = 1),
    shape = guide_legend(order = 2, title = NULL)
  ) + 
  coord_fixed()

## With animal silhouettes
# fig1a <-
#   fig1a +
#   lapply(sp_type, function(s){
#     delta <- (filter(bycatch_df, species %in% sp_type))$delta
#     fe <- (filter(bycatch_df, species %in% sp_type))$fe
#     j <- (bycatch_df %>% mutate(n = row_number()) %>% filter(species==s))$n
#     z <- (bycatch_df %>% filter(species==s))$silhouette
#     img <- readPNG(paste0("Figures/AnimalSilhouettes/",z,"-silhouette.png"))
#     g_img <- rasterGrob(img, interpolate=FALSE)
#     lapply(j, function(i) {
#       annotation_custom(g_img, xmin=delta[i]-0.025, xmax=delta[i]+0.025, ymin=fe[i]-0.025, ymax=fe[i]+0.025)
#     })
#   })
# fig1a +
#   xlim(-0.2, 0) +
#   ylim(0, 0.2)


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
  geom_sf(data = fao_sf, mapping = aes(fill = avpctmey/100), lwd = 0.25) +
  scale_fill_gradientn(
    name = "Reduction in fishing effort (MEY vs. 2010-2012)",
    colours = rev(brewer_pal(palette = "Spectral")(11)), #trans = "reverse",
    labels = percent, limits=c(min(fao_sf$avpctmey)/100, 1)
    ) +
  # scale_fill_viridis(
  #   name = "Reduction in fishing effort (MEY vs. 2010-2012)",
  #   labels = percent
  #   )  +
  guides(
    fill=guide_colourbar(barwidth=21, label.position="bottom", title.position="top")
  ) +
  theme(
    legend.title = element_text(), ## Turn legend text back on
    legend.position = "bottom",
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major = element_line(colour = "white")
  )

# fig1b

#### Composite Fig. 1 ####

fig1 <-
  ggdraw() +
  # draw_plot(figureName, xpos, ypos, width, height) +
  draw_plot(fig1a, 0, 0.5, 1, 0.5) +
  draw_plot(fig1b, 0, 0, 1, 0.5) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.475), size = 15)
# fig1

save_plot("Figures/fig-1.png", fig1,
          base_height = 10,
          base_aspect_ratio = 1/1.6#1/1.3
          )
save_plot("Figures/PDFs/fig-1.pdf", fig1,
          base_height = 10,
          base_aspect_ratio = 1/1.6#1/1.3
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
lh_nesters <- 
  read_sf("Data/Shapefiles/NW_Atl_Loggerhead/NW_Loggerhead_nesters.shp") %>%
  st_set_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(proj_string)

## Extract points from lh_nesters SF object and create normal data frame 
## Need to do this, since currently unable to change the size of points with geom_sf
## See: https://github.com/tidyverse/ggplot2/issues/2037
lh_nesters_df <- 
  as_data_frame(lh_nesters) %>%
  bind_cols(
    do.call(rbind, unclass(st_geometry(lh_nesters))) %>% 
      as_data_frame() %>%
      rename(x=V1, y=V2)
    )

## Bounding box
# # extent <- st_bbox(st_buffer(lh_rmus, 1000000))
# extent <- st_bbox(lh_rmus)
#       xmin       ymin       xmax       ymax 
# -9042320.1   463591.8  1603802.3  5358890.1 
## Manually adjust limits
extent <- c(-10000000.0, -500000.0, 2100000.0, 6250000.0)

fig2a <-
  ggplot() +
  geom_sf(data = countries, fill="bisque3", col="bisque3") +
  geom_sf(data = lh_rmus, aes(fill="Range"), col="dodgerblue", alpha = 0.7) + 
  scale_fill_manual(values = c("Range" = "dodgerblue")) +  
  scale_shape_manual(values = c("Nesting sites" = 21)) 

fig2a_inset <-
  fig2a +
  geom_point(data = lh_nesters_df, aes(x = x, y = y, shape = "Nesting sites"), fill="black", alpha = 0.2, size = 0.1) +
  geom_rect(
    data = data.frame(),
    aes(xmin = extent[[1]], xmax = extent[[3]], ymin = extent[[2]], ymax = extent[[4]]),
    colour = "red", fill = NA
    ) +
  ggthemes::theme_map() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "grey75")
  )

fig2a <-
  fig2a +
  geom_point(data = lh_nesters_df, aes(x = x, y = y, shape = "Nesting sites"), fill="black", alpha = 0.25, size = 2) +
  coord_sf(xlim = c(extent[[1]], extent[[3]]),ylim = c(extent[[2]], extent[[4]])) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.5, size = 3))) +
  theme(
    legend.text=element_text(size=14),
    legend.position = "bottom",
    axis.line=element_blank(), 
    axis.ticks=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank(), 
    axis.title.x=element_blank(), axis.title.y=element_blank(),
    panel.grid.major = element_line(colour = "grey75")
  )

# fig2a

#### Fig 2.B (Heatmap) ####

fig2b <-
  crossing(delta = seq(-10,0, length.out = 100), fe = seq(0,10, length.out = 100)) %>%
  mutate(z = abs(delta/fe)*100) %>%
  mutate(z = ifelse(z>100, 100, z)) %>%
  mutate_all(funs(./100)) %>%
  ggplot(aes(delta, fe, fill = z)) +
  geom_raster(interpolate = T) +
  scale_fill_gradientn(
    name = expression(Delta/~italic(F)[e]),#bquote(atop("Reduction in"~italic(F)[e], "to halt decline ")),
    colours = brewer_pal(palette = "Spectral")(11),
    trans = "reverse",
    labels = percent
  ) +
  geom_polygon(
    data=data_frame(delta=c(-10,-10,0)/100, fe=c(0,10,0)/100), 
    fill="#F2F2F2FF", col="#F2F2F2FF", lwd=0.75
    ) +
  labs(
    x = expression(Rate~of~population~change~(Delta)),
    y = expression(Bycatch~mortality~rate~(italic(F)[e]))
  ) +
  theme(legend.title = element_text())

fig2b <-
  fig2b +
  lapply(all_species[1], function(s){
    delta <- (filter(bycatch_df, species %in% sp_type))$delta
    fe <- (filter(bycatch_df, species %in% sp_type))$fe
    j <- (bycatch_df %>% mutate(j = row_number()) %>% filter(species==s))$j
    z <- (bycatch_df %>% filter(species==s))$silhouette
    img <- readPNG(paste0("Figures/AnimalSilhouettes/",z,"-silhouette.png"))
    g_img <- rasterGrob(img, interpolate=FALSE)
    lapply(j, function(i) {
      geom_img <- annotation_custom(g_img, xmin=delta[i]-0.004, xmax=delta[i]+0.004, ymin=fe[i]-0.004, ymax=fe[i]+0.004) 
      # geom_deltas <- geom_vline(xintercept = delta*c(.75, 1.25), lty = 2) ## Don't like vline extending above plot
      # geom_fes <- geom_hline(yintercept = fe*c(.75, 1.25), lty = 2) ## Don't like hline extending beyond plot
      geom_deltal <- geom_segment(y=-Inf, yend=.1, x=delta*.75, xend=delta*.75, lty=2)
      geom_deltau <- geom_segment(y=-Inf, yend=.1, x=delta*1.25, xend=delta*1.25, lty=2)
      geom_fel <- geom_segment(x=-Inf, xend=0, y=fe*.75, yend=fe*.75, lty=2) 
      geom_feh <- geom_segment(x=-Inf, xend=0, y=fe*1.25, yend=fe*1.25, lty=2)
      return(list(geom_img, geom_deltal, geom_deltau, geom_fel, geom_feh))
    }) 
  }) +
  scale_x_continuous(breaks=seq(-0.1, 0, by=0.02)) +
  scale_y_continuous(breaks=seq(0, 0.1, by=0.02)) 


#### Fig 2.C (Target species) ####
fig2c <- 
  stockselect_func(sp_type) %>%
  samples_plot() + ## Note log scale
  theme(strip.text = element_text(size = 14))

#### Fig 2.D (Bycatch reduction disb) ####
fig2d <- bycatchdist_plot(all_dt %>% filter(species==sp_type), "MEY") 

#### Fig 2.E (Cost disb) ####
fig2e <- cost_plot(all_dt %>% filter(species==sp_type), "MEY")

#### Fig 2.F (Targeting disb) ####
fig2f <- targeting_plot(all_dt %>% filter(species==sp_type), "MEY")

#### Composite Fig. 2 ####

## Extract legend
legend_fig2 <- g_legend(fig2d) 

### Tweak plots before putting theme together in composite figure
fig2d <- fig2d + theme(strip.text.x = element_blank(), legend.position = "none")
fig2e <- fig2e + theme(strip.text.x = element_blank(), legend.position = "none")
fig2f <- fig2f + theme(strip.text.x = element_blank(), legend.position = "none")

### Now, draw the figure (without legend)
fig2 <-
  ggdraw() +
  # draw_plot(fig, xpos,  ypos, width, height) +
  draw_plot(fig2a, 0.025, 0.67, 0.475, 0.33) +
  draw_plot(fig2b, 0.55,  0.67, 0.45,  0.33) +
  draw_plot(fig2c, 0,     0.34, 1,     0.33) +
  draw_plot(fig2d, 0,     0,    0.33,  0.33) +
  draw_plot(fig2e, 0.33,  0,    0.33,  0.33) +
  draw_plot(fig2f, 0.66,  0,    0.33,  0.33) +
  draw_plot_label(c("A", "B", "C", "D", "E", "F"), 
                  c(0, 0.525, 0, 0, 0.33, 0.66), 
                  c(1, 1, 0.66, 0.33, 0.33, 0.33), 
                  size = 15)

# fig2

save_plot("Figures/fig-2.png", fig2,
          base_height = 10,
          base_aspect_ratio = 1
          )
save_plot("Figures/PDFs/fig-2.pdf", fig2,
          base_height = 10,
          base_aspect_ratio = 1
          )

rm(fig2, fig2a, fig2a_inset, fig2b, fig2c, fig2d, fig2e)
dev.off()


###############################
### Fig. 3 (Tradeoff plots) ###
###############################

fig_3mey <- tradeoffs_plot(results_summary, "MEY")
fig_3mey + ggsave(paste0("Figures/fig-3-mey.png"), width=8, height=8)
fig_3mey + ggsave(paste0("Figures/PDFs/fig-3-mey.pdf"), width=8, height=8, device = cairo_pdf)
rm(fig_3mey)
dev.off()

fig_3msy <- tradeoffs_plot(results_summary, "MSY")
fig_3msy + ggsave(paste0("Figures/fig-3-msy.png"), width=8, height=8)
fig_3msy + ggsave(paste0("Figures/PDFs/fig-3-msy.pdf"), width=8, height=8, device = cairo_pdf)
rm(fig_3msy)
dev.off()


#############################################
########### SUPPLEMENTARY FIGURES ########### 
#############################################

#############################################################
##### Fig S1 (Upsides by FAO region & taxonomic group) #####
#############################################################

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
  geom_sf(data = fao_tax_sf, mapping = aes(fill = avpctmey/100), lwd = 0.25) +
  scale_fill_viridis(
    name = "Reduction in fishing effort (MEY vs. 2010-2012)",
    labels = percent
  )  +
  guides(
    fill=guide_colourbar(barwidth=21, label.position="bottom", title.position="top")
  ) +
  facet_wrap(~taxonomy, ncol=2) +
  theme(
    legend.title = element_text(), ## Turn legend text back on
    legend.position = "bottom",
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_line(colour = "white")
  )
fig_s1 + ggsave("Figures/fig-S1.png", width = 7, height = 7)
fig_s1 + ggsave("Figures/PDFs/fig-S1.pdf", width = 7, height = 7)
rm(fig_s1)
dev.off()

##############################################################
##### Fig S2 (Combined bycatch reduction distributions) #####
##############################################################
fig_s2 <- 
  bycatchdist_plot(all_dt) +
  facet_wrap(~species, ncol = 3, scales = "free_x") 
fig_s2 + ggsave("Figures/fig-S2.png", width = 10, height = 13)
fig_s2 + ggsave("Figures/PDFs/fig-S2.pdf", width = 10, height = 13, device = cairo_pdf)
rm(fig_s2)
dev.off()

#################################################
##### Fig S3 (Combined cost distributions) #####
#################################################
fig_s3 <- 
  cost_plot(all_dt) +
  facet_wrap(~species, ncol = 3, scales = "free_x")
fig_s3 + ggsave("Figures/fig-S3.png", width = 10, height = 13)
fig_s3 + ggsave("Figures/PDFs/fig-S3.pdf", width = 10, height = 13, device = cairo_pdf)
rm(fig_s3)
dev.off()

##############################################################
##### Fig S4 (Combined targeting change distributions) #####
##############################################################
fig_s4 <- 
  targeting_plot(all_dt) +
  facet_wrap(~species, ncol = 3, scales = "free_x")
fig_s4 + ggsave("Figures/fig-S4.png", width = 10, height = 13)
fig_s4 + ggsave("Figures/PDFs/fig-S4.pdf", width = 10, height = 13, device = cairo_pdf)
rm(fig_s4)
dev.off()

#################################################
#### Fig. S5 (theoretical alpha sensitivity) ####
#################################################
fig_s5 <-
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
fig_s5 + ggsave("Figures/fig-S5.png", width = 4, height = 4)
fig_s5 + ggsave("Figures/PDFs/fig-S5.pdf", width = 4, height = 4)
rm(fig_s5)
dev.off()

################################################
#### Figs. S6 and S7 (Sensitivity analysis) ####
################################################

## Sensitivity plots (possible analyses to include)

# 1. Main run: base year is 2010-2012 (geom. mean) 
#              alpha = 1
#              no uncertainty in Fe or delta
#              all stocks at MEY
df_run1 <- read_csv("Results/bycatch_summary_results.csv")

# 2. Estimated current F cut in half for all Catch-MSY stocks (to estimate effects of possible bias)
df_run2 <- read_csv("Results/bycatch_summary_results_fcorrected.csv")

# 3. Only conservation concern stocks (as defined by Costello et al. 2016)
#              rebuilt to MEY. 
df_run3 <- read_csv("Results/bycatch_summary_results_conservation.csv")

# 4. Main run with alpha = 0.5
df_run4 <- read_csv("Results/bycatch_summary_results_alpha=05.csv")

# 5. Main run with alpha = 2
df_run5 <- read_csv("Results/bycatch_summary_results_alpha=2.csv")

# 6. Main run with 'nei' stocks removed
df_run6 <- read_csv("Results/bycatch_summary_results_nonei.csv")

# 7. Main run with a 25% sensitivity range on Fe and delta
df_run7 <- read_csv("Results/bycatch_summary_results_sensrange25.csv")

# 8. Main run with a 25% sensitivity range on weights for target stock groups 
#              (e.g. demersals)
df_run8 <- read_csv("Results/bycatch_summary_results_weights.csv")

# 9. Main run with 2012 only as the base year (rather than 2010-2012)
df_run9 <- read_csv("Results/bycatch_summary_results_2012only.csv")


#### Fig. S6 (main run plus sensitivity analyses 2-5) ####

fig_s6 <-
  bind_rows(
    df_run1 %>% mutate(sens = "A"),
    df_run2 %>% mutate(sens = "B"),
    df_run3 %>% mutate(sens = "C"),
    df_run4 %>% mutate(sens = "D"),
    df_run5 %>% mutate(sens = "E")
    ) %>%
  tradeoffs_plot("MEY") +
  scale_x_continuous(expand = c(0.075, 0)) +
  facet_grid(clade~sens, scales = "free_y", space = "free", switch = "y") +
  theme(
    strip.text.x = element_text(face="bold", size = 16, hjust=0),
    panel.background = element_rect(fill = "#F2F2F2FF", colour = "#F2F2F2FF"),
    panel.spacing = unit(1, "lines")
    )

fig_s6 + ggsave("Figures/fig-S6.png", width=16, height=8)
fig_s6 + ggsave("Figures/PDFs/fig-S6.pdf", width=16, height=8, device=cairo_pdf)

rm(fig_s6); dev.off()

#### Fig. S7 (main run plus sensitivity analyses 6-9) ####

fig_s7 <-
  bind_rows(
    df_run1 %>% mutate(sens = "A"),
    df_run6 %>% mutate(sens = "B"),
    df_run7 %>% mutate(sens = "C"),
    df_run8 %>% mutate(sens = "D"),
    df_run9 %>% mutate(sens = "E")
    ) %>% 
  tradeoffs_plot("MEY") +
  scale_x_continuous(expand = c(0.075, 0)) +
  facet_grid(clade~sens, scales = "free_y", space = "free", switch = "y") +
  theme(
    strip.text.x = element_text(face="bold", size = 16, hjust=0),
    panel.background = element_rect(fill = "#F2F2F2FF", colour = "#F2F2F2FF"),
    panel.spacing = unit(1, "lines")
    )

fig_s7 + ggsave("Figures/fig-S7.png", width=16, height=8)
fig_s7 + ggsave("Figures/PDFs/fig-S7.pdf", width=16, height=8, device=cairo_pdf)

rm(fig_s7); dev.off()



##########################################################
#### Replicas of Figs. S2 and S3 for sensitivity runs ####
##########################################################

## Figs. S2 and S3 already made for main run (1) ##
## Remaining sensitivty runs 2-9 as described above
sensitivity_runs <- 
  c("fcorrected", "conservation", "alpha=05", "alpha=2", "nonei", "sensrange25", "weights", "2012only")

## Plot the figures over all sensitivity runs
lapply(
  sensitivity_runs, function(s){
    s_df <- read_csv(paste0("Results/bycatch_results_", s, ".csv"), col_types=c("ddddc"))
    
    bycatchdist_plot(s_df) +
      facet_wrap(~species, ncol = 3, scales = "free_x") + 
      ggsave(paste0("Figures/SensitivityS2S3/fig-S2-", s, ".png"), width = 10, height = 13)
    
    cost_plot(s_df) +
      facet_wrap(~species, ncol = 3, scales = "free_x") +
    ggsave(paste0("Figures/SensitivityS2S3/fig-S3-", s, ".png"), width = 10, height = 13)
    
    Sys.sleep(4)
  }
)

#### List of species categories ####
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