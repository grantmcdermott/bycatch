library(rgeos)
# library(tidyverse) ## Already loaded
library(rgdal)
library(broom)
library(viridis)

## Load FAO regions shapefile
fao_regs <- readOGR(dsn = "Data/FAO_AREAS", layer = "FAO_AREAS")
## Filter out only FAO major regions to map 
fao_major <- subset(fao_regs, fao_regs@data$F_LEVEL == "MAJOR")
## Make data frame of FAO major regions polygons and join with MEY/MSY reductions
## from the upsides model
fao_df <- 
  tidy(fao_major, region = "F_AREA") %>%
  rename(regionfao = id) %>%
  left_join(fao_red) %>%
  as_data_frame()

## Plot
fig1b <-
  fao_df %>%
  ggplot(aes(x = long, y = lat, group = group, fill = wtpctmey)) +
  geom_polygon() +
  scale_fill_viridis(
    name = "Reduction in fishing effort (MEY vs. 2012)",
    labels = percent,
    direction=-1)  +
  # scale_fill_gradientn(
  #   name = "Reduction in fishing effort (MEY vs. 2012)",
  #   # trans = "reverse",
  #   colours = rev(brewer_pal(palette = "Spectral")(11)), #limits = c(0,1), #brewer_pal(palette = "PRGn")(11)[1:5],
  #   labels = percent
  # ) +
  coord_map("gilbert") +
  coord_equal(ratio=1) + ## Uncommment for flat map projection
  guides(
    fill=guide_colourbar(barwidth=15, label.position="bottom", title.position="top")
    ) +
  theme(
    legend.title = element_text(), ## Turn legend text back on
    legend.position = "bottom",
    axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
    )
