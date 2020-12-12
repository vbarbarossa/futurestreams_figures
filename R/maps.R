source('config.R');

library(ggplot2); library(RColorBrewer); library(sf); library(foreach); library(raster); library(dplyr)


# now reading from fishsuit, needs to be changed to new layers <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
lyr <- readRDS('../fishsuit/proc/pcrglobwb_anomaly_ensemble_perc_change.rds')[["3.2"]][['Tma']]

# # flow acc layer, in case needed
# fa <- raster('../fishsuit/data/pcrglobwb_hydrography/flowAcc_5min.tif')

# average discharge layer
q <- raster('../fishsuit/proc/ipsl/pcrglobwb_processed/merged/Qav_hist.tif')

# base layers
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
crs_custom <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_transform(crs_custom)

# discharge cutoff for main map
cutoff <- 10 #m3/s

# rescale factor for main map
res_fac <- 1

#--------------------------------------------------------------------------------------
# one plot for each variable at 4 different warming targets

# create mask layer based on q
qm <- q
qm[qm < cutoff] <- NA
qm[qm >= cutoff] <- 1

crs(qm) <- crs(lyr)

main <- as(lyr %>% mask(.,qm) %>% aggregate(fact = res_fac, fun = max) %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
  as.data.frame(.) %>%
  mutate(var = 'Tma',varlong = 'max. weekly water temp.') %>%
  dplyr::select(x,y,value = layer,var,varlong) %>%
  filter(!is.na(value)) %>%
  mutate(value = round(value*100,0))

main$value[main$value > 5] <- 6
main$value[main$value < -5] <- -6

# dd$value[dd$value == Inf] <- max(dd$value[!is.infinite(dd$value)])
# dd$value[dd$value == -Inf] <- min(dd$value[!is.infinite(dd$value)])

# define different insets here
# inset <- as(lyr, "SpatialPixelsDataFrame") %>%
#   as.data.frame(.) %>%
#   mutate(warming = deg[i], var = 'Tma',varlong = 'max. weekly water temp.') %>%
#   dplyr::select(x,y,value = layer,warming,var,varlong) %>%
#   filter(!is.na(value)) %>%
#   mutate(value = round(value*100,0))

# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "grey80", lwd = 0.1) +
  geom_sf(data = world, fill = 'grey90', lwd = NA) + # '#e4cead' yellow ochre
  geom_tile(data = main, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_gradient2(mid = 'grey90',
                       low = scales::muted('blue'),
                       high = scales::muted('red'),
                       midpoint = 0,
                       breaks = seq(-6,6,1),
                       labels = paste0(c('<5',seq(-5,0,1),paste0('+',seq(1,5,1)),'>5'),'%'),
                       na.value = 'transparent') +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        panel.spacing = unit(0.1, "lines"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.title = element_blank()
  )
p

ggsave(paste0('figs/maps_pcrglobwb_anomaly_Qzf.jpg'),p,
       width = 200,height = 130,dpi = 600,units = 'mm')
