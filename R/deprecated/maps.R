library(ggplot2); library(RColorBrewer); library(sf); library(foreach); library(raster); library(dplyr)


# now reading from fishsuit, needs to be changed to new layers <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
lyr <- readRDS('../fishsuit/proc/pcrglobwb_anomaly_ensemble_perc_change.rds')[["3.2"]][['Tma']]

# # flow acc layer, in case needed
# fa <- raster('../fishsuit/data/pcrglobwb_hydrography/flowAcc_5min.tif')

# average discharge layer
q <- raster('../fishsuit/proc/ipsl/pcrglobwb_processed/merged/Qav_hist.tif')

# base layers
# crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
crs_custom <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# discharge cutoff for main map
cutoff <- 10 #m3/s

# rescale factor for main map
res_fac <- 2

# cropping of main map
crop_main <- c(-180, 180, -60, 90)

# bbox for 3 insets
insets_c <- data.frame(
  x = c(-57,20,88.5),
  y = c(-2.5,45,24.5)
)


insets_ts <- data.frame(
  x = c(-57.257,20.132,88.380),
  y = c(-2.622,45.210,24.376)
)


insets_side = 4

insets_ext <- foreach(i = 1:nrow(insets_c)) %do% extent(c(insets_c[i,'x']-insets_side, insets_c[i,'x']+insets_side,
                                                          insets_c[i,'y']-insets_side, insets_c[i,'y']+insets_side))
#--------------------------------------------------------------------------------------
# one plot for each variable at 4 different warming targets
cropping_poly <- st_bbox(extent(crop_main), crs = st_crs(4326)) %>% st_as_sfc(.)
inset_poly <- foreach(j = 1:length(insets_ext),.combine = 'c') %do% {st_bbox(insets_ext[[j]], crs = st_crs(4326)) %>% st_as_sfc(.) %>% st_transform(crs_custom)}

# base layers
world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)

# create mask layer based on q
qm <- q
qm[qm < cutoff] <- NA
qm[qm >= cutoff] <- 1

crs(qm) <- crs(lyr)

main <- as(lyr %>% mask(.,qm) %>% crop(crop_main) %>% aggregate(fact = res_fac, fun = max) %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
  as.data.frame(.) %>%
  dplyr::select(x,y,value = layer) %>%
  filter(!is.na(value)) %>%
  mutate(value = round(value*100,0))

main$value[main$value > 4] <- 4
main$value[main$value < -4] <- -4

# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = world, fill = 'black', lwd = NA) + # '#e4cead' yellow ochre
  geom_tile(data = main, aes(x=x, y=y, fill=value), alpha=0.8) +
  geom_sf(data = inset_poly, lwd = 1, color = 'black', fill = 'transparent') +
  scale_fill_gradient2(
    mid = 'white',
    low = 'blue',
    high = 'red',
    # mid = scales::muted('purple'),
    # low = scales::muted('blue'),
    # high = scales::muted('red'),
    midpoint = 0,
    breaks = seq(-4,4,1),
    labels = paste0(c('<3',seq(-3,0,1),paste0('+',seq(1,3,1)),'>3'),'%'),
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


inset <- foreach(j = 1:length(insets_ext),.combine = 'rbind') %do% {
  as(lyr %>% mask(.,qm) %>% crop(insets_ext[[j]]), "SpatialPixelsDataFrame") %>%
    as.data.frame(.) %>%
    mutate(inset_no = j) %>%
    dplyr::select(x,y,value = layer,inset_no) %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value*100,0))
}

inset$inset_no <- as.factor(inset$inset_no)
p_in <- ggplot() +
  geom_tile(data = inset, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_gradient2(    mid = 'white',
                           low = 'blue',
                           high = 'red',
                           midpoint = 0,
                           breaks = seq(-4,4,1),
                           labels = paste0(c('<3',seq(-3,0,1),paste0('+',seq(1,3,1)),'>3'),'%'),
                           na.value = 'transparent') +
  facet_wrap('inset_no',nrow = 1, scales = 'free') +
  coord_cartesian(expand = F) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.background = element_rect(fill = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        # axis.text = element_blank(),
        axis.ticks = element_line(color = 'black'),
        axis.title = element_blank(),
        legend.position = 'none',
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        aspect.ratio = 1
  )
# p_in

library(ggpubr)
ps <- ggarrange(p+theme(legend.margin = margin(-15,1,0,1),plot.margin = unit(c(0,-1,0,-1),'cm')),
                p_in+theme(plot.margin = unit(c(0,0.5,0,0.5),'cm')),
                nrow = 2,heights = c(1.8,1))
# ps

ggsave('figs/maps2.jpg',ps,
       width = 200,height = 175,dpi = 600,units = 'mm')

