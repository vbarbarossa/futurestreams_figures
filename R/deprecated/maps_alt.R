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
res_fac <- 2

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
#temporary resample to handle the raster

writeRaster(lyr,'proc/lyr_tmp.tif')

src.file <- 'proc/lyr_tmp.tif'
dst.file <- 'proc/lyr_tmp_res.tif'

r_orig <- raster(src.file)
target.res <- res(r_orig)*res_fac
crs.file <- crs(r_orig)@projargs

if(!file.exists(dst.file)) gdalUtils::gdalwarp(srcfile = src.file,dstfile = dst.file,
                                               s_srs = crs.file,t_srs = crs.file,
                                               tr = target.res,
                                               # srcnodata = '-1',dstnodata = '-1',
                                               r='max',overwrite = TRUE)
### end temporary resampling

r_res <- raster(dst.file)
r_res <- crop(r_res,extent(c(-180,180,-60,84)))

r_res <- lyr %>% aggregate(fact = res_fac, fun = max) %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb)
r_res <- crop(r_res,extent(c(-180,180,-60,84)))

length = 8.5
dpi = 600
labs = 0.5


wi = length
he = length*13/18
un = 'in'

insets_ext <-list(extent(c(-65,-61,-15,-11)),
                  extent(c(29,33,17,21)),
                  extent(c(103,107,9,13)))

r_upscaled <- r_res*100

points <- data.frame(matrix(ncol = 2,nrow = length(insets_ext)))
colnames(points) <- c('lon','lat')
#create shapefile from points
for(i in 1:length(insets_ext)) points[i,] <- c((insets_ext[[i]]@xmin + insets_ext[[i]]@xmax)/2,
                                               (insets_ext[[i]]@ymin + insets_ext[[i]]@ymax)/2)
p <- SpatialPoints(points)
crs(p) <- crs(r_orig)

library(tmap)
#MAIN MAP
main <- tm_shape(r_upscaled) +
  # tm_grid(projection="longlat", labels.size = labs, col = 'light grey',
  #         lwd = 0.01,alpha = 1,labels.inside.frame = FALSE) +
  
  tm_raster(palette="YlGnBu", contrast = c(0, 1), n = 9,
            auto.palette.mapping = FALSE, #saturation = 0.8,
            legend.show = TRUE, breaks = -3:6,title = expression('streamflow'~'['*m^{3}%.%s^{-1}*']')) + 
  
  tm_layout(legend.outside = TRUE, legend.outside.position = 'bottom',
            legend.position = c('left','bottom'),
            legend.format = list(fun = function(x) format(10**x,scientific = TRUE),
                                 text.separator = '-')
  ) +
  
  tm_shape(p) +  tm_dots(col = 'black',size = 0.3,shape = 0,border.lwd = 2)

#3 INSETS
insets <- list()
for(i in 1:length(insets_ext)){
  
  rc <- log10(crop(r_orig,insets_ext[[i]]))
  
  insets[[i]] <- tm_shape(rc) + 
    
    tm_grid(projection="longlat", n.x = 2,n.y = 2,
            labels.size = labs, col = 'light grey',
            lwd = 0.01,alpha = 1,labels.inside.frame = FALSE) +
    
    tm_raster(palette="YlGnBu", contrast = c(0, 1), n = 9,
              auto.palette.mapping = FALSE, #saturation = 0.8,
              legend.show = FALSE, breaks = -3:6) +
    
    tm_layout(outer.margins = c(0.1,0,0,0))
  
}
inset_size = 0.3
dist_insets = 0.25

library(grid)
tmap_save(main,insets_tm = insets,
          insets_vp=list(viewport(x= (0.6 - dist_insets), y= 0.25, width= inset_size, height= inset_size),
                         viewport(x=  0.6,                y= 0.25, width= inset_size, height= inset_size),
                         viewport(x= (0.6 + dist_insets), y= 0.25, width= inset_size, height= inset_size)), 
          filename="try1.jpg",width = wi, height = he, units = un,
          outer.margins=c(0.1,0.03,0,0.03), dpi = dpi)



### save maps for figure 1
scale_factor = 60
target.res <- res(r_orig)*scale_factor
dst.file <- paste0(wd,'maps/FLO1K.qav.long.term.19602015_resampledmaxval',scale_factor,'.tif')

if(!file.exists(dst.file)) gdalUtils::gdalwarp(srcfile = src.file,
                                               dstfile = dst.file,
                                               s_srs = crs.file,t_srs = crs.file,
                                               tr = target.res,
                                               srcnodata = '-1',dstnodata = '-1',
                                               r='max',overwrite = TRUE)

r_res <- raster(dst.file)
r_res <- crop(r_res,extent(c(-180,180,-60,84)))
r_upscaled <- log10(r_res)
r_upscaled[][is.infinite(r_upscaled[])] <- NA



main <- tm_shape(r_upscaled) +
  
  tm_raster(palette="YlGnBu", contrast = c(0, 1), n = 9,
            auto.palette.mapping = FALSE,breaks = -3:6,legend.show = FALSE) +
  tm_layout(bg.color = NA,frame = FALSE)

save_tmap(main, 
          filename="/vol/milkun3/Valerio/FLO1K/graphs/map.small.pdf",width = 6.5, height = 6.5/2.2, units = un,
          outer.margins=c(0.1,0.03,0,0.03))

save_tmap(main, 
          filename="/vol/milkun3/Valerio/FLO1K/graphs/map.small.jpg",width = 6.5, height = 6.5/2.2, units = un,
          outer.margins=c(0.1,0.03,0,0.03),dpi=300)
