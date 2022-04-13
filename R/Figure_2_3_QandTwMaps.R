library(ggplot2); library(RColorBrewer); library(sf); library(foreach); library(raster); library(dplyr)
sf_use_s2(FALSE)

#<<<<<<<<<<<< try a different point shape, maybe one that is much smaller?

# library(ggplot2)
# df <- data.frame(x = 1:10, y = 1:10, s = 1:10)
# ggplot(df) +
#   geom_point(aes(x=x,y=y,size = s),shape = 20) +
#   scale_size(range = c(0,2)) +
#   theme_minimal()

# now reading from fishsuit, needs to be changed to new layers <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# lyr <- readRDS('../fishsuit/proc/pcrglobwb_anomaly_ensemble_perc_change.rds')[["3.2"]][['Tma']]

tf <- list.files(path = 'data/lyrs/',pattern = c('waterTemp.+rcp'), full.names = T) %>%
  lapply(raster) %>% brick() %>% calc(mean,na.rm=T)

th <- list.files(path = 'data/lyrs/',pattern = c('waterTemp.+hist'), full.names = T) %>%
  lapply(raster) %>% brick() %>% calc(mean,na.rm=T)

td <- tf - th

qf <- list.files(path = 'data/lyrs/',pattern = c('discharge.+rcp'), full.names = T) %>%
  lapply(raster) %>% brick() %>% calc(mean,na.rm=T)

qh <- list.files(path = 'data/lyrs/',pattern = c('discharge.+hist'), full.names = T) %>%
  lapply(raster) %>% brick() %>% calc(mean,na.rm=T)

qd <- qf - qh
qdl <- log10(qf) - log10(qh)

# average discharge layer
q <- qh

# base layers
crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# crs_custom <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# discharge cutoff for main map
cutoff <- 50 #m3/s

# cropping of main map
crop_main <- c(-180, 180, -60, 90)
# crop_main <- c(-90, -30, -60, 10)


# bbox for 3 insets
insets_c <- data.frame(
  x = c(-57,20,88.5),
  y = c(-2.5,45,24.5)
)


insets_ts <- data.frame(
  x = c(-57.2083,20.1250,88.3750),
  y = c(-2.6250,45.2083,24.3750),
  inset_no = 1:3
)

#1 -57.257 -2.622	## -57.2083, -2.625 Amazone
#2 20.132 45.210	## 20.125, 45.2083 Danube
#3 88.380 24.376	## 88.375 24.375 Ganges

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

# world <- read_sf('~/surfdrive/data/naturalearth/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp')[,1] %>%
#   st_crop(cropping_poly) %>%
#   st_transform(crs_custom)
# bb <- read_sf('~/surfdrive/data/naturalearth/ne_110m_graticules_all/ne_110m_wgs84_bounding_box.shp') %>%
#   st_crop(cropping_poly) %>%
#   st_transform(crs_custom)
# graticules <- read_sf('~/surfdrive/data/naturalearth/ne_110m_graticules_all/ne_110m_graticules_30.shp') %>%
#   st_crop(cropping_poly) %>%
#   st_transform(crs_custom)

# create mask layer based on q
qm <- q
qm[qm < cutoff] <- NA
qm[qm >= cutoff] <- 1
crs(qm) <- 4326

main <- rasterToPoints(q %>% mask(.,qm) %>% crop(crop_main) %>% mask(.,bb)) %>%
  as.data.frame() %>%
  rename(qav = layer) %>%
  st_as_sf(coords = c('x','y'))

st_crs(main) <- 4326

main$q <- raster::extract(qd,main)
main$ql <- raster::extract(qdl,main)
main$t <- raster::extract(td,main)

main <- main %>% st_transform(crs_custom)

main <- main %>%  
  bind_cols(data.frame(st_coordinates(.))) %>%
  as_tibble() %>%
  select(-geometry)


# main <- rasterToPoints(q %>% mask(.,qm) %>% crop(crop_main) %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb)) %>%
#   as.data.frame() %>%
#   rename(qav = Qav_hist)
# 
# main$q <- raster::extract(qd,main[,c('x','y')])
# main$t <- raster::extract(td,main[,c('x','y')])

# main$s <- log10(main$qav)
# main$s[main$s]

# WATER TEMPERATURE ###############################################################################

Tbreaks = 0.25
Tlim_up = 5
Tlim_lo = -5
# PlotBreaks = 1
PlotFactor = 0.035

br <- seq((Tlim_lo-1),(Tlim_up+1),1)

main$t[main$t > Tlim_up] <- Tlim_up+1
main$t[main$t < Tlim_lo] <- Tlim_lo-1

main$s <- main$qav %>% log10
main$s[main$s < 2.5] <- 1
main$s[main$s < 3 & main$s >= 2.5] <- 2
main$s[main$s < 4 & main$s >= 3] <- 3
main$s[main$s >= 4] <- 4

# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = world, fill = 'black', lwd = NA) + # '#e4cead' yellow ochre
  geom_point(data = main, aes(x=X, y=Y, size = qav, color=t), shape=15, stroke = 0) +
  geom_sf(data = inset_poly, lwd = 2, color = 'white', fill = 'transparent') +
  scale_color_gradientn(
    colors = inlmisc::GetColors(length(Tlim_lo:Tlim_up),scheme='sunset'),
    breaks = seq((Tlim_lo-1),(Tlim_up+1),1),
    labels = paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,1),paste0('+',seq(1,Tlim_up,1)),paste0('>',Tlim_up)),'°C'),
    na.value = 'transparent') +
  scale_size(range = c(0,2), trans = 'log10', guide = 'none') +
  theme_minimal(base_size = 50,
                base_line_size = 50,
                base_rect_size = 50) +
  theme(
    # base_size = 50,
    # text = element_text(size = 50),
    panel.grid.major = element_line(color=NA),
    panel.spacing = unit(0.01, "lines"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.key.width = unit(25,'line'),
    legend.key.height = unit(6,'line'),
    strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_text(angle = 0, vjust = -1, size = 12),
    legend.title = element_blank())
# p
# ggsave('figs/maps_try7.jpg',p,
#        width = 200, height = 125, dpi = 600,units = 'mm', scale = 6)

# create mask layer based on q
cutoff = 10
qm <- q
qm[qm < cutoff] <- NA
qm[qm >= cutoff] <- 1
crs(qm) <- 4326

inset <- foreach(j = 1:length(insets_ext),.combine = 'rbind') %do% {
  as(td %>% mask(.,qm) %>% crop(insets_ext[[j]]), "SpatialPixelsDataFrame") %>%
    as.data.frame(.) %>%
    mutate(inset_no = j) %>%
    dplyr::select(x,y,value = layer,inset_no) %>%
    filter(!is.na(value))
}

inset$value[inset$value >= max(br)] <- max(br)
inset$value[inset$value <= min(br)] <- min(br)


inset$inset_no <- as.factor(inset$inset_no)
p_in <- ggplot() +
  geom_tile(data = inset, aes(x=x, y=y, fill=value), alpha=0.8) +
  geom_point(data = insets_ts, aes(x=x, y=y ), color = 'white', shape = 3, size = 10, stroke=6) +
  scale_fill_gradientn(
    colors = inlmisc::GetColors(length(Tlim_lo:Tlim_up),scheme='sunset'),
    breaks = seq((Tlim_lo-1),(Tlim_up+1),1),
    labels = paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,1),paste0('+',seq(1,Tlim_up,1)),paste0('>',Tlim_up)),'°C'),
    na.value = 'transparent') +
  facet_wrap('inset_no',nrow = 1, scales = 'free') +
  coord_cartesian(expand = F) +
  ylab(' ') +
  theme_minimal(
    base_size = 50,
    base_rect_size = 50
  ) +
  theme(
    # text = element_text(size = 50),
    panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    # axis.text = element_blank(),
    axis.ticks = element_line(color = 'black'),
    axis.title.x = element_blank(),
    legend.position = 'none',
    strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    aspect.ratio = 1
  )
# p_in

# create time series tables for the three rivers
climate_models <- c('GFDL','HadGEM','IPSL','MIROC','NorESM')
rivers <- c('Amazon','Danube','Ganges')
rcps <- c('hist','2p6','4p5','6p0','8p5')

ts <- foreach(var = c('discharge','waterTemp'),.combine = 'rbind') %do%{
  foreach(riv = rivers,.combine = 'rbind') %do%{
    foreach(cl = climate_models,.combine = 'rbind') %do%{
      foreach(rc = rcps,.combine = 'rbind') %do%{
        
        v <- list.files(path = 'data/timeseries/',full.names = T, 
                        pattern = paste0(var,'.+',cl,'.+',rc,'.+',riv),ignore.case = T) %>%
          ncdf4::nc_open() %>% ncdf4::ncvar_get()
        
        if(rc != 'hist'){
          t <- data.frame(
            year = seq(2006,2100-1/52,1/52),
            value = v,
            GCM = cl,
            RCP = rc,
            river = riv,
            var = var
          )
          # correct year 2006 by pasting year 2007
          if(var == 'waterTemp') t$value[t$year >= 2006 & t$year < 2007] <- t$value[t$year >= 2007 & t$year < 2008]
          
        }else{
          t <- data.frame(
            year = seq(1976,2006-1/52,1/52),
            value = v,
            GCM = cl,
            RCP = rc,
            river = riv,
            var = var
          )
        }
        return(t)
      }
    }
  }
}

# # moving average function
# ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}

# calculate annual means
ts_annual <- foreach(var = c('discharge','waterTemp'),.combine = 'rbind') %do%{
  foreach(riv = rivers,.combine = 'rbind') %do%{
    foreach(cl = climate_models,.combine = 'rbind') %do%{
      foreach(rc = rcps,.combine = 'rbind') %do%{
        
        v <- list.files(path = 'data/timeseries/',full.names = T, 
                        pattern = paste0(var,'.+',cl,'.+',rc,'.+',riv),ignore.case = T) %>%
          ncdf4::nc_open() %>% ncdf4::ncvar_get()
        
        if(rc != 'hist'){
          
          v_hist <- list.files(path = 'data/timeseries/',full.names = T, 
                               pattern = paste0(var,'.+',cl,'.+','hist','.+',riv),ignore.case = T) %>%
            ncdf4::nc_open() %>% ncdf4::ncvar_get()
          
          t <- data.frame(
            year = rep(1997:2099,each = 52),
            value = c(v_hist[(52*(30-9)+1):length(v_hist)],v)
          )
          # correct year 2006 by pasting year 2007
          if(var == 'waterTemp') t$value[t$year == 2006] <- t$value[t$year == 2007]
          
          # compute annual means
          t <- t %>% group_by(year) %>% summarize(value = mean(value)) %>%
            mutate(GCM = cl,
                   RCP = rc,
                   river = riv,
                   var = var)
          
          
          
        }else{
          t <- data.frame(
            year = rep(1976:2005,each = 52),
            value = v
          )
          
          # compute annual means
          t <- t %>% group_by(year) %>% summarize(value = mean(value)) %>%
            mutate(GCM = cl,
                   RCP = rc,
                   river = riv,
                   var = var)
          
        }
        return(t)
      }
    }
  }
}

ts_annual_hist <- ts_annual %>% filter(RCP == 'hist')
ts_annual_fut <- ts_annual %>% filter(RCP != 'hist')
ts_annual <- rbind(
  ts_annual_fut %>% group_by(var,river,RCP,GCM) %>% summarise(value = zoo::rollmean(value,10),year = 2001:2094)
  ,ts_annual_hist %>% group_by(var,river,RCP,GCM) %>% summarise(value = zoo::rollmean(value,10),year = 1980:2000)
)

ts$RCP <- factor(as.factor(ts$RCP), levels = c("hist","2p6","4p5","6p0","8p5"))
p_ts <- ggplot(ts %>% filter(var == 'waterTemp') %>% 
                 mutate(value = (value-273.15))) +
  geom_line(aes(x = year,y = value, color = GCM, linetype = RCP), size = 1) +
  # geom_line(aes(x = year,y = zoo::rollmean(value,52*10,na.pad=T), color = GCM, linetype = RCP), size = 5) +
  geom_line(data = ts_annual%>% filter(var == 'waterTemp') %>% mutate(value = (value-273.15))
            ,aes(x = year,y = value, color = GCM, linetype = RCP)
            ,show.legend = F, size = 5) +
  # geom_smooth(aes(x = year,y = value, color = GCM, linetype = RCP), show.legend = F,
  #             method =  'gam',se = FALSE, size = 5) +
  # ylim(c(15,45)) +
  ylab('Water temperature [°C]') +
  facet_wrap('river', nrow=1) +
  theme_bw(
    base_size = 50
    # ,base_line_size = 50
    # ,base_rect_size = 50
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         linetype = guide_legend(override.aes = list(size = 5))) +
  
  theme(
    legend.position = 'bottom',
    legend.direction = 'vertical',
    legend.box = 'horizontal',
    legend.justification = 'left',
    panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.key.width = unit(5,"cm"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1,
    panel.spacing = unit(5,"lines")
  )
# p_ts

library(ggpubr)
ps <- ggarrange(
  p+theme(legend.margin = margin(-15,1,0,1),plot.margin = unit(c(0,-1,0,-1),'cm')),
  p_in+theme(plot.margin = unit(c(0,0.5,0,0.5),'cm')),
  p_ts + theme(axis.title.x = element_blank(),plot.margin = unit(c(0.5,1.2,0.5,1.2),'cm')),
  nrow = 3,
  heights = c(
    1.8,
    1,
    1.4
  )
)

ggsave('figs/maps_waterTemp2.jpg',ps,
       width = 200,height = 263,dpi = 300,units = 'mm', scale = 5, limitsize = F)

# DISCHARGE ###############################################################################

main$t <- main$ql # to reuse the script

Tbreaks = 0.1
Tlim_up = 0.5
Tlim_lo = -0.5
# PlotBreaks = 1
PlotFactor = 0.035

main$t[main$t > Tlim_up] <- Tlim_up+Tbreaks
main$t[main$t < Tlim_lo] <- Tlim_lo-Tbreaks

main$s <- main$qav %>% log10
main$s[main$s < 2.5] <- 1
main$s[main$s < 3 & main$s >= 2.5] <- 2
main$s[main$s < 4 & main$s >= 3] <- 3
main$s[main$s >= 4] <- 4

color_scheme <- inlmisc::GetColors(length(seq(Tlim_lo,Tlim_up,Tbreaks)),scheme='sunset',reverse = T)
br <- seq((Tlim_lo-Tbreaks),(Tlim_up+Tbreaks),Tbreaks)
lab <- paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,Tbreaks),paste0('+',seq(Tbreaks,Tlim_up,Tbreaks)),paste0('>',Tlim_up)),'')

# and draw
p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = world, fill = 'black', lwd = NA) + # '#e4cead' yellow ochre
  geom_point(data = main, aes(x=X, y=Y, size = qav, color=t), shape=15, stroke = 0) +
  geom_sf(data = inset_poly, lwd = 2, color = 'white', fill = 'transparent') +
  scale_color_gradientn(
    colors = color_scheme,
    breaks = br,
    labels = lab,
    na.value = 'transparent') +
  scale_size(range = c(0,2), trans = 'log10', guide = 'none') +
  theme_minimal(base_size = 50,
                base_line_size = 50,
                base_rect_size = 50) +
  theme(
    # base_size = 50,
    # text = element_text(size = 50),
    panel.grid.major = element_line(color=NA),
    panel.spacing = unit(0.01, "lines"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.key.width = unit(25,'line'),
    legend.key.height = unit(6,'line'),
    strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_text(angle = 0, vjust = -1, size = 12),
    legend.title = element_blank())
p
# ggsave('figs/maps_try7.jpg',p,
#        width = 200, height = 125, dpi = 600,units = 'mm', scale = 6)

# create mask layer based on q
cutoff = 10
qm <- q
qm[qm < cutoff] <- NA
qm[qm >= cutoff] <- 1
crs(qm) <- 4326

inset <- foreach(j = 1:length(insets_ext),.combine = 'rbind') %do% {
  as(qdl %>% mask(.,qm) %>% crop(insets_ext[[j]]), "SpatialPixelsDataFrame") %>%
    as.data.frame(.) %>%
    mutate(inset_no = j) %>%
    dplyr::select(x,y,value = layer,inset_no) %>%
    filter(!is.na(value))
}

inset$value[inset$value >= max(br)] <- max(br)
inset$value[inset$value <= min(br)] <- min(br)

inset$inset_no <- as.factor(inset$inset_no)
p_in <- ggplot() +
  geom_tile(data = inset, aes(x=x, y=y, fill=value), alpha=0.8) +
  geom_point(data = insets_ts, aes(x=x, y=y ), color = 'white', shape = 3, size = 10, stroke=6) +
  scale_fill_gradientn(
    colors = color_scheme,
    breaks = br,
    labels = lab,
    na.value = 'transparent') +
  facet_wrap('inset_no',nrow = 1, scales = 'free') +
  coord_cartesian(expand = F) +
  ylab(' ') +
  theme_minimal(
    base_size = 50,
    base_rect_size = 50
  ) +
  theme(
    # text = element_text(size = 50),
    panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    # axis.text = element_blank(),
    axis.ticks = element_line(color = 'black'),
    axis.title.x = element_blank(),
    legend.position = 'none',
    strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    aspect.ratio = 1
  )
# p_in

p_ts <- ggplot(ts %>% filter(var == 'discharge')) +
  geom_line(aes(x = year,y = value, color = GCM, linetype = RCP), size = 1) +
  geom_line(data = ts_annual%>% filter(var == 'discharge')
            ,aes(x = year,y = value, color = GCM, linetype = RCP)
            ,show.legend = F, size = 5) +
  ylab('Discharge [m3/s]') +
  facet_wrap('river', nrow=1,scales = 'free_y') +
  theme_bw(
    base_size = 50
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         linetype = guide_legend(override.aes = list(size = 5))) +
  theme(
    legend.position = 'bottom',
    legend.direction = 'vertical',
    legend.box = 'horizontal',
    legend.justification = 'left',
    panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.key.width = unit(5,"cm"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1,
    panel.spacing = unit(5,"lines")
  )
p_ts

library(ggpubr)
ps <- ggarrange(
  p+theme(legend.margin = margin(-15,1,0,1),plot.margin = unit(c(0,-1,0,-1),'cm')),
  p_in+theme(plot.margin = unit(c(0,0.5,0,0.5),'cm')),
  p_ts + theme(axis.title.x = element_blank(),plot.margin = unit(c(0.5,1.2,0.5,1.2),'cm')),
  nrow = 3,
  heights = c(
    1.8,
    1,
    1.4
  )
)

ggsave('figs/maps_discharge2.jpg',ps,
       width = 200,height = 263,dpi = 300,units = 'mm', scale = 5, limitsize = F)
