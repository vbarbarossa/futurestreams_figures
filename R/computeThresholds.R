require(ncdf4)

start_time <- 1976
nmodels <- 20

serie_length <- 2100-start_time+1

readNC <- function(fileName){
  NC = nc_open(fileName)
  data = tail(ncvar_get(NC), serie_length)
  nc_close(NC)
  return(data)
}

runningMean <- function(data, years = 30){
  output = filter(data, rep(1/30,30))
  return(output)
}

files = list.files(path = 'data/thresholds',pattern = "*.nc", full.names = T)

count = 1
model = array()
scen = array()
series = matrix(NA, 20, serie_length)
for(file in files){
  series[count,] = readNC(file)
  model[count] = strsplit(file, "_")[[1]][4]
  scen[count] = strsplit(file, "_")[[1]][5]
  count = count + 1
}

plot(1,1, xlim=c(start_time,2100), ylim=c(-0.4,5), type="n")
for(m in 1:4){
  lines(start_time:2100, series[m,]-series[m,1], col=ceiling(m/4))
}

smoothSeries = apply(series, 1, runningMean)
plot(1,1, xlim=c(start_time,2100), ylim=c(-0.4,5), type="n")
for(m in 1:nmodels){
  lines(start_time:2100, smoothSeries[,m]-smoothSeries[15,m], col=ceiling(m/4))
}


exceedence = matrix(NA, nmodels, 4)

for(m in 1:nmodels){
  levCount = 1
  for(level in c(1.5, 2, 3.2,4.5)){
    sel = which(smoothSeries[15,m] + level < smoothSeries[,m])[1]
    exceedence[m,levCount] = c(start_time:2100)[sel]
    levCount = levCount + 1
  }
}

climate_models <- c('gfdl','hadgem','ipsl','miroc','noresm')
scenarios <- c('rcp2p6','rcp4p5','rcp6p0','rcp8p5')
warming_targets <- as.character(format( c(1.5, 2, 3.2,4.5), nsmall = 1))


rownames(exceedence) <- paste(rep(climate_models,each=4), rep(scenarios,5), sep = '_')
colnames(exceedence) <- warming_targets

# write.csv(as.data.frame(exceedence), "thresholdYears_4targets.csv", row.names = T)

# # alternatively, calculate delta first and then smooth 
# # does not make snese
# deltaSeries <- apply(series,1,function(x) {x - x[1]}) < wrong

# plot temperature anomaly
climate_models <- c('GFDL','HadGEM','IPSL','MIROC','NorESM')
scenarios <- c('2.6','4.5','6.0','8.5')

tab_ss <- as.data.frame(smoothSeries)
colnames(tab_ss) <- paste0(rep(climate_models,each = length(scenarios)),scenarios)



climI <- seq(1,nmodels,4)

library(foreach)
tab <- foreach(c = 1:length(climate_models),.combine = 'rbind') %do% {
  t <- tab_ss[,climI[c]:(climI[c]+3)]
  foreach(s = 1:length(scenarios),.combine = 'rbind') %do% {
    
    tt <- cbind(data.frame(year = start_time:2100),data.frame(Value = (t[,s] - t[15,s])))
    tt$GCM <- climate_models[c]
    tt$RCP <- scenarios[s]
    return(tt)
  }
  
}

library(ggplot2)
p <- ggplot(tab) +
  geom_line(aes(x = year, y = Value, color = GCM, linetype = RCP)) +
  geom_hline(yintercept = c(0,1.5,2,3.2,4.5),linetype = 4) +
  scale_y_continuous(breaks = c(0,1.5,2,3.2,4.5)) +
  ylab('Global Mean Air Temperature difference') +
  theme_bw() +
  theme(panel.grid = element_blank())

# ggsave(filename = 'figs/plot_airtemp_anomaly.jpg',p,width = 460,height = 400,dpi = 1000,unit = 'mm',scale = 0.35)

# do it only for the period 1976-2005
library(ggplot2)
p <- ggplot(tab) +
  geom_line(aes(x = year, y = Value, color = GCM, linetype = RCP)) +
  geom_hline(yintercept = 0:5,linetype = 2) +
  scale_y_continuous(breaks = 0:5) +
  scale_x_continuous(breaks = c(1976,2006,seq(2040,2100,30))) +
  ylab('Global Mean Air Temperature difference') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        # legend.position = 'left',
        legend.position = c(0.01,0.75),
        # legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.box.background = element_rect(color='white'),
        legend.justification = 'left')
p
ggsave(filename = 'figs/plot_airtemp_anomaly_1976.jpg',p,width = 130,height = 100,dpi = 300,unit = 'mm',scale = 1)
