require(ncdf4)
require(fields)

plotTemperature <- function(Temp, lon, lat, discharge, fileName = "averageWaterTemperature_annual.png", legendTitle = "Water temperature (Celcius)", plotHeight = 1500, plotWidth=2500, pointsize = 36, Qlim = 50, Tlim = 35, Tbreaks=0.25, PlotBreaks = 5, PlotFactor = 0.035, cols=colorRampPalette(rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))), sides=1, unit="C", ylim=c(-60,90), xlim=c(-168.5,191.5), size=1.3){
  mask = matrix(NA, length(lon), length(lat))
  mask[discharge > Qlim] = 1
  if(ylim[2] < 90){
    selLon = which(lon >= xlim[1] & lon <= xlim[2])
    selLat = which(lat >= ylim[1] & lat <= ylim[2])
    mask[1:selLon[1]] = NA
    mask[tail(selLon[2],1):length(lon)] = NA
    mask[1:selLat[1]] = NA
    mask[tail(selLat[2],1):length(lat)] = NA
  }
  
  demPlot = matrix(1, length(lon), length(lat))
  demPlot[is.na(Temp)] = NA
  
  toPlot = Temp
  toPlot[is.na(mask)] = NA
  toPlot[toPlot > Tlim] = Tlim
  if(sides == 2){
    toPlot[toPlot < -Tlim] = -Tlim
  }
  coordSel = which(is.na(toPlot) == FALSE, arr.ind=T)
  pointSel = which(is.na(toPlot) == FALSE)
  
  breaks = seq(0,Tlim,Tbreaks)
  Legendbreaks = seq(0,Tlim,PlotBreaks)
  if(sides == 2){
    breaks = seq(-Tlim,Tlim,Tbreaks)
    Legendbreaks = seq(-Tlim,Tlim,PlotBreaks)
  }
  cols = cols
  ncols = length(breaks)-1
  
  png(fileName, height=plotHeight, width = plotWidth, pointsize=pointsize)
  A = matrix(1, 9, 1)
  A[9,] = 2
  layout(A)
  par(mar=c(0.5,0.5,0.5,0.5), bg=NA)
  image(lon, lat ,demPlot, col = 1,
        main="", ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="")
  image(lon+360, lat ,demPlot, col = 1,
        main="", axes=FALSE, xlab="", ylab="", add=T)
  print(range(toPlot[pointSel], na.rm=T))
  col = cols(ncols)[as.numeric(cut(toPlot[pointSel], breaks=ncols))]
  points(lon[coordSel[,1]], lat[coordSel[,2]], col=col, 
         cex=(discharge[pointSel])^PlotFactor-1, pch=15)
  points(lon[coordSel[,1]]+360, lat[coordSel[,2]], col=col, 
         cex=(discharge[pointSel])^PlotFactor-1, pch=15)
  
  par(mar=c(0.2,0.2,0.2,0.2))
  plot(1, 1, xlim=c(0,1), ylim=c(0,1), type="n", axes=FALSE, xaxs="i", yaxs="i")
  symbols(seq(0.1, 0.9,length=ncols), rep(0.5,ncols), rectangles = matrix(rep(c(1/(ncols),0.4)), ncols, 2, byrow=TRUE), add=TRUE, fg=cols(ncols), bg=cols(ncols), inches=FALSE)
  text(seq(0.1,0.9, length=length(Legendbreaks)),rep(0.15, length(Legendbreaks)), paste(Legendbreaks, unit, sep=""), cex=size)
  text(0.5,0.85, legendTitle, cex=size)
  
  dev.off()
}

require(raster)

tf <- raster('../fishsuit/proc/ipsl/pcrglobwb_processed/merged/Tma_rcp8p5_4.5C_2074.tif')
th <- raster('../fishsuit/proc/ipsl/pcrglobwb_processed/merged/Tma_hist.tif')
td <- tf - th

q <- raster('../fishsuit/proc/ipsl/pcrglobwb_processed/merged/Qav_hist.tif')


str(as.matrix(td))
hist(td)

plotTemperature(Temp = t(as.matrix(td))[,1680:1],lon = seq(-180+1/24,180,1/12), lat = seq(-56+1/24,84,1/12),discharge = t(as.matrix(q))[,1680:1],
                ylim=c(-60,90), xlim=c(-180,180), sides=2, Tlim=10, legendTitle="Watertemperatue difference  (Celcius)", fileName = "figs/averageWaterTemperature_difference_rcp8p5_try.png")
