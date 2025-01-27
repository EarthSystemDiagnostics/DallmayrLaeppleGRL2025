## Aim: Read B31 and B32; plot the simulated isotope profile 
## Status/belongs to:  Topography hypothesis manuscript
## Result: 
## Next step: 

basedrive="/Users/tlaepple/data/" 
library(dplyr)
library(ggplot2)
library(sf)
library(stars)
library(deming)
library(zoo)
library(RColorBrewer)
setwd(paste(basedrive,"topohypothesis",sep=""))



## Read the slope, altitude and velocity fields
sl_proxy <-  stars::read_stars("./data/topography/slope_250km_5kmTransect_20pts.tif", proxy = T)
altitude_proxy <- stars::read_stars("./data/topography/REMA_200m_dem_filled.tif",proxy=T)
vel_proxy <- stars::read_stars("./data/topography/antarctic_ice_vel_phase_map_v01.tif", proxy = T)    


### Functions to extract 'records' from the 2D fields


##' @title Create a 'slope' and 'altitude' ice core record according to the provided time vector
##' @param pnt Point of the 'ice core'
##' @param sl_proxy stars field of the slope
##' @param vel_proxy stars field of the velocity
##' @param altitude_proxy stars field of the altitude 
##' @param outTime time vector on which to interpolate the results
##' @param maxSpeed maximum speed above NA's are returned
##' @return ts object containing a timeseries object wth two timeseries (time unit years) of the slope and altitude from which the ice originated; or NA's if the speed is too fast
##' @author Thomas Laepple 
createRecord <- function(pnt,sl_proxy,vel_proxy,altitude_proxy,outTime=seq(from=0,to=10000,by=50),maxSpeed=50,adjust=1)
{
  
  ##Extract the velocity at the starting point
  start <- st_point(pnt) %>% st_sfc() %>% st_set_crs(4326) %>% st_transform(st_crs(vel_proxy))
  vel <- stars::st_extract(vel_proxy, start)  %>% unlist()
  speed = sqrt(vel[1]^2 + vel[2]^2)
  
  ##if the speed is too fast, return a time-series of NA
  if (speed > maxSpeed) return(ts(rep(NA,length(outTime)),start=min(outTime),deltat=diff(outTime)[1]))
  
  ##Calculate the length needed if the starting speed + 20%
  traj_length_m<-max(outTime) * speed * 1.2
  
  #Caclulate the sampling distance.. minimum 200m as there is no more information on the DEM
  sampling_distance_m<-max(200,diff(outTime)[1]*speed)
  
  
  result<-topoExtract.2param(pnt,vel_proxy,sl_proxy,altitude_proxy,traj_length_m, sampling_distance_m, interpolate = FALSE,namefield1 = 'slope',namefield2 = 'altitude')
  timevector <- c(0,(cumsum(1/(result$speed*adjust))*diff(result$dist)[1])[-length(result$speed)])
  
  if (sum(!is.na(result$slope))<1) return(NA)
  
  return(ts(cbind(approx(timevector,result$slope,outTime)$y,approx(timevector,result$altitude,outTime)$y),start=min(outTime),deltat=diff(outTime)[1]))
  
}




topoExtract.2param <- function(p,   vel_proxy, field1_proxy, field2_proxy,traj_length_m, sampling_distance_m, interpolate = FALSE,namefield1 = 'field1',namefield2='field2') {
  
  invisible(lapply(c("sf", "stars", "dplyr", "geosphere"), require, 
                   character.only = TRUE))
  invisible(sf_use_s2(FALSE))
  
  start <- st_point(p) %>% st_sfc() %>% st_set_crs(4326) 
  
  ## start point
  res <-   c(start %>% st_coordinates(),
             (st_extract(field1_proxy, start %>% st_transform(st_crs(field1_proxy)), bilinear = interpolate) %>% unlist())[1],
             (st_extract(field2_proxy, start %>% st_transform(st_crs(field2_proxy)), bilinear = interpolate) %>% unlist())[1],
             (st_extract(vel_proxy, start %>% st_transform(st_crs(vel_proxy)), bilinear = interpolate) %>% unlist())) %>% 
    matrix(ncol = 6, nrow = 1) %>% as_tibble() %>% setNames(c("lon", "lat", namefield1, namefield2, "vx", "vy")) %>%
    mutate(dist = 0, dir = ((atan2(vx, vy) * (180/pi) + 360) +180) %% 360, speed = sqrt(vx^2 + vy^2))
  
  ## forward loop
  repeat{
    nextP <- destPoint(res[nrow(res), c("lon", "lat")] %>% as.matrix(), res$dir[nrow(res)], sampling_distance_m) %>% 
      st_point() %>% st_sfc() %>% st_set_crs(4326) 
    
    res <- res %>% rbind(
      c(nextP %>% st_coordinates(),
        (st_extract(field1_proxy, nextP %>% st_transform(st_crs(field1_proxy)), bilinear = interpolate) %>% unlist())[1],
        (st_extract(field2_proxy, nextP %>% st_transform(st_crs(field2_proxy)), bilinear = interpolate) %>% unlist())[1],
        (st_extract(vel_proxy, nextP %>% st_transform(st_crs(vel_proxy)), bilinear = interpolate) %>% unlist())) %>% 
        matrix(ncol = 6, nrow = 1) %>% as_tibble() %>% setNames(c("lon", "lat", namefield1, namefield2, "vx", "vy")) %>%
        mutate(dist = res$dist[nrow(res)] + sampling_distance_m, dir = ((atan2(vx,vy) * (180/pi) + 360) + 180) %% 360, speed = sqrt(vx^2 + vy^2))
    )
    
    if(res$dist[nrow(res)]>=traj_length_m) break
  }
  res
}


#Extract the topography (slope and altitude) from the upstream trajectory including the uncertainty in flow speed

Extract <- function(coord,speeduncertainty=c(1/1.3,1.3))
{
  
  rt<-createRecord(matrix(coord,ncol=2,nrow=1),sl_proxy,vel_proxy,altitude_proxy)
  rt.s<-createRecord(matrix(coord,ncol=2,nrow=1),sl_proxy,vel_proxy,altitude_proxy,adjust=speeduncertainty[1])
  rt.f<-createRecord(matrix(coord,ncol=2,nrow=1),sl_proxy,vel_proxy,altitude_proxy,adjust=speeduncertainty[2])
  return(list(mean=rt,slow=rt.s,fast=rt.f))
}

sr <- function(tss) return(tss - mean(window(tss,0,500),na.rm=TRUE)) ##Remove the mean of the first 500 years

#Translate to climate anomalies using the range of slope / climate and altitude / climate relationships
Translate2Climate <- function(rt,slope2climate=c(0.39,0.5,0.7),altitude2climate=c(-0.0025)) 
{
  x<-c(time(rt$mean[,1]))
  y1<-sr(rt$slow[,1]*1000*slope2climate[3])
  y2<-sr(rt$slow[,1]*1000*slope2climate[1])
  y3<-sr(rt$fast[,1]*1000*slope2climate[3])
  y4<-sr(rt$fast[,1]*1000*slope2climate[1])
  y<-sr(rt$mean[,1]*1000*slope2climate[2])
  
  ymax <- y
  ymin <- y
  ymax[]<-apply(cbind(y1,y2,y3,y4),1,max)
  ymin[]<-apply(cbind(y1,y2,y3,y4),1,min)
  
  yaltitude <- sr(rt$mean[,2]*altitude2climate)
  return(list(x=x,ymin=ymin,ymax=ymax,y=y,yaltitude=yaltitude))
}
###

#Read B31 and B32 ice-core data

b31.raw<-read.table("./data/icecore/B31_Isotope_Age.txt",sep="\t",header=TRUE)

b31<-ts(rev(b31.raw$d18O.H2O..per.mil.SMOW.),start=last(b31.raw$Age..year.AD.))
b31.accum<-ts(rev(b31.raw$Acc.rate.ice..kg.m..2.a.), start=last(b31.raw$Age..year.AD.))


b32.raw<-read.table("./data/icecore/B32_Isotope_Age.txt",sep="\t",header=TRUE)
b32<-ts(b32.raw$d18O.H2O..per.mil.SMOW.,start=last(b32.raw$Age..year.AD.))
b32.accum<-ts(b32.raw$Acc.rate.ice..kg.m..2.a., start=last(b32.raw$Age..year.AD.))



coord.B31<-matrix(c(-3.4303,-75.5815),ncol=2,nrow=1)
coord.B32<-matrix(c(-0.0073,-75.0025),ncol=2,nrow=1)



#Some defs for a common plotting style

colors <- brewer.pal(4,"Set1")
colors[1]<-"black"
SC = 1.75 #Scaling of the plots to get the right font size

rgb_values <- col2rgb(colors[3]) / 255  
col.transparent <- rgb(rgb_values[1], rgb_values[2], rgb_values[3], alpha = 0.3)

ylab.accum = expression("Accum. anomaly [kg" ~ m^-2 ~ a^-1 ~"]")
xlab.time <- "Time [yr BP]"
ylab.d18O <- expression(delta^{18}*O~"anomaly ["*"\u2030"*"]")


pdf(file = "./plots/Figure3_left.pdf", width = 3.75*SC, height = 4.2*SC)

result<-topoExtract.2param(coord.B31,vel_proxy,sl_proxy,altitude_proxy,20000, 200, namefield1 = 'slope',namefield2 = 'altitude')

timevector <- c(0,(cumsum(1/(result$speed))*diff(result$dist)[1])[-length(result$speed)])
yearticks <- seq(0,8000,by=500)
at=approx(timevector,result$dist,yearticks)$y/1000
xlimval <- c(0,approx(timevector,result$dist,3000)$y/1000) #Show only last 3000 years



par(mfcol=c(2,1))
#increase the margin on the right side
par(mai=c(0.8,0.8,0.8,0.8))

plot(result$dist/1000, result$altitude, type="l", col=colors[2], ylab="Elevation [m]", xlab="Distance [km]",ylim=c(2650,2693),xlim=xlimval,axes=FALSE)
box()
axis(1)
axis(2,at=c(2670,2680,2690),col=colors[2],col.axis=colors[2])

# Set up the second y-axis (slope) on the right
par(new=TRUE)

plot(result$dist/1000, result$slope*1000, type="l", col=colors[3], axes=FALSE, xlab="", ylab="",ylim=c(-2.5,8),xlim=xlimval)
axis(4,at=c(-2,0,2,4),col=colors[3],col.axis=colors[3])  # Add an axis on the right side

# Label for the second y-axis with a red text color
mtext("Slope [m/km]", side=4, line=3, col="black") 
abline(h=0,lty=2)

axis(3,at=at,labels=round(yearticks,0))
mtext(xlab.time, side=3, line=3)  # Label for the second y-axis


#####


result<-topoExtract.2param(coord.B32,vel_proxy,sl_proxy,altitude_proxy,20000, 200, namefield1 = 'slope',namefield2 = 'altitude')


timevector <- c(0,(cumsum(1/(result$speed))*diff(result$dist)[1])[-length(result$speed)])
yearticks <- seq(0,8000,by=500)
at=approx(timevector,result$dist,yearticks)$y/1000

xlimval <- c(0,approx(timevector,result$dist,3000)$y/1000) #Show only last 3000 years

#increase the margin on the right side
par(mai=c(0.8,0.8,0.8,0.8))

plot(result$dist/1000, result$altitude, type="l", col=colors[2], ylab="Elevation [m]", xlab="Distance [km]",ylim=c(2850,2893),xlim=xlimval,axes=FALSE)
box()
axis(1)
axis(2,at=c(2670,2680,2690)+200,col=colors[2],col.axis=colors[2])

# Set up the second y-axis (slope) on the right
par(new=TRUE)

plot(result$dist/1000, result$slope*1000, type="l", col=colors[3], axes=FALSE, xlab="", ylab="",ylim=c(-2.5,8),xlim=xlimval)
axis(4,at=c(-2,0,2,4),col=colors[3],col.axis=colors[3])  # Add an axis on the right side

mtext("Slope [m/km]", side=4, line=3, col="black") 
abline(h=0,lty=2)


axis(3,at=at,labels=round(yearticks,0))
mtext(xlab.time, side=3, line=3)  

####
dev.off()


b31.upstream <- Extract(coord.B31,speeduncertainty = c(1/1.3,1.3)) #Assume 30% speed uncertainty
b32.upstream <- Extract(coord.B32,speeduncertainty = c(0.2,1))    #Assume that it could be much (5x slower) and max is the value from the satellite

b31.correct<-Translate2Climate(b31.upstream,slope2climate=c(0.24,0.42,0.6),altitude2climate = -3.39/1000)
b32.correct<-Translate2Climate(b32.upstream,slope2climate=c(0.24,0.42,0.6),altitude2climate = -3.39/1000)



pdf(file = "./plots/Figure3_right.pdf", width = 3.75*SC, height = 4.2*SC)

par(mfcol=c(2,1))
par(mai=c(0.8,0.8,0.8,0.8))
### Plot B31

#Switch to years BP
isotope.obs <- b31
isotope.obs.pTs <- ts(rev(c(isotope.obs)), start=rev(1950 - c(time(isotope.obs)))[1])

# plot of isotope timeseries
plot(sr(isotope.obs.pTs), type="l", xlim=c(0, 3000), ylim=c(-2, 4), 
     col=colors[[1]], lwd=1, xlab=xlab.time, ylab=ylab.d18O, main="B31")

# Add the shaded area for uncertainty range (ymin to ymax); use the same normalization as the main predicted record
mv <- mean(window(b31.correct$y,0,500),na.rm=TRUE)


polygon(c((b31.correct$x), rev((b31.correct$x))),
        c((b31.correct$ymax-mv), rev((b31.correct$ymin-mv))),
        col=col.transparent, border=NA)  # Semi-transparent red shading

# Add lines for corrected data
lines(sr(b31.correct$yaltitude), col=colors[2], lwd=2)   
lines(sr(b31.correct$y), col=colors[3], lwd=2)    
  


  ### Now for B32
  
  ###
#Switch to years BP
  isotope.obs <- b32
  isotope.obs.pTs <- ts(rev(c(isotope.obs)), start=rev(1950 - c(time(isotope.obs)))[1])
  
  # plot of isotope timeseries
  plot(sr(isotope.obs.pTs), type="l", xlim=c(0, 3000), ylim=c(-2, 4), 
       col=colors[[1]], lwd=1, xlab=xlab.time, ylab=ylab.d18O, main="B32")
  
  # Add the shaded area for uncertainty range (ymin to ymax); use the same normalization as the main predicted record
  mv <- mean(window(b32.correct$y,0,500),na.rm=TRUE)
  
  
  polygon(c((b32.correct$x), rev((b32.correct$x))),
          c((b32.correct$ymax-mv), rev((b32.correct$ymin-mv))),
          col=col.transparent, border=NA)  # Semi-transparent red shading
  
  # Add lines for corrected data
  lines(sr(b32.correct$yaltitude), col=colors[2], lwd=2)   
  lines(sr(b32.correct$y), col=colors[3], lwd=2)    
  
  
  
legend("topright", bty="n",legend=c(expression(delta^{18}*O), "slope effect", "slope effect unc.","elevation effect"), 
       col=c(colors[[1]], colors[[3]], col.transparent,col=colors[[2]]), lty=c(1,1,NA), lwd=c(2,2,NA), fill=c(NA, NA, col.transparent),cex=1)

dev.off()

##### Now for accumulation



b31.correct.accum<-Translate2Climate(b31.upstream,slope2climate=c(-5.68,-5.6,-5.52))
b32.correct.accum<-Translate2Climate(b32.upstream,slope2climate=c(-5.68,-5.6,-5.52))

pdf(file = "./plots/SI_Figure_right.accum.pdf", width = 3.75*SC, height = 4.2*SC)

par(mfcol=c(2,1))
#reduce the lower margin
par(mai=c(0.8,1,0.8,0.8))


accum.obs <- b31.accum
accum.obs.pTs <- ts(rev(c(accum.obs)), start=rev(1950 - c(time(accum.obs)))[1])
plot(sr(accum.obs.pTs),main="B31",ylab = ylab.accum,xlab=xlab.time,col=colors[1],xlim=c(0, 3000))
lines(sr(b31.correct.accum$y),col=colors[[3]],lwd=2)


mv <- mean(window(b31.correct.accum$y,0,500),na.rm=TRUE)
polygon(c((b31.correct.accum$x), rev((b31.correct.accum$x))),
        c((b31.correct.accum$ymax-mv), rev((b31.correct.accum$ymin-mv))),
        col=col.transparent, border=NA)  # Semi-transparent red shading



accum.obs <- b32.accum
accum.obs.pTs <- ts(rev(c(accum.obs)), start=rev(1950 - c(time(accum.obs)))[1])
plot(sr(accum.obs.pTs),main="B32",ylab=ylab.accum,xlab=xlab.time,col=colors[1],xlim=c(0, 3000))
lines(sr(b32.correct.accum$y),lwd=2,col=colors[3])


mv <- mean(window(b32.correct.accum$y,0,500),na.rm=TRUE)
polygon(c((b32.correct.accum$x), rev((b32.correct.accum$x))),
        c((b32.correct.accum$ymax)-mv, rev((b32.correct.accum$ymin))-mv),
        col=col.transparent, border=NA)  # Semi-transparent red shading

legend("topright", bty="n",legend=c("reconstr. accumulation", "slope effect", "slope effect unc."), 
       col=c(colors[[1]], colors[[3]], col.transparent), lty=c(1,1,NA), lwd=c(2,2,NA), fill=c(NA, NA, col.transparent),cex=1)


dev.off()
