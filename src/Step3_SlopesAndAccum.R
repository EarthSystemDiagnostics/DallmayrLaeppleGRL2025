## Aim: Extract the slopes along the main line of the traverse and read in accumulation
## Sensitivity study on the spatial scale of the slope vs. accumulation 
## For Figure 2:  Plot Slopes vs. Accum
## For Figure 2: Plot Slopes vs. Isotopes (along the main line of the traverse)
## SI Figure 1
## Get the relation of the accumulation rate and d18O (only for the compensation of the isotope to accum relationship)
## Status/belongs to: Topography hypothesis manuscript
## Next step: 


library(dplyr)
library(ggplot2)
library(sf)
library(stars)
library(RColorBrewer)
library(zoo)
library(geosphere)
library(stringr)


basedrive="/Users/tlaepple/data/" 
setwd(paste(basedrive,"topohypothesis",sep=""))


load("./save/data.dat",verbose=TRUE) #Isotope and coordinate data




sl_proxy <-  stars::read_stars("./data/topography/slope_250km_5kmTransect_20pts.tif", proxy = T)
dir_proxy <- read_stars("./data/topography/dir_200m.tif", proxy = T)
altitude_proxy <- stars::read_stars("./data/topography/REMA_200m_dem_filled.tif",proxy=T)



## Extract the slope every 50 points of the radar data (about 0.5km) and interpolate to the resolution of the radar data

slopeExtract <-  function(p, dm_proxy, dir_proxy, transect_distance = 5000, sample = 20, sf = TRUE) {
  require(geosphere)
  start <- p
  pProj <- start %>% st_transform(st_crs(dir_proxy))
  
  p_lonlat <- st_transform(p, crs = 4326)
  
  # Extract the coordinates as a numeric vector
  p_coords <- as.numeric(st_coordinates(p_lonlat))
  
  
  transect <- st_linestring(
    matrix(c(destPoint(p_coords,  st_extract(dir_proxy, pProj)[[1]],               transect_distance/2),
             destPoint(p_coords, (st_extract(dir_proxy, pProj)[[1]] + 180) %% 360, transect_distance/2)), ncol = 2, byrow = TRUE)) %>%
    st_sfc(crs = 4326) %>% st_transform(st_crs(dm_proxy)) %>%
    st_sample(size = sample, type = "regular") %>% st_cast("POINT") %>% st_sf()
  
  if(sf) {
    transect %>% mutate(dm    = st_extract(dm_proxy, transect)[[1]],
                        dir   = st_extract(dir_proxy, transect)[[1]],
                        slope = summary(lm(dm ~ seq(0, transect_distance, length = sample)))$coefficients[2,1]) %>%
      dplyr::select(dm, dir, slope)
  } else {
    summary(lm(st_extract(dm_proxy, transect)[[1]] ~ c(1:nrow(transect))))$coefficients[2,1]
  }
}



GetSlopes<-function(tdistance)
{
  # Ensure the first and last points are included
  indices <- c(1, seq(50, nrow(radar_sf) - 50, by = 50), nrow(radar_sf))
  
  # Extract slopes at the selected points
  slopes_sampled <- sapply(indices, function(i) {
    point <- radar_sf[i, ]
    slopeExtract(point, altitude_proxy, dir_proxy,transect_distance=tdistance)$slope[1]
  })
  
  # Perform linear interpolation
  interp <- approx(x = indices, y = slopes_sampled, xout = 1:nrow(radar_sf), rule = 2)
  
  # Interpolated slopes with same length as original
  interp$y
}



#Read the Accumulation (Radar) data from https://doi.pangaea.de/10.1594/PANGAEA.935129
radar <- read.table("./data/accumulation/Dallmayr-etal_2021_SMB.tab",
                   skip = 23, header = TRUE, sep = "\t") %>%
  mutate(
    x=Distance..km.,
    accum=SMB..kg.m..2.a.,
    lat=Latitude,
    lon=Longitude 
  ) %>%
  select(x,accum,lat,lon) 

# Extract slope values

radar_sf <- radar %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%  st_transform(radar_sf, crs = st_crs(sl_proxy)) 

extracted_values <- stars::st_extract(sl_proxy, radar_sf)

# Bind extracted values to the original radar data
radar <- cbind(radar, extracted_values) 
radar$slope <- radar$slope_250km_5kmTransect_20pts.tif*1000 #Slope im m/km

 
merged<-data  %>% filter(str_detect(name, "^D\\d+$")) 

# Find the closest x in radar for each x in merged and add the corresponding accum
merged_with_accum <- merged %>%
  rowwise() %>%
  mutate(
    closest_index = which.min(abs(radar$x - x)),
    accum = radar$accum[closest_index]         
  ) %>%
  ungroup() %>%
  select(-closest_index) 



### Plotting preparations
#41 points = ~500m running mean for the figures to minimize aliasing
rm41 <- function(x) return(rollmean(x,41))



colors <- brewer.pal(4,"Set1")
colors[1]<-"black"
ylab.accum = expression("Accumulation [kg" ~ m^-2 ~ yr^-1 ~"]")
SC=2


#Figure 2a 

pdf(file="./plots/Fig2a.pdf",width=4.5*SC,height=3*SC)
par(mfcol=c(1,1))
par(mai=c(1,0.8,0.8,0.8))
plot(rm41(radar$x),rm41(radar$slope),type="l",xlab="Distance from B32 [km]",ylab="Slope [m/km]",axes=FALSE,ylim=c(3.5,-1.7),lwd=2)
axis(1)
axis(2)

par(new=TRUE)

plot(rm41(radar$x),rm41(radar$accum),col=colors[3],ylim=c(33,72),type="l",lwd=2,axes=FALSE, xlab="", ylab="")
axis(4,col=colors[3],col.axis=colors[3])  # Add an axis on the right side
box()
mtext(ylab.accum, side=4, line=3, col="black") 
dev.off()

### Figure 2b
#Figure 2a 
pdf(file="./plots/Fig2b.pdf",width=4.5*SC,height=3*SC)
par(mai=c(1,0.8,0.8,0.8))
plot(rm41(radar$x),rm41(radar$slope),type="l",xlab="Distance from B32 [km]",ylab="Slope [m/km]",axes=FALSE,ylim=c(3.5,-1.7),lwd=2)
axis(1)
axis(2)
box()

par(new=TRUE)
plot(merged_with_accum$x,merged_with_accum$d18O,type="p",pch=19,col="red",axes=FALSE, xlab="", ylab="",ylim=c(-43.5,-46.5))
arrows(merged_with_accum$x, merged_with_accum$d18O + 2 * merged_with_accum$d18O.sde, merged_with_accum$x, merged_with_accum$d18O - 2 * merged_with_accum$d18O.sde,
       angle = 90, code = 3, length = 0.05, col = "red")

spline_interpolation <- spline(merged_with_accum$x, merged_with_accum$d18O, n = 500) # Increase n for smoother interpolation

# Add the interpolated spline to the plot
lines(spline_interpolation$x, spline_interpolation$y, col = "red", lwd = 2)


axis(4,col="red",col.axis="red") # Add an axis on the right side

mtext(ylab.accum, side=4, line=3, col="black") 

dev.off()


cor.test(radar$slope,radar$accum)$estimate #0.69
summary(lm(radar$accum~radar$slope)) #-5.6kg/yr+- (0.039) per slope (m/km)
summary(lm(merged_with_accum$accum~merged_with_accum$d18O)) #-5.8kg/yr per permille 



## Sensitivity experiment for the SI; Extract the slopes on a 500m to 12km scale
slopes<-list()
tdistance<-seq(from=500,to=12000,by=500)
for (i in 1:length(tdistance)) {
  slopes[[i]]<-GetSlopes(tdistance=tdistance[i])
  print(i)
}

## Figure S1
pdf(file="./plots/FigS1.pdf",width=6,height=6)
save<-vector()
for (i in 1:length(slopes)) save[i]<-cor(radar$accum,slopes[[i]])
plot(tdistance/1000,save,type="b",pch=19,ylab="Correlation of Slope with Accumulation",xlab="Length Used for Fitting Topographic Slope (km)")
dev.off()


