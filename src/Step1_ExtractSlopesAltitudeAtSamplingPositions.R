
## Aim: Extract the slopes and altitude incl. uncertainty at the sampling grid
## Status/belongs to: Topography hypothesis manuscript
## Result: saves 'topo.dat' in the temporary save folder


library(dplyr)
library(ggplot2)
library(sf)
library(stars)
basedrive="/Users/tlaepple/data/" 
setwd(paste(basedrive,"topohypothesis",sep=""))



## Extract the altitude and slope at all sites
sl_proxy <-  stars::read_stars("./data/topography/slope_250km_5kmTransect_20pts.tif", proxy = T)
altitude_proxy <- stars::read_stars("./data/topography/REMA_200m_dem_filled.tif",proxy=T)

temp<-read.csv("./data/coords/Basispositionen_LatLon_WGS 84.csv")
coord<-tibble(name=temp[,1],lon=temp[,2],lat=temp[,3])

lnd_geom  <- st_as_sf(coord, coords = c("lon", "lat"), crs = 4326) %>% st_transform(st_crs(sl_proxy))

coord$slope <- stars::st_extract(sl_proxy, lnd_geom)[[1]] * 1000

#Also extract the spatial standard deviation
coord$slope.sd<- (st_extract(sl_proxy, lnd_geom %>% st_buffer(1000), FUN = sd))[[1]]*1000


coord$altitude <- stars::st_extract(altitude_proxy , lnd_geom)[[1]] 

save(coord,file="./save/coords.dat")


