## Aim: Prepare the isotope data
## Status/belongs to: Topography hypothesis manuscript
## Result: 
## Next step: 
basedrive="/Users/tlaepple/data/" 

library(dplyr)
library(ggplot2)
library(sf)
library(stars)
setwd(paste(basedrive,"topohypothesis",sep=""))

load(file="./save/coords.dat",verbose=TRUE)

#Load Snow sampling tool samples
#Samples have 2-10 locations per site and three depths per location
temp<-read.csv("./data/isotope/KohnenQK_Isotopes_PANGAEA_R_20250120.csv",sep=",",dec=".") %>% mutate(position=Position,name = paste(Site,Site.number,sep="")) %>%
  select(name,position,d18O)
#First step, summarize across depths to d18O.location

samplingtool <- temp %>%
    group_by(name,position) %>%
  summarise(d18OAvg = mean(d18O,na.rm=TRUE))


#Read the Liner data from Hirsch et al., PANGAEA
temp <- read.table("./data/isotope/2018_Kohnen_isotopes.tab",
                    skip = 44, header = TRUE, sep = "\t") %>%
  mutate(
    depth = Depth.ice.snow.top..m.,
    d18O = δ18O.H2O....SMOW.,
    position = Position..m.,
    name = Loc.ID
  ) %>%
  select(depth, d18O, position, name) 


# Interpolate d18O over depth on 1cm resolution before averaging
interpolate_profile <- function(data) {
  interpolated <- approx(
    x = data$depth,
    y = data$d18O,
    xout = seq(min(data$depth), max(data$depth), by = 0.01),
    method = "linear"
  )
  
  
  tibble(
    depth = interpolated$x,
    d18O = interpolated$y
  )
}

# Group by Name and Position, interpolate, and calculate averages
liner <- temp %>%
  group_by(name, position) %>%
  group_modify(~ interpolate_profile(.x)) %>%
  ungroup() %>%
  group_by(name, position) %>%
  summarize(d18OAvg = mean(d18O, na.rm = TRUE), .groups = "drop")

combined <- bind_rows(
  mutate(liner, type = "liner"),
  mutate(samplingtool, type = "samplingtool")
)



#Second step, take means, standard deviation and N across locations
data.linerANDremi <- combined  %>%
  group_by(name) %>%
  summarise(d18O = mean(d18OAvg, na.rm = TRUE),d18O.sd = sd(d18OAvg,na.rm=TRUE),d18O.n=sum(!is.na(d18OAvg))) %>%
  mutate(d18O.sde = d18O.sd/sqrt(d18O.n))

data.sde.mean <- mean(data.linerANDremi$d18O.sde)
data.sd.mean <- mean(data.linerANDremi$d18O.sd)

#Load mean data
#Mean data has two replicates of three depths averaged from 10 locations per site
temp<-read.csv("./data/isotope/KohnenQK_Isotopes_PANGAEA_M_20250120.csv",sep=",",dec=".")

data.mean <- temp %>% mutate(name = paste(Site,Site.number,sep=""))  %>%
  group_by(name) %>%
  dplyr::summarise(d18O = mean(d18O,na.rm=TRUE),d18O.sde=data.sd.mean/sqrt(10))

#Sensitivity test of using 0.9 instead of 1.2m
#data.mean09 <- temp %>% filter(depth.bottom <= 0.9) %>% mutate(name = paste(Site,Site.number,sep=""))  %>%
#  group_by(name) %>%
#  dplyr::summarise(d18O = mean(d18O,na.rm=TRUE),d18O.sde=data.remi.sd.mean/sqrt(10))

#plot(data.mean$d18O,data.mean09$d18O)
#t.test(data.mean$d18O,data.mean09$d18O)
#mean(data.mean$d18O)-mean(data.mean09$d18O) #0.07 permille difference


data<-right_join(coord,rbind(data.linerANDremi %>% select(name,d18O,d18O.sde),data.mean),by="name") %>% select(name,lon,lat,slope,slope.mean,slope.sd,altitude,d18O,d18O.sde)


save(data,file="./save/data.dat")
