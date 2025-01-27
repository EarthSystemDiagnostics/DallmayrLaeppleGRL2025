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


#Load Snow sampling tool samples https://doi.pangaea.de/10.1594/PANGAEA.974749 
#Samples have 2-10 locations per site and three depths per location

temp <- read.table("./data/isotope/KohnenQK_Isotopes_R.tab",
                   skip = 34, header = TRUE, sep = "\t") %>%
  mutate(
    depth = Depth.ice.snow..m...top.,
    d18O = δ18O.H2O....SMOW.,
    position = Position,
    site = Site..label.,
    site.nr = Site..number.,
    name = paste(site,site.nr,sep=""),
  ) %>%
  select(depth, d18O, position, name) 




#First step, summarize across depths to d18O.location

samplingtool <- temp %>%
  group_by(name,position) %>%
  summarise(d18OAvg = mean(d18O,na.rm=TRUE))


#Read the Liner data from Hirsch et al., PANGAEA https://doi.pangaea.de/10.1594/PANGAEA.956273

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

data.linerANDremi <- combined %>%
  group_by(name) %>%
  summarise(
    d18O = mean(d18OAvg, na.rm = TRUE),
    d18O.sd = sd(d18OAvg, na.rm = TRUE),
    d18O.n = sum(!is.na(d18OAvg)),
    liner.count = sum(type == "liner", na.rm = TRUE), # Count 'liner' type samples
    samplingtool.count = sum(type == "samplingtool", na.rm = TRUE) # Count 'samplingtool' type samples
  ) %>%
  mutate(d18O.sde = d18O.sd / sqrt(d18O.n))


data.sde.mean <- mean(data.linerANDremi$d18O.sde)
data.sd.mean <- mean(data.linerANDremi$d18O.sd)

#Load mean data https://doi.pangaea.de/10.1594/PANGAEA.974778
#Mean data has two replicates of three depths averaged from 10 locations per site

temp <- read.table("./data/isotope/KohnenQK_Isotopes_M.tab",
                   skip = 34, header = TRUE, sep = "\t") %>%
  mutate(
    depth = Depth.ice.snow..m...top.,
    d18O = δ18O.H2O....SMOW.,
    position = Position,
    site = Site..label.,
    site.nr = Site..number.,
    repl = Repl,
    name = paste(site,site.nr,sep="")
  ) %>%
  select(depth, d18O, position, name,repl) 

### Check the difference between replicate samples
data.mean.r <- temp   %>%
  group_by(name,repl) %>%
  dplyr::summarise(d18O = mean(d18O,na.rm=TRUE))

r1 <- data.mean.r %>% filter(repl==1)
r2 <- data.mean.r %>% filter(repl==2)

sqrt(mean((r1$d18O-r2$d18O)^2))

####

data.mean <- temp   %>%
  group_by(name) %>%
  dplyr::summarise(d18O = mean(d18O,na.rm=TRUE),d18O.sde=data.sd.mean/sqrt(10)) %>% mutate(mean.count=10)

#Sensitivity test of using 0.9 instead of 1.2m
#data.mean09 <- temp %>% filter(depth.bottom <= 0.9) %>% mutate(name = paste(Site,Site.number,sep=""))  %>%
#  group_by(name) %>%
#  dplyr::summarise(d18O = mean(d18O,na.rm=TRUE),d18O.sde=data.remi.sd.mean/sqrt(10))

#plot(data.mean$d18O,data.mean09$d18O)
#t.test(data.mean$d18O,data.mean09$d18O)
#mean(data.mean$d18O)-mean(data.mean09$d18O) #0.07 permille difference

temp <- bind_rows(data.linerANDremi %>% select(name,d18O,d18O.sde,d18O.sd,liner.count,samplingtool.count),data.mean)

data<-right_join(coord,temp,by="name") %>% select(name,lon,lat,slope,slope.sd,altitude,d18O,d18O.sd,d18O.sde,x,liner.count,mean.count,samplingtool.count)

data <- data %>%
  mutate(
    liner.count = replace_na(liner.count, 0),
    mean.count = replace_na(mean.count, 0),
    samplingtool.count = replace_na(samplingtool.count, 0)
  )

save(data,file="./save/data.dat")
