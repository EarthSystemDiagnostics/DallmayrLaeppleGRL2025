## Aim: Predicted effect as map
## Status/belongs to: Topography hypothesis paper



library(dplyr)
library(ggplot2)
library(sf)
library(stars)
library(deming)
library(zoo)
library(RColorBrewer)
basedrive="/Users/tlaepple/data/" 
setwd(paste(basedrive,"topohypothesis",sep=""))


### prepare maps
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name == "Antarctica") %>% st_geometry() %>%
  st_transform("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% 
  st_buffer(25) %>% st_union()
mask <- world %>% st_buffer(100) %>% st_union() %>% st_buffer(-75000)

### Grid
bufferKohnen <- st_point(c(0, -75 )) %>% st_sfc(crs = 4326) %>% st_transform(st_crs(world)) %>%
  st_buffer(700000) %>% st_intersection(mask)

grid <- bufferKohnen %>% st_as_stars(dx = 5000, dy = 5000) %>% st_as_sf() %>% st_centroid() %>%
  filter(st_intersects(., bufferKohnen, sparse = FALSE))

map_grid <- world %>% st_intersection(bufferKohnen %>% st_buffer(2.25e5) 
                                      %>% st_bbox() %>% st_as_sfc())
           

#Load precomputed timeseries at every gridpoint
load("./data/topography/ts_list.rda",verbose=TRUE)

#data is in 50yr steps; thus 60 = 3kyr
sm <- apply(ts_list, 1, function(x) diff(range(window(x,0,60), na.rm = T))) ### here you can make a selection of sampling depth
sm[is.infinite(sm) ] <- NA

sm_rast  <- grid %>% mutate(vals = (sm)) %>% dplyr::select(vals) %>% 
  st_rasterize(., (st_bbox(map_grid) %>% st_as_stars(values = NA, dx = 5000, dy = 5000))) * 0.42 * 1000  





num_colors <- 9  # Number of custom colors you want
custom_colors <-  rev(brewer.pal(11, "Spectral"))[c(1,5,6,8,10,10,10,10)]

# Define your breaks (bin edges)
breaks <- c(0, 0.25,  0.5, 1, 1.5, 2, 3, 5,10)

# Labels for breaks (one per break)
labels <- c("0", "0.25", "0.5", "0.75", "1", "1.5", "2", "3", "5")

sm_rast[sm_rast>10]<-10

# Define coordinates for B31 and B32 for plotting
coords_sf  <- data.frame(
  lat = c(-75.5815, -75.0025),
  lon = c(-3.4303, -0.0073),
  label = c("B31", "B32")
) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)


# Load the 2K coordinates file 
coord.2k_sf <- read.csv("./data/coords/coord2k_all.csv") %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)


map_crs <- st_crs(map_grid)
coordinates_transformed <- st_transform(coord.2k_sf, crs = map_crs)
coordinates_within_map <- st_intersection(coordinates_transformed, map_grid)


lab.d18O <- expression("peak-peak"~delta^{18}*O~"["*"\u2030"*"]")
# ggplot with corrected binned scale

pdf("./plots/Figure4.pdf",width=6,height=6)
ggplot() +
  geom_sf(data = map_grid) +
  geom_stars(data = sm_rast) +
  geom_sf(data = coordinates_within_map, color = "black", size = 1) +  # Add points
  geom_sf(data = coords_sf, color = "green", size = 1) +  # Add points
  scale_fill_stepsn(
    colors = custom_colors,
    breaks = breaks,       # Define bin edges
    labels = function(x) {x}, # Just display the breaks
    show.limits = TRUE,
    na.value = "transparent",
    name = lab.d18O
  ) +
  xlab('') + ylab('') +
  ggspatial:::annotation_scale() +
  theme_light()
dev.off()
