#Table 1: Sample label, longitude, latitude,  for samples on the main track, Kilometers from B31, Elevation, Slope, Slope SE



library(dplyr)


library(flextable)
library(officer)

dataForTable <- data  %>%
  mutate(distance = x)  %>% mutate(
    slope = sprintf("%.2f", slope),
    slope.sd = sprintf("%.2f", slope.sd),
    lon = sprintf("%.3f", lon),
    lat = sprintf("%.3f", lat),
    distance = round(distance*10)/10,
    altitude = sprintf("%.0f", altitude),
  ) %>% select(name,lon,lat,distance,altitude,slope,slope.sd)

# Convert tibble to a flextable
ft <- flextable(dataForTable)



ft <- ft %>%
   set_header_labels(
    name = "Site Name",
    lon = "Longitude (degE)",
    lat = "Latitude (degN)",
    distance = "Distance from B31 (km)",
    altitude = "Elevation (masl)",
    slope = "Slope (m/km)",
    slope.sd = "SE Slope (m/km)"
  ) %>%
  
  autofit() %>%  # Automatically adjust column widths
  theme_vanilla()  # Apply a nice pre-defined theme

# View the table in RStudio Viewer (optional)
ft


ft <- ft %>%
  set_table_properties(width = 1, layout = "autofit")  # Scale the table to 100% of the page width

# Create a Word document
doc <- read_docx()

# Add the table to the document
doc <- doc %>%
  body_add_flextable(ft) %>%
  body_add_par("")  # Add an empty paragraph for spacing

# Save the Word document
print(doc, target = "table_output.docx")


### Second table



#Table 1:
Sample Label, #Snow sampling tool samples #Mean samples #Liner #Total Nr. #d18O, d18O sd, d18O se

