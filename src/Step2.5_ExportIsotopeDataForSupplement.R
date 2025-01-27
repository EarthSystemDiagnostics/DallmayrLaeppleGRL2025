#Code to create the table for the supplements
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
    distance = "Distance from B31 (km)*",
    altitude = "Elevation (masl)",
    slope = "Slope (m/km)",
    slope.sd = "sd Slope (m/km)"
  ) %>%
  
  autofit() %>%  # 
  theme_vanilla()  %>%
  set_table_properties(width = 0.8, layout = "autofit")  

### Second table




dataForTable <- data  %>%
  mutate(distance = x)  %>% mutate(
    d18O.sd = sprintf("%.2f", d18O.sd),
    d18O.sde = sprintf("%.2f", d18O.sde),
    d18O = sprintf("%.2f", d18O),
  ) %>% select(name,d18O,d18O.sd,d18O.sde,samplingtool.count,mean.count,liner.count) 



# Convert tibble to a flextable
ft2 <- flextable(dataForTable)


ylab.d18O <- expression(delta^{18}*O~"["*"\u2030"*"]")


ft2 <- ft2 %>%
  set_header_labels(
    name = "Site Name",
    d18O = "mean d18O (‰)",
    d18O.sd = "sd d18O (‰)*",
    d18O.sde = "se d18O (‰)**",
    samplingtool.count = "#R samples",
    mean.count = "#M samples",
    liner.count = "#L samples"

  ) %>%
  
  autofit() %>%  # Automatically adjust column widths
  theme_vanilla() %>%  
 set_table_properties(width = 0.8, layout = "autofit")  # Scale the table to 100% of the page width


# Create a Word document
doc <- read_docx()

# Add the table to the document
doc <- doc %>%
  body_add_flextable(ft) %>%
  body_add_par("")  # Add an empty paragraph for spacing


# Add the table to the document
doc <- doc %>%
  body_add_flextable(ft2) %>%
  body_add_par("")  # Add an empty paragraph for spacing
 

# Save the Word document
print(doc, target = "table_output.docx")

