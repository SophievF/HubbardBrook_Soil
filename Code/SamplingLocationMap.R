## Hubbard Brook archived soil samples project ##
## Map - Sampling locations ##
## Matthew Monroe and Sophie von Fromm ##
## 2025-03-11 ##

#Load Packages
library(maps)
library(sf)

#Create North America map
map("world", region = c("usa", "canada", "mexico"), xlim = c(-140, -40), 
    ylim = c(20, 55), 
    col = "black", fill = FALSE, lwd = 1.6)
map("state", add=TRUE, col="black")
map("world", region=Provinces, col="black", add=TRUE)
map("state", region="new hampshire", col="gray", fill=TRUE, add=TRUE)
map.scale(x = -75, y = 30, relwidth = 0.2, metric = TRUE, ratio = FALSE, 
          col = "black", cex=0.65)

#Make Hubbard Brook shp files into a plottable map data
HB_border_shape <- "./Data/hbef_boundary/hbef_boundary.shp"
HB_border_data <- st_read(HB_border_shape)
HB_border_data <- st_transform(HB_border_data, crs = 4326)

HB_watershed_shape <- "./Data/hbef_wsheds/hbef_wsheds.shp"
HB_watershed_data <- st_read(HB_watershed_shape)
HB_watershed_data <- st_transform(HB_watershed_data, crs=4326)

watershed6 <- HB_watershed_data[HB_watershed_data$WS == "WS6", ]
watershed6 <- st_transform(watershed6, crs=4326)

HB_hydro_shape <- "./Data/hbef_hydro/hbef_hydro.shp"
HB_hydro_data <- st_read(HB_hydro_shape, crs=26919)
HB_hydro_data <- st_transform(HB_hydro_data, crs=4326)

#Create New Hampshire map
map("state", region="new hampshire", col="black", fill=FALSE, lwd = 1.8)
map("county", region="new hampshire", add=TRUE, col="black")
plot(st_geometry(HB_border_data), add = TRUE, border = "black", col ="gray")
map.scale(x = -72.5, y = 45, relwidth = 0.3, metric = TRUE, ratio = FALSE, 
          col = "black", cex=0.8)

#Create Hubbard Brook map
plot(st_geometry(HB_border_data), col ="white")
plot(st_geometry(HB_watershed_data), add=TRUE, border="black", col="mediumorchid1")
plot(st_geometry(HB_hydro_data), add=TRUE, border="black", col="black")
map.scale(x = -71.723, y = 43.925, relwidth = 0.2, metric = TRUE, ratio = FALSE, 
          col = "black", cex=0.5)
points(x=-71.735321, y=43.947937, col = "blue", pch = 16, cex = 1)
legend("bottom", legend = c("Sampling Location"), col = "blue", pch = 16, 
       pt.cex = 1, bty = "y", text.col = "black", cex=0.75)

#Create watershed map
plot(st_geometry(HB_watershed_data), border="black", col="white", lwd = 1.8)
plot(st_geometry(watershed6), add=TRUE, border="black", col="gray")
plot(st_geometry(HB_hydro_data), add=TRUE, border="black", col="black")
map.scale(x = -71.7227, y = 43.9535, relwidth = 0.17, metric = TRUE, ratio = FALSE, 
          col = "black", cex=0.6)
points(x=-71.735321, y=43.947937, col = "blue", pch = 16, cex = 1.5)
legend("bottom", legend = c("1 (1969-2018)", "2 (1998-2023)"), col = "blue", 
       pch = 16, pt.cex = 1.2, bty = "y", text.col = "black", cex = 0.9, 
       title = "Dataset")

