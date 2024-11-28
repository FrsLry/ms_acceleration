## Figures
library(sf)
library(terra)
library(mgcv)
library(geodata)
library(tidyverse)
library(randomForest)

rm(list = ls())

## Load the shapefile
routes <- st_read("data/routes_shp/routes_selected_lines.shp") %>% arrange(routes)

# bounding_box
bbox <- st_bbox(routes)
# extend the bbox by 2 degrees each direction
bbox["xmin"] <- bbox$xmin-2
bbox["ymin"] <- bbox$ymin-2
bbox["xmax"] <- bbox$xmax+2
bbox["ymax"] <- bbox$ymax+2


# CROPLAND DENSITY AND ITS TEMPORAL CHANGE (2003-2019) -------------------------
# downloaded data are at a 30 sec resolution (~1 km at equator)

crops03 <- geodata::cropland(source="GLAD", year=2003, path="interpreting_patterns/rasters")
crops03 <- terra::crop(crops03, bbox)
crops03 <- terra::aggregate(crops03, fact = 10, fun = sum)

crops19 <- geodata::cropland(source="GLAD", year=2019, path="interpreting_patterns/rasters")
crops19 <- terra::crop(crops19, bbox)
crops19 <- terra::aggregate(crops19, fact = 10, fun = sum)

delta.crops <- crops19 - crops03

plot(crops03, col=map.pal("viridis", 100))
plot(delta.crops, col=map.pal("differences", 100))

# extract the values from the rasters to routes
routes$crops03 <- terra::extract(crops03, routes, fun = sum)[,2]
routes$crops19 <- terra::extract(crops19, routes, fun = sum)[,2]
routes$delta.crops <- routes$crops19 - routes$crops03




# HUMAN FOOTPRINT AND ITS TEMPORAL CHANGE (1993-2009) --------------------------
# downloaded data are at a 30 sec resolution (~1 km at equator)

footprint93 <- geodata::footprint(year = 1993, path = "interpreting_patterns/rasters")
footprint93 <- terra::crop(footprint93, bbox)
footprint93 <- terra::aggregate(footprint93, fact = 10, fun = mean)

footprint09 <- geodata::footprint(year = 2009, path = "interpreting_patterns/rasters")
footprint09 <- terra::crop(footprint09, bbox)
footprint09 <- terra::aggregate(footprint09, fact = 10, fun = mean)

delta.footprint <- footprint09 - footprint93

plot(footprint93, col=map.pal("viridis", 100))
plot(delta.footprint, col=map.pal("differences", 100))
lines(routes)

hist(delta.footprint)

# extract the values from the rasters to routes
routes$footprint93 <- terra::extract(footprint93, routes, fun = mean)[,2]
routes$footprint09 <- terra::extract(footprint09, routes, fun = mean)[,2]
routes$delta.footprint <- terra::extract(delta.footprint, routes, fun = mean)[,2]




# POPULATION DENSITY, GPW dataset ----------------------------------------------

# load the 2.5 min raster with human population density
pop2000 <-rast("interpreting_patterns/rasters/gpw_v4_population_count_rev11_2000_2pt5_min.tif")

# extract the population density values, summing the values of 2.5 min pixels touching a route
routes$pop2000 <- terra::extract(pop2000, routes, fun = sum)[,2]


## Load the models' output
for(file in list.files("mixed_model/outputs/", full.names = T)){
  assign(gsub("^mixed_model/outputs/|\\.rds$","",file),
         readRDS(file = file))
}

## Add values of beta1 (slope) to the shapefile
routes$growth_rate <- overall_g_output$mean$beta1


# ------------------------------------------------------------------------------
# THE GAM MODEL ----------------------------------------------------------------
# extract latitude and longitude of route centroids
coords <- data.frame(st_centroid(routes) %>% st_coordinates() )
names(coords) <- c("lon","lat")
routes <- cbind(routes, coords)

# fit the gam model
growth_rate_gam <- gam(growth_rate ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                       data = routes)
routes$delta_g <- growth_rate_gam$fitted.values
# ------------------------------------------------------------------------------

plot(footprint93 ~ log10(crops03), data = routes)
plot(footprint93 ~ log10(pop2000), data = routes)

ggplot(routes, aes(x = pop2000, y=delta_g)) +
  geom_point() +
  scale_x_log10() +
  geom_hline(yintercept=0) +
  geom_smooth(method='lm', formula= y~x)


ggplot(routes, aes(x = footprint93, y=delta_g)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_smooth(method='lm', formula= y~x)


ggplot(routes, aes(x = footprint09, y=delta_g)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_smooth(method='lm', formula= y~x)


ggplot(routes, aes(x = delta.footprint, y=delta_g)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_smooth(method='lm', formula= y~x)


ggplot(routes, aes(x = crops03, y=delta_g)) +
  geom_point() +
  geom_hline(yintercept=0) +
  scale_x_log10() +
  geom_smooth(method='lm', formula= y~x)


ggplot(routes, aes(x = delta.crops, y=delta_g)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_smooth(method='lm', formula= y~x)

# ------------------------------------------------------------------------------

rf <- randomForest(delta_g ~ footprint93 +
                             delta.footprint +
                             crops03 +
                             delta.crops,
             data = routes)
rf

varImpPlot(rf)












