## Figures
library(sf)
library(terra)
library(mgcv)
library(geodata)
library(tidyverse)
library(randomForest)
library(psych)

rm(list = ls())

## Load the shapefile
routes <- st_read("data/routes_shp/routes_selected_lines.shp") %>% arrange(routes)

# LAND COVER CLASS AREA --------------------------------------------------------

grass <- terra::rast("interpreting_patterns/rasters/grassland.tif")
routes$grass <- sqrt(terra::extract(grass, routes, fun = sum)[,2])

shrubs <- terra::rast("interpreting_patterns/rasters/shrubs.tif")
routes$shrubs <- sqrt(terra::extract(shrubs, routes, fun = sum)[,2])

trees <- terra::rast("interpreting_patterns/rasters/trees.tif")
routes$trees <- sqrt(terra::extract(trees, routes, fun = sum)[,2])

built <- terra::rast("interpreting_patterns/rasters/built.tif")
routes$built <- sqrt(terra::extract(built, routes, fun = sum)[,2])

wetland <- terra::rast("interpreting_patterns/rasters/wetland.tif")
routes$wetland <- sqrt(terra::extract(wetland, routes, fun = sum)[,2])

water <- terra::rast("interpreting_patterns/rasters/water.tif")
routes$water <- sqrt(terra::extract(water, routes, fun = sum)[,2])

# NPP (MOD17) ------------------------------------------------------------------
NPP <- terra::rast("interpreting_patterns/rasters/NPP.tif")
routes$NPP <- terra::extract(NPP, routes, fun = mean)[,2]

# WorldClim PRECIPITATION AND TEMPERATURE --------------------------------------
temp <- terra::rast("interpreting_patterns/rasters/temp.tif")
routes$temp <- terra::extract(temp, routes, fun = mean)[,2]

precip <- terra::rast("interpreting_patterns/rasters/precip.tif")
routes$precip <- terra::extract(precip, routes, fun = mean)[,2]

# TEMPERATURE TREND (1991-2020) ------------------------------------------------
temp91to20 <- terra::rast("interpreting_patterns/rasters/T_trend.tif")
routes$temp91to20 <- terra::extract(temp91to20, routes, fun = mean)[,2]

# N FERTILIZER USAGE -----------------------------------------------------------
Nfertilizer <- terra::rast("interpreting_patterns/rasters/Nfertilizer.tif")
routes$Nfertilizer <- terra::extract(Nfertilizer, routes, fun = mean)[,2]


# CROPLAND DENSITY AND ITS TEMPORAL CHANGE (2003-2019) -------------------------
crops03 <- terra::rast("interpreting_patterns/rasters/crops03.tif")
crops19 <- terra::rast("interpreting_patterns/rasters/crops19.tif")

routes$crops03 <- terra::extract(crops03, routes, fun = sum)[,2]
routes$crops19 <- terra::extract(crops19, routes, fun = sum)[,2]
routes$crops03to19 <- routes$crops19 - routes$crops03

# NDVI CHANGE (1982-2012) ------------------------------------------------------
NDVI82to12 <- terra::rast("interpreting_patterns/rasters/NDVI82to12.tif")
routes$NDVI82to12 <- terra::extract(NDVI82to12, routes, fun = mean)[,2]

# HUMAN FOOTPRINT AND ITS TEMPORAL CHANGE (1993-2009) --------------------------
footprint93 <- terra::rast("interpreting_patterns/rasters/footprint93.tif")
footprint09 <- terra::rast("interpreting_patterns/rasters/footprint09.tif")
delta.footprint <- footprint09 - footprint93

plot(footprint93, col=map.pal("viridis", 100))
plot(delta.footprint, col=map.pal("differences", 100))
lines(routes)

# extract the values from the rasters to routes
routes$footprint93 <- terra::extract(footprint93, routes, fun = mean)[,2]
routes$footprint09 <- terra::extract(footprint09, routes, fun = mean)[,2]
routes$footprint93to09 <- terra::extract(delta.footprint, routes, fun = mean)[,2]


# POPULATION DENSITY, GPW dataset ----------------------------------------------

# load the 2.5 min raster with human population density
pop2000 <-rast("interpreting_patterns/rasters/pop2000.tif")

# extract the population density values, summing the values of 2.5 min pixels touching a route
routes$log10pop2000 <- log10(terra::extract(pop2000, routes, fun = sum)[,2] + 1)


# ------------------------------------------------------------------------------
## Load the models' output
for(file in list.files("mixed_model/outputs/", full.names = T)){
  assign(gsub("^mixed_model/outputs/|\\.rds$","",file),
         readRDS(file = file))
}


# ------------------------------------------------------------------------------
# THE GAM MODEL ----------------------------------------------------------------
# extract latitude and longitude of route centroids
coords <- data.frame(st_centroid(routes) %>% st_coordinates() )
names(coords) <- c("lon","lat")
routes <- cbind(routes, coords)

## Add raw (not smoothed) values delta_g to the dataset
routes$delta_g <- overall_g_output$mean$beta1

# fit the gam model
growth_rate_gam <- gam(delta_g ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                       data = routes)
routes$delta_g_gam <- growth_rate_gam$fitted.values

# ------------------------------------------------------------------------------

# preparing the data for correlation matrix plot
data.frame(names(routes))
dat <- as.data.frame(routes[, 24:ncol(routes)])
dat <- select(dat, -c("geometry")) # remove the last "geometry" column

# ------------------------------------------------------------------------------

# correlation matrix plot from the "psych" package

png("interpreting_patterns/all_pairs.png", width= 2500, height = 2500)

psych::pairs.panels(dat,
                    ellipses = FALSE,
                    density = TRUE,
                    smooth = TRUE,
                    ci = TRUE,
                    stars = TRUE,
                    scale = TRUE)

dev.off()






plot(footprint93 ~ log10(crops03), data = routes)
plot(footprint93 ~ log10(pop2000), data = routes)

ggplot(routes, aes(x = pop2000, y=delta_g)) +
  geom_point() +
  scale_x_log10() +
  geom_hline(yintercept=0) +
  geom_smooth(method='lm', formula= y~x)


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

plot(predict(rf), routes$delta_g)
abline(a=0, b=1)












