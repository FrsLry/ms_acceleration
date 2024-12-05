## Figures
library(sf)
library(terra)
library(mgcv)
library(geodata)
library(rnaturalearth)
library(tidyverse)
library(randomForest)
library(psych)

rm(list = ls())

## Load the shapefile
routes <- st_read("data/routes_shp/routes_selected_lines.shp") %>% arrange(routes)


## get the states polygon
states <-
  ne_states(returnclass = "sf") %>%
  filter(adm0_a3 == "CAN" | adm0_a3 == "USA" | adm0_a3 == "MEX") %>%
  filter(woe_name != "Hawaii") %>%
  st_crop(st_bbox(c(xmin = -165,
                    xmax = -60,
                    ymax = 62,
                    ymin = 16)))

# LAND COVER CLASS AREA --------------------------------------------------------

grass <- sqrt(terra::rast("interpreting_patterns/rasters/grassland.tif"))
routes$grass <- terra::extract(grass, routes, fun = mean)[,2]

shrubs <- sqrt(terra::rast("interpreting_patterns/rasters/shrubs.tif"))
routes$shrubs <- terra::extract(shrubs, routes, fun = mean)[,2]

trees <- sqrt(terra::rast("interpreting_patterns/rasters/trees.tif"))
routes$trees <- terra::extract(trees, routes, fun = mean)[,2]

built <- sqrt(terra::rast("interpreting_patterns/rasters/built.tif"))
routes$built <- terra::extract(built, routes, fun = mean)[,2]

wetland <- sqrt(terra::rast("interpreting_patterns/rasters/wetland.tif"))
routes$wetland <- terra::extract(wetland, routes, fun = mean)[,2]

water <- sqrt(terra::rast("interpreting_patterns/rasters/water.tif"))
routes$water <- terra::extract(water, routes, fun = mean)[,2]

# NPP (MOD17) ------------------------------------------------------------------
NPP <- terra::rast("interpreting_patterns/rasters/NPP.tif")
routes$NPP <- terra::extract(NPP, routes, fun = mean)[,2]

# WorldClim PRECIPITATION AND TEMPERATURE --------------------------------------
temp <- terra::rast("interpreting_patterns/rasters/temp.tif")
routes$temp <- terra::extract(temp, routes, fun = mean)[,2]

precip <- terra::rast("interpreting_patterns/rasters/precip.tif")
routes$precip <- terra::extract(precip, routes, fun = mean)[,2]

# ELEVATION --------------------------------------------------------------------

elevat <- terra::rast("interpreting_patterns/rasters/elev.tif")
routes$elevat <- terra::extract(elevat, routes, fun = mean)[,2]


# TEMPERATURE TREND (1991-2020) ------------------------------------------------
temp91to20 <- terra::rast("interpreting_patterns/rasters/T_trend.tif")
routes$temp91to20 <- terra::extract(temp91to20, routes, fun = mean)[,2]

# N FERTILIZER USAGE -----------------------------------------------------------
Nfertilizer <- terra::rast("interpreting_patterns/rasters/Nfertilizer.tif")
routes$Nfertilizer <- terra::extract(Nfertilizer, routes, fun = mean)[,2]


# CROPLAND DENSITY AND ITS TEMPORAL CHANGE (2003-2019) -------------------------
crops03 <- sqrt(terra::rast("interpreting_patterns/rasters/crops03.tif"))
crops19 <- sqrt(terra::rast("interpreting_patterns/rasters/crops19.tif"))

routes$crops03 <- terra::extract(crops03, routes, fun = mean)[,2]
routes$crops19 <- terra::extract(crops19, routes, fun = mean)[,2]
routes$crops03to19 <- routes$crops19 - routes$crops03

# NDVI CHANGE (1982-2012) ------------------------------------------------------
NDVI82to12 <- terra::rast("interpreting_patterns/rasters/NDVI82to12.tif")
routes$NDVI82to12 <- terra::extract(NDVI82to12, routes, fun = mean)[,2]

# HUMAN FOOTPRINT AND ITS TEMPORAL CHANGE (1993-2009) --------------------------
footprint93 <- terra::rast("interpreting_patterns/rasters/footprint93.tif")
footprint09 <- terra::rast("interpreting_patterns/rasters/footprint09.tif")
delta.footprint <- footprint09 - footprint93

# extract the values from the rasters to routes
routes$footprint93 <- terra::extract(footprint93, routes, fun = mean)[,2]
routes$footprint09 <- terra::extract(footprint09, routes, fun = mean)[,2]
routes$footprint93to09 <- terra::extract(delta.footprint, routes, fun = mean)[,2]


# POPULATION DENSITY, GPW dataset ----------------------------------------------

# load the 2.5 min raster with human population density
pop2000 <- log10(rast("interpreting_patterns/rasters/pop2000.tif") + 1)

# extract the population density values, summing the values of 2.5 min pixels touching a route
routes$log10pop2000 <- terra::extract(pop2000, routes, fun = mean)[,2]


# ------------------------------------------------------------------------------
# PLOT THE RASTERS

png("interpreting_patterns/rasters.png", width= 3000, height = 3000, res = 200)


par(mfrow=c(5,4))

plot(footprint09 - footprint93, col=map.pal("differences", 100),
     main = "Human footprint change (1993-2009)",
     background = "grey")
lines(states, col = "white")

plot(NDVI82to12, col=map.pal("differences", 100),
     main = "NDVI change (1982-2012)",
     background = "grey")
lines(states, col = "white")

plot(temp91to20, col=map.pal("differences", 100),
     main = "Temperature change (1991-2020)",
     background = "grey")
lines(states, col = "white")

plot(crops19 - crops03, col=map.pal("differences", 100),
     main = "Cropland area change (2003-2019)",
     background = "grey")
lines(states, col = "white")

# ---------

plot(footprint93, col=map.pal("viridis", 100),
     main = "Human footprint (1993)",
     background = "grey")
lines(states, col = "white")

plot(NPP, col=map.pal("viridis", 100),
     main = "mean NPP",
     background = "grey")
lines(states, col = "white")

plot(temp, col=map.pal("viridis", 100),
     main = "mean Annual Temperature",
     background = "grey")
lines(states, col = "white")

plot(crops03, col=map.pal("viridis", 100),
     main = "sqrt(Cropland area in 2003)",
     background = "grey")
lines(states, col = "white")

# -----------

plot(pop2000, col=map.pal("viridis", 100),
     main = "log10 Population (2000)",
     background = "grey")
lines(states, col = "white")

plot(precip, col=map.pal("viridis", 100),
     main = "mean Annual Precipitation",
     background = "grey")
lines(states, col = "white")

plot(elevat, col=map.pal("viridis", 100),
     main = "mean Elevation",
     background = "grey")
lines(states, col = "white")

plot(Nfertilizer, col=map.pal("viridis", 100),
     main = "mean N fertilizer use",
     background = "grey")
lines(states, col = "white")

# -----------

plot(grass, col=map.pal("viridis", 100),
     main = "sqrt(grassland area)",
     background = "grey")
lines(states, col = "white")

plot(shrubs, col=map.pal("viridis", 100),
     main = "sqrt (shrubland area)",
     background = "grey")
lines(states, col = "white")

plot(trees, col=map.pal("viridis", 100),
     main = "sqrt (tree-covered area)",
     background = "grey")
lines(states, col = "white")

plot(built, col=map.pal("viridis", 100),
     main = "sqrt (built-up area)",
     background = "grey")
lines(states, col = "white")

plot(wetland, col=map.pal("viridis", 100),
     main = "sqrt (wetland area)",
     background = "grey")
lines(states, col = "white")

plot(water, col=map.pal("viridis", 100),
     main = "sqrt (water area)",
     background = "grey")
lines(states, col = "white")

# --------------

dev.off()



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

## Add raw (not smoothed) values the dataset
routes$delta_N <- overall_N_output$mean$beta1
routes$delta_g <- overall_g_output$mean$beta1
routes$delta_l <- overall_l_output$mean$beta1
routes$delta_r <- overall_r_output$mean$beta1


# fit the gam models
routes$delta_N_gam <- gam(delta_N ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                          data = routes)$fitted.values
routes$delta_g_gam <- gam(delta_g ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                          data = routes)$fitted.values
routes$delta_l_gam <- gam(delta_l ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                          data = routes)$fitted.values
routes$delta_r_gam <- gam(delta_r ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                          data = routes)$fitted.values

# ------------------------------------------------------------------------------

# preparing the data for analyses
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
                    scale = TRUE,
                    hist.col = "grey",
                    method = "spearman") # note the Spearman coefficient!

dev.off()

# ------------------------------------------------------------------------------
# RANDOM FOREST ANALYSIS
rf.g <-randomForest(delta_g_gam ~ grass + trees + shrubs + water +
                                  NPP + temp +   Nfertilizer + footprint93 +
                                  NDVI82to12 + temp91to20 + footprint93to09 + crops03to19,
                    data = na.omit(dat)) # NOTE THE NA OMIT!!

rf.g
varImpPlot(rf.g)
partialPlot(rf.g, pred.data=na.omit(dat), x.var="Nfertilizer")

plot(na.omit(dat)$delta_g_gam, predict(rf.g))
abline(a = 0, b = 1)


rf.N <-randomForest(delta_N_gam ~ grass + trees + shrubs + water +
                      NPP + temp +   Nfertilizer + footprint93 +
                      NDVI82to12 + temp91to20 + footprint93to09 + crops03to19,
                    data = na.omit(dat)) # NOTE THE NA OMIT!!
rf.N
varImpPlot(rf.N)
