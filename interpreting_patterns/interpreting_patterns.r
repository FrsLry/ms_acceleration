## Figures
library(sf)
library(stars)
library(terra)
library(mgcv)

rm(list = ls())

## Load the shapefile
routes_shp <- st_read("data/routes_shp/routes_selected_lines.shp") %>% arrange(routes)

## Load the models' output
for(file in list.files("mixed_model/outputs/", full.names = T)){
  assign(gsub("^mixed_model/outputs/|\\.rds$","",file),
         readRDS(file = file))
}

## Add values of beta1 (slope) to the shapefile
routes_shp$growth_rate <- overall_g_output$mean$beta1

## Add lat and long of route centroid
routes_shp <- routes_shp %>%
  mutate(lon = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
         lat = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(Y))

# fit the gam model
growth_rate_gam <- gam(growth_rate ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                       data = routes_shp)
routes_shp$delta_g <- growth_rate_gam$fitted.values

# load the 2.5 min raster with human population density
pop2000 <-rast("interpreting_patterns/rasters/gpw_v4_population_count_rev11_2000_2pt5_min.tif")
pop2000star <-st_as_stars(pop2000)

# extract the population density values, summing the values of 2.5 min pixels touching a route
pop2000 <- st_as_sf(st_extract(pop2000star, routes_shp, FUN = sum))
routes_shp$pop2000 <- pop2000$gpw_v4_population_count_rev11_2000_2pt5_min
routes_shp$log10pop2000 <- log10(routes_shp$pop2000)
plot(routes_shp[c("log10pop2000", "delta_g")], lwd=2)

# plot the relationship
ggplot(routes_shp, aes(x = pop2000, y=delta_g)) +
  geom_point() +
  scale_x_log10() +
  geom_hline(yintercept=0) +
  geom_smooth()



