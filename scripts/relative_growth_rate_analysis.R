## Analyse of relative growth rate, as relative to time 0 (1987)
# Output analysis
rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(tmap)
library(tmaptools)
library(ggallin)
library(rnaturalearth)
library(ggrepel)
library(ggridges)
## Analyse summaries
file.list <-list.files(path = "data/summary_models/", full.names = T)
## Get the species names
species <- gsub("^summary_|\\.rds$","",list.files(path = "data/summary_models/"))
## Load the route names
routes <- sf::read_sf("data/routes_shp/routes_selected.shp") %>% arrange(routes) %>% pull(routes)
for(i in 1:length(file.list)){
  assign(species[i],
         readRDS(file.list[i]))
  # print(i)
}
## Params to work with
params <- c('beta.time',
            'beta.temp',
            'beta.sky[2]',
            'beta.sky[3]',
            'beta.sky[4]',
            'beta.sky[5]',
            'beta.sky[6]',
            'beta.sky[7]',
            'beta.wind[2]',
            'beta.wind[3]',
            'beta.wind[4]',
            'beta.wind[5]',
            'beta.wind[6]',
            'beta.wind[7]')

## Let's check the Rhats
rhats <- data.frame(species = as.character(), param = as.character(), Rhat = as.numeric())
for(sp in species){
  d <- get(sp)
  for(param in params){
    tmp <- d[rownames(d) == param,]
    rhats <- rbind(rhats,
                   c(sp, param, as.numeric(tmp$Rhat)))

  }
  # print(sp)
}
colnames(rhats) <- c("species", "param", "rhat")
rhats$rhat <- as.numeric(rhats$rhat)

### Take off the species with high values of rhats ####
ok_sp <-
  rhats %>% filter(!grepl("mean.", .$param)) %>%
  group_by(species) %>%
  filter(all(rhat <= 3)) %>%
  distinct(species) %>% pull()

#### Working with the raw data ######
## Modifying the summary to have matrix of N, S, R
d_list <- list()
for(sp in species){
  d <- get(sp)
  ## take off sepcies for which sd > mean for params > 10
  if(nrow(d[(abs(d$mean) > 10 & (abs(d$mean) < abs(d$sd))) ,]) != 0){
    next
  }else{
    ## abundance
    N <- d %>%
      filter(grepl("N", rownames(d))) %>%
      pull(mean) %>%
      matrix(nrow = 1033, ncol = 35)

    d_list[[which(species == sp)]] <- list("N"=N)
    names(d_list)[which(species == sp)] <- sp
  }
  print(sp)
}
## Keep only the species for which rhat are decent
d_list <- d_list[names(d_list) %in% ok_sp]
## Cleaning
rm(list = ls()[ls() %in% species])
gc()

## Now compute the overall abundance
overall_ab <- matrix(nrow = 1033, ncol = 35)
for(i in 1:length(d_list)){
  overall_ab <- ifelse(is.na(overall_ab), 0, overall_ab)
  overall_ab <- overall_ab + d_list[[i]][["N"]]
}
## Add the year as colname and the road
colnames(overall_ab) <- 1987:2021
overall_ab <- cbind.data.frame(routes, overall_ab)

# Rates from T0
overall_abPercent <- matrix(nrow = nrow(overall_ab), ncol = ncol(overall_ab))
for(i in 1:nrow(overall_ab)){
  for(j in 1:(ncol(overall_ab)-1)){
    overall_abPercent[i,j+1] <- ((overall_ab[i,j+1] - overall_ab[i,2])/overall_ab[i,2])
  }
}
colnames(overall_abPercent) <- c("route", 1987:2021)
overall_abPercent <- cbind(routes, overall_abPercent[,-1])

## Compute the trend g to t0
gT0_trend <-
  overall_abPercent %>%
  as.data.frame() %>%
  select(-c("1987")) %>%
  pivot_longer(cols = -routes,
               names_to = "year",
               values_to = "ab") %>%
  mutate(year = as.numeric(year),
         ab = as.numeric(ab)) %>%
  group_by(routes) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(ab ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })

## Add to shapefile
routes_shp <- st_read("data/routes_shp/routes_selected_lines.shp")
routes_shp <-
  routes_shp %>%
  left_join(gT0_trend %>% select(routes,slope) %>% rename(gT0_trend = slope), by = "routes")
## get the states polygon
states_shp <-
  ne_states(returnclass = "sf") %>%
  filter(adm0_a3 == "CAN" | adm0_a3 == "USA" | adm0_a3 == "MEX") %>%
  filter(woe_name != "Hawaii") %>%
  st_crop(st_bbox(c(xmin = -165,
                    xmax = -60,
                    ymax = 62,
                    ymin = 16)))

## Fit GAM
## Add lat and long of route centroid
routes_shp <- routes_shp %>%
  mutate(lon = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
         lat = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(Y))
## k is the number of basis functions used the build the spline. This is somehow the degree of freedom
gT0_gp <- mgcv::gam(gT0_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                    data = routes_shp)
## add the smoothed values to the shp
routes_shp$gT0_trend_gam <- gT0_gp$fitted.values

## Plots
pdf("figures/gT0.pdf",
    width = 5.83, height = 4.13)

gT0_trend %>%
  ggplot()+
  geom_histogram(aes(slope), bins = 40)+
  geom_vline(aes(xintercept = median(slope)), color="red")+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_bw()

overall_abPercent %>%
  as.data.frame() %>%
  pivot_longer(cols = -routes, names_to = "year", values_to = "gT0") %>%
  mutate(year = as.numeric(year), gT0 = as.numeric(gT0)) %>%
  ggplot()+
  geom_line(aes(year, gT0, group = routes), color = "grey")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  stat_summary(aes(year, gT0), fun.y = median, geom = "line", color = "blue",
               linewidth = 1, linetype = 2)+
  theme_light()

routes_shp %>%
  ggplot()+
  geom_histogram(aes(gT0_trend_gam), bins = 40)+
  geom_vline(aes(xintercept = median(gT0_trend_gam)), color="red")+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_bw()

dev.off()

## Maps
pdf("figures/gT0_maps.pdf",
    width=11, height=8.5)

print(
  routes_shp %>%
    st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
    st_centroid() %>%
    ggplot()+
    geom_sf(data = states_shp)+
    geom_sf(aes(color = gT0_trend), size = 5, show.legend = FALSE)+
    scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
    geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
    ggtitle("gT0")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18))
)

plot <-
  routes_shp %>%
  st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
  st_centroid() %>%
  ggplot()+
  geom_sf(data = states_shp)+
  geom_sf(aes(color = gT0_trend), size = 5)+
  scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
  geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
  ggtitle("gT0")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18))

grid::grid.newpage()
grid::grid.draw(cowplot::get_legend(plot))

print(
  routes_shp %>%
    st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
    st_centroid() %>%
    ggplot()+
    geom_sf(data = states_shp)+
    geom_sf(aes(color = gT0_trend_gam), size = 5, show.legend = FALSE)+
    scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
    geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
    ggtitle("gT0_gam")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18))
)

plot <-
  routes_shp %>%
  st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
  st_centroid() %>%
  ggplot()+
  geom_sf(data = states_shp)+
  geom_sf(aes(color = gT0_trend_gam), size = 5)+
  scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
  geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
  ggtitle("gT0_gam")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18))

grid::grid.newpage()
grid::grid.draw(cowplot::get_legend(plot))

dev.off()
