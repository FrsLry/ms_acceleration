## The following script contains the code to reproduce the spatial analysis of our paper.
## The data used is the output of the Dail-Madsen model
#######################################

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
library(mgcv)
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

## Plot the distributions of Rhat for each parameter
jpeg("figures/rhats.jpg",
     width = 8.27, height = 5.83, units = "in", res = 500)

ggplot(data = rhats)+
  geom_histogram(aes(x=rhat), bins =  100)+
  facet_wrap(vars(param), scales = "free")+
  # ggtitle("Distributions of Rhat values for the different parameters and for the 564 species")+
  theme_bw()

dev.off()
### Take off the species with rhats => 3 ####
ok_sp <-
  rhats %>% filter(!grepl("mean.", .$param)) %>%
  group_by(species) %>%
  filter(all(rhat <= 3)) %>%
  distinct(species) %>% pull()

#### Working with the raw data ######
## Modifying the output of the DM model to have matrices of N, S, R and L
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
      # floor() %>%
      matrix(nrow = 1033, ncol = 35)
    ## survival
    S <-
      d %>%
      filter(grepl("S", rownames(d))) %>%
      pull(mean) %>%
      # floor() %>%
      ## 34 because no data for year 1
      matrix(nrow = 1033, ncol = 34)
    ## recruitment
    R_raw <- d %>%
      filter(grepl("R", rownames(d))) %>%
      pull(mean) %>%
      matrix(nrow = 1033, ncol = 34)
    ## extinction from survival and abundance
    L_raw <- N[,-ncol(N)] - S
    # Abundance
    N <- cbind(N[,1], S+R_raw)

    d_list[[which(species == sp)]] <- list("N"=N,
                                           "S"=S,"R_raw"=R_raw,"L_raw"=L_raw)
    names(d_list)[which(species == sp)] <- sp
  }
  print(sp)
}
## Keep only the species for which rhat are <= 3
d_list <- d_list[names(d_list) %in% ok_sp]
## Memory Cleaning
rm(list = ls()[ls() %in% species])
gc()

## Now compute the overall abundance, survival, recruitment for each year and route
overall_ab <- matrix(nrow = 1033, ncol = 35)
overall_surv <- matrix(nrow = 1033, ncol = 34)
overall_lossRaw <- matrix(nrow = 1033, ncol = 34)
overall_recRaw <- matrix(nrow = 1033, ncol = 34)
for(i in 1:length(d_list)){
  overall_ab <- ifelse(is.na(overall_ab), 0, overall_ab)
  overall_ab <- overall_ab + d_list[[i]][["N"]]
  overall_surv <- ifelse(is.na(overall_surv), 0, overall_surv)
  overall_surv <- overall_surv + d_list[[i]][["S"]]
  overall_lossRaw <- ifelse(is.na(overall_lossRaw), 0, overall_lossRaw)
  overall_lossRaw <- overall_lossRaw + d_list[[i]][["L_raw"]]
  overall_recRaw <- ifelse(is.na(overall_recRaw), 0, overall_recRaw)
  overall_recRaw <- overall_recRaw + d_list[[i]][["R_raw"]]
}
## Add the year as colname and the road
colnames(overall_ab) <- 1987:2021
overall_ab <- cbind.data.frame(routes, overall_ab)
colnames(overall_surv) <- 1988:2021
overall_surv <- cbind.data.frame(routes, overall_surv)
colnames(overall_lossRaw) <- 1988:2021
overall_lossRaw <- cbind.data.frame(routes, overall_lossRaw)
colnames(overall_recRaw) <- 1988:2021
overall_recRaw <- cbind.data.frame(routes, overall_recRaw)

## Compute the trend
#Abundance
ab_trend <-
  overall_ab %>%
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
# Survival
surv_trend <-
  overall_surv %>%
  pivot_longer(cols = -routes,
               names_to = "year",
               values_to = "surv") %>%
  mutate(year = as.numeric(year),
         surv = as.numeric(surv)) %>%
  group_by(routes) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(surv ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })
# Number of Recruitment
recRaw_trend <-
  overall_recRaw %>%
  pivot_longer(cols = -routes,
               names_to = "year",
               values_to = "recRaw") %>%
  mutate(year = as.numeric(year),
         recRaw = as.numeric(recRaw)) %>%
  group_by(routes) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(recRaw ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })
# Number of Extinction
lossRaw_trend <-
  overall_lossRaw %>%
  pivot_longer(cols = -routes,
               names_to = "year",
               values_to = "lossRaw") %>%
  mutate(year = as.numeric(year),
         lossRaw = as.numeric(lossRaw)) %>%
  group_by(routes) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(lossRaw ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })
## Plot all the abundance trends trends together
pdf("figures/trends_numbers.pdf",
    width = 5.83, height = 4.13)
overall_ab %>%
  pivot_longer(cols = -routes, names_to = "year", values_to = "ab") %>%
  mutate(year = as.numeric(year), ab = as.numeric(ab)) %>%
  ggplot()+
  geom_line(aes(year, ab, group = routes), color = "grey")+
  stat_summary(aes(year, ab), fun.y = median, geom = "line", color = "blue",
               linewidth = 1, linetype = 2)+
  theme_light()

overall_lossRaw %>%
  pivot_longer(cols = -routes, names_to = "year", values_to = "loss_number") %>%
  mutate(year = as.numeric(year), loss_number = as.numeric(loss_number)) %>%
  ggplot()+
  scale_y_continuous(trans=pseudolog10_trans)+
  # ylim(0, 500)+
  geom_line(aes(year, loss_number, group = routes), color = "grey")+
  stat_summary(aes(year, loss_number), fun.y = median, geom = "line", color = "blue",
               linewidth = 1, linetype = 2)+
  theme_light()

overall_recRaw %>%
  pivot_longer(cols = -routes, names_to = "year", values_to = "rec_number") %>%
  mutate(year = as.numeric(year), rec_number = as.numeric(rec_number)) %>%
  ggplot()+
  scale_y_continuous(trans=pseudolog10_trans)+
  geom_line(aes(year, rec_number, group = routes), color = "grey")+
  stat_summary(aes(year, rec_number), fun.y = median, geom = "line", color = "blue",
               linewidth = 1, linetype = 2)+
  theme_light()

## frequencies
## Plot histograms
for(metrics in c("abundance", "loss_number", "recruitment_number")){
  print(ab_trend %>% mutate(metric = "abundance") %>%
          rbind(surv_trend %>% mutate(metric = "survival")) %>%
          rbind(lossRaw_trend %>% mutate(metric = "loss_number")) %>%
          rbind(recRaw_trend %>% mutate(metric = "recruitment_number")) %>%
          select(routes, slope, metric) %>%
          group_by(metric) %>%
          mutate(median = mean(slope)) %>%
          ungroup() %>%
          filter(metric == metrics) %>%
          ggplot()+
          geom_histogram(aes(slope), binwidth = ifelse(metrics == "abundance", 5, 0.2))+
          geom_vline(aes(xintercept = median), color="red")+
          ylab(metrics)+
          # facet_wrap(~metric, scales = "free")+
          geom_vline(aes(xintercept = 0), linetype = "dashed")+
          theme_bw())
}

dev.off()
######

# Growth rate, Extinction and Recruitment rate
overall_abPercent <- matrix(nrow = nrow(overall_ab), ncol = ncol(overall_ab))
overall_lossPercent <- matrix(nrow = nrow(overall_lossRaw), ncol = ncol(overall_lossRaw))
overall_survPercent <- matrix(nrow = nrow(overall_surv), ncol = ncol(overall_surv))
overall_recPercent <- matrix(nrow = nrow(overall_recRaw), ncol = ncol(overall_recRaw))
for(i in 1:nrow(overall_lossRaw)){
  for(j in 2:ncol(overall_lossRaw)){
    overall_abPercent[i,j+1] <- ((overall_ab[i,j+1] - overall_ab[i,j])/overall_ab[i,j])
    overall_lossPercent[i,j] <- (overall_lossRaw[i,j]/overall_ab[i,j])
    overall_survPercent[i,j] <- (overall_surv[i,j]/overall_ab[i,j])
    overall_recPercent[i,j] <- (overall_recRaw[i,j]/overall_ab[i,j])
  }
}
colnames(overall_abPercent) <- c("route", 1987:2021)
overall_abPercent <- cbind(routes, overall_abPercent[,-1])
colnames(overall_lossPercent) <- c("route", 1988:2021)
overall_lossPercent <- cbind(routes, overall_lossPercent[,-1])
overall_lossPercent <- as.data.frame(overall_lossPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
colnames(overall_survPercent) <- c("route", 1988:2021)
overall_survPercent <- cbind(routes, overall_survPercent[,-1])
overall_survPercent <- as.data.frame(overall_survPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
colnames(overall_recPercent) <- c("route", 1988:2021)
overall_recPercent <- cbind(routes, overall_recPercent[,-1])
overall_recPercent <- as.data.frame(overall_recPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))

## Compute the trends for percentages
#Abundance
abPercent_trend <-
  as.data.frame(overall_abPercent) %>%
  select(-"1987") %>% # because NAs
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
# rec ratio
recPercent_trend <-
  overall_recPercent %>%
  pivot_longer(cols = -routes,
               names_to = "year",
               values_to = "rec_percent") %>%
  mutate(year = as.numeric(year)) %>%
  group_by(routes) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(rec_percent ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })
# ext ratio
lossPercent_trend <-
  overall_lossPercent %>%
  pivot_longer(cols = -routes,
               names_to = "year",
               values_to = "loss_percent") %>%
  mutate(year = as.numeric(year)) %>%
  group_by(routes) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(loss_percent ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })

## Plot the ab, rec, ext and surv ratio
pdf("figures/trends_ratio.pdf",
    width = 5.83, height = 4.13)

overall_abPercent %>%
  as.data.frame() %>%
  pivot_longer(cols = -routes, names_to = "year", values_to = "growth_rate") %>%
  mutate(year = as.numeric(year), growth_rate = as.numeric(growth_rate)) %>%
  ggplot()+
  geom_line(aes(year, growth_rate, group = routes), color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  scale_y_continuous(trans= ssqrt_trans)+
  stat_summary(aes(year, growth_rate), fun.y = median, geom = "line", color = "blue",
               linewidth = 1, linetype = 2)+
  theme_light()

overall_lossPercent %>%
  pivot_longer(cols = -routes, names_to = "year", values_to = "loss_ratio") %>%
  mutate(year = as.numeric(year), ext_ratio = as.numeric(loss_ratio)) %>%
  ggplot()+
  geom_line(aes(year, loss_ratio, group = routes), color = "grey")+
  scale_y_continuous(trans= "log10")+
  stat_summary(aes(year, loss_ratio), fun.y = median, geom = "line", color = "blue",
               linewidth = 1, linetype = 2)+
  # geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  scale_y_continuous(trans= "log10")+
  theme_light()

overall_recPercent %>%
  pivot_longer(cols = -routes, names_to = "year", values_to = "rec_ratio") %>%
  mutate(year = as.numeric(year), rec_ratio = as.numeric(rec_ratio)) %>%
  ggplot()+
  geom_line(aes(year, rec_ratio, group = routes), color = "grey")+
  scale_y_continuous(trans= "log10")+
  stat_summary(aes(year, rec_ratio), fun.y = median, geom = "line", color = "blue",
               linewidth = 1, linetype = 2)+
  theme_light()

## histograms
for(metrics in c("growth rate","recruitment rate","loss rate")){
  print(abPercent_trend %>% mutate(metric = "growth rate") %>%
          rbind(recPercent_trend %>% mutate(metric = "recruitment rate")) %>%
          rbind(lossPercent_trend %>% mutate(metric = "loss rate")) %>%
          group_by(metric) %>%
          mutate(median = mean(slope)) %>%
          ungroup() %>%
          filter(metric == metrics) %>%
          ggplot()+
          geom_histogram(aes(slope), binwidth = 0.5*10^(-4))+
          geom_vline(aes(xintercept = median), color="red")+
          ylab(metrics)+
          # facet_wrap(~metric, scales = "free")+
          geom_vline(aes(xintercept = 0), linetype = "dashed")+
          theme_bw())
}

dev.off()


#### Load shapefile #####
routes_shp <- st_read("data/routes_shp/routes_selected_lines.shp")
proj <- "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs"
# Smoothing maps using GAMs
## Add lat and long of route centroid
routes_shp <- routes_shp %>%
  mutate(lon = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
         lat = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(Y))
## Join the slopes of ab, rec, ext, surv to the shapefile
routes_shp <- routes_shp %>%
  left_join(ab_trend %>% select(routes,slope) %>% rename(ab_trend = slope), by = "routes") %>%
  left_join(lossRaw_trend %>% select(routes,slope) %>% rename(loss_trend = slope), by = "routes") %>%
  left_join(recRaw_trend %>% select(routes,slope) %>% rename(rec_trend = slope), by = "routes")
## Join the slopes of ab, rec, ext rates to the shapefile
routes_shp <-
  routes_shp %>%
  left_join(abPercent_trend %>% select(routes,slope) %>% rename(growthRate_trend = slope), by = "routes") %>%
  left_join(lossPercent_trend %>% select(routes,slope) %>% rename(lossRate_trend = slope), by = "routes") %>%
  left_join(recPercent_trend %>% select(routes,slope) %>% rename(recRate_trend = slope), by = "routes")

## k is the number of basis functions used the build the spline. This is somehow the degree of freedom
ab_gp <- gam(ab_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
             data = routes_shp)
rec_gp <- gam(rec_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
              data = routes_shp)
loss_gp <- gam(loss_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
               data = routes_shp)
growthRate_gp <- gam(growthRate_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                     data = routes_shp)
recRate_gp <- gam(recRate_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                  data = routes_shp)
lossRate_gp <- gam(lossRate_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                   data = routes_shp)

## add the smoothed values to the shp
routes_shp$ab_trend_gam <- ab_gp$fitted.values
routes_shp$rec_trend_gam <- rec_gp$fitted.values
routes_shp$loss_trend_gam <- loss_gp$fitted.values

## add the smoothed values to the shp
routes_shp$growthRate_trend_gam <- growthRate_gp$fitted.values
routes_shp$recRate_trend_gam <- recRate_gp$fitted.values
routes_shp$lossRate_trend_gam <- lossRate_gp$fitted.values

## add the rec vs loss hue
routes_shp <- routes_shp %>% mutate(rec_vs_loss = case_when(
  (rec_trend > 0 & loss_trend > 0 & abs(rec_trend) > abs(loss_trend)) ~ "#b3df8afa",
  (rec_trend > 0 & loss_trend > 0 & abs(rec_trend) < abs(loss_trend)) ~ "#f5a582fa",
  (rec_trend < 0 & loss_trend > 0 & abs(rec_trend) < abs(loss_trend)) ~ "#ca0020fa",
  (rec_trend < 0 & loss_trend > 0 & abs(rec_trend) > abs(loss_trend)) ~ "#e66101fa",
  (rec_trend < 0 & loss_trend < 0 & abs(rec_trend) > abs(loss_trend)) ~ "#feb963fa",
  (rec_trend < 0 & loss_trend < 0 & abs(rec_trend) < abs(loss_trend)) ~ "#a6cee3fa",
  (rec_trend > 0 & loss_trend < 0 & abs(rec_trend) < abs(loss_trend)) ~ "#1f78b5fa",
  (rec_trend > 0 & loss_trend < 0 & abs(rec_trend) > abs(loss_trend)) ~ "#33a02cfa"),
  recRate_vs_lossRate = case_when(
    (recRate_trend > 0 & lossRate_trend > 0 & abs(recRate_trend) > abs(lossRate_trend)) ~ "#b3df8afa",
    (recRate_trend > 0 & lossRate_trend > 0 & abs(recRate_trend) < abs(lossRate_trend)) ~ "#f5a582fa",
    (recRate_trend < 0 & lossRate_trend > 0 & abs(recRate_trend) < abs(lossRate_trend)) ~ "#ca0020fa",
    (recRate_trend < 0 & lossRate_trend > 0 & abs(recRate_trend) > abs(lossRate_trend)) ~ "#e66101fa",
    (recRate_trend < 0 & lossRate_trend < 0 & abs(recRate_trend) > abs(lossRate_trend)) ~ "#feb963fa",
    (recRate_trend < 0 & lossRate_trend < 0 & abs(recRate_trend) < abs(lossRate_trend)) ~ "#a6cee3fa",
    (recRate_trend > 0 & lossRate_trend < 0 & abs(recRate_trend) < abs(lossRate_trend)) ~ "#1f78b5fa",
    (recRate_trend > 0 & lossRate_trend < 0 & abs(recRate_trend) > abs(lossRate_trend)) ~ "#33a02cfa"),
  rec_vs_loss_gam = case_when(
    (rec_trend_gam > 0 & loss_trend_gam > 0 & abs(rec_trend_gam) > abs(loss_trend_gam)) ~ "#b3df8afa",
    (rec_trend_gam > 0 & loss_trend_gam > 0 & abs(rec_trend_gam) < abs(loss_trend_gam)) ~ "#f5a582fa",
    (rec_trend_gam < 0 & loss_trend_gam > 0 & abs(rec_trend_gam) < abs(loss_trend_gam)) ~ "#ca0020fa",
    (rec_trend_gam < 0 & loss_trend_gam > 0 & abs(rec_trend_gam) > abs(loss_trend_gam)) ~ "#e66101fa",
    (rec_trend_gam < 0 & loss_trend_gam < 0 & abs(rec_trend_gam) > abs(loss_trend_gam)) ~ "#feb963fa",
    (rec_trend_gam < 0 & loss_trend_gam < 0 & abs(rec_trend_gam) < abs(loss_trend_gam)) ~ "#a6cee3fa",
    (rec_trend_gam > 0 & loss_trend_gam < 0 & abs(rec_trend_gam) < abs(loss_trend_gam)) ~ "#1f78b5fa",
    (rec_trend_gam > 0 & loss_trend_gam < 0 & abs(rec_trend_gam) > abs(loss_trend_gam)) ~ "#33a02cfa"),
  recRate_vs_lossRate_gam = case_when(
    (recRate_trend_gam > 0 & lossRate_trend_gam > 0 & abs(recRate_trend_gam) > abs(lossRate_trend_gam)) ~ "#b3df8afa",
    (recRate_trend_gam > 0 & lossRate_trend_gam > 0 & abs(recRate_trend_gam) < abs(lossRate_trend_gam)) ~ "#f5a582fa",
    (recRate_trend_gam < 0 & lossRate_trend_gam > 0 & abs(recRate_trend_gam) < abs(lossRate_trend_gam)) ~ "#ca0020fa",
    (recRate_trend_gam < 0 & lossRate_trend_gam > 0 & abs(recRate_trend_gam) > abs(lossRate_trend_gam)) ~ "#e66101fa",
    (recRate_trend_gam < 0 & lossRate_trend_gam < 0 & abs(recRate_trend_gam) > abs(lossRate_trend_gam)) ~ "#feb963fa",
    (recRate_trend_gam < 0 & lossRate_trend_gam < 0 & abs(recRate_trend_gam) < abs(lossRate_trend_gam)) ~ "#a6cee3fa",
    (recRate_trend_gam > 0 & lossRate_trend_gam < 0 & abs(recRate_trend_gam) < abs(lossRate_trend_gam)) ~ "#1f78b5fa",
    (recRate_trend_gam > 0 & lossRate_trend_gam < 0 & abs(recRate_trend_gam) > abs(lossRate_trend_gam)) ~ "#33a02cfa"))

## histograms for gams
pdf("figures/hist_gam.pdf",
    width = 5.83, height = 4.13)
for(metrics in c("ab_trend_gam", "rec_trend_gam", "loss_trend_gam",
                 "growthRate_trend_gam", "recRate_trend_gam", "lossRate_trend_gam")){
  print(
    routes_shp %>% rename(var = metrics) %>%
      ggplot()+
      geom_histogram(aes(var), bins = 40)+
      geom_vline(aes(xintercept = mean(var)), color="red")+
      ylab(metrics)+
      geom_vline(aes(xintercept = 0), linetype = "dashed")+
      theme_bw()
  )
}
dev.off()

##########
## get the states polygon
states_shp <-
  ne_states(returnclass = "sf") %>%
  filter(adm0_a3 == "CAN" | adm0_a3 == "USA" | adm0_a3 == "MEX") %>%
  filter(woe_name != "Hawaii") %>%
  st_crop(st_bbox(c(xmin = -165,
                    xmax = -60,
                    ymax = 62,
                    ymin = 16)))

pdf("figures/trend_maps_gam.pdf",
    width=11, height=8.5)
for(metric in c("ab_trend_gam",
                "rec_trend_gam",
                "loss_trend_gam",
                "growthRate_trend_gam",
                "recRate_trend_gam",
                "lossRate_trend_gam")){
  print(
    routes_shp %>%
      st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
      st_centroid() %>%
      rename(variable = metric) %>%
      ggplot()+
      geom_sf(data = states_shp)+
      geom_sf(aes(color = variable), size = 5, show.legend = FALSE)+
      scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
      geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
      ggtitle(metric)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 18))
  )

  plot <-
    routes_shp %>%
    st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
    st_centroid() %>%
    rename(variable = metric) %>%
    ggplot()+
    geom_sf(data = states_shp)+
    geom_sf(aes(color = variable), size = 5)+
    scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695",
                          labels = scales::scientific_format())+
    geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
    ggtitle(metric)+
    labs(colour = metric)+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18))

  grid::grid.newpage()
  grid::grid.draw(cowplot::get_legend(plot))

}

# for(metric in c("rec_vs_loss_gam", "recRate_vs_lossRate_gam")){
#   tmp <- routes_shp %>%
#     st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
#     st_centroid() %>%
#     rename(variable = metric) %>%
#     mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))
#   print(tmp %>%
#           ggplot()+
#           geom_sf(data = states_shp)+
#           geom_sf(aes(color = variable), size = 5, show.legend = F)+
#           scale_color_manual(values = unique(sort(tmp$variable)))+
#           geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
#           ggtitle(metric)+
#           theme_bw()+
#           # labs(colour = paste0("log(",metric,")"))+
#           theme(plot.title = element_text(hjust = 0.5),
#                 text = element_text(size = 18)))
# }
dev.off()

pdf("figures/trend_maps.pdf",
    width=11, height=8.5)
for(metric in c("ab_trend",
                "rec_trend",
                "loss_trend",
                "growthRate_trend",
                "recRate_trend",
                "lossRate_trend")){
  print(
    routes_shp %>%
      st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
      st_centroid() %>%
      rename(variable = metric) %>%
      ggplot()+
      geom_sf(data = states_shp)+
      geom_sf(aes(color = variable), size = 5, show.legend = F)+
      scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
      geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
      ggtitle(metric)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 18))
  )

  plot <- routes_shp %>%
    st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
    st_centroid() %>%
    rename(variable = metric) %>%
    ggplot()+
    geom_sf(data = states_shp)+
    geom_sf(aes(color = variable), size = 5)+
    scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
    geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
    ggtitle(metric)+
    theme_bw()+
    labs(colour = metric)+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18))

  grid::grid.newpage()
  grid::grid.draw(cowplot::get_legend(plot))

}

for(metric in c("rec_vs_loss", "recRate_vs_lossRate")){
  tmp <- routes_shp %>%
    st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
    st_centroid() %>%
    rename(variable = metric) %>%
    mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))
  print(tmp %>%
          ggplot()+
          geom_sf(data = states_shp)+
          geom_sf(aes(color = variable, shape = nSign), size = 5, show.legend = F)+
          scale_color_manual(values = unique(sort(tmp$variable)))+
          geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
          ggtitle(metric)+
          theme_bw()+
          # labs(colour = paste0("log(",metric,")"))+
          theme(plot.title = element_text(hjust = 0.5),
                text = element_text(size = 18)))
}
dev.off()

####### Rec vs Loss
svg("figures/rec_vs_loss_notsmooth.svg",
    height = 4.13, width = 5.83)
my_palette <- c("#1f78b5fa",
                "#33a02cfa",
                "#a6cee3fa",
                "#b3df8afa",
                "#feb963fa",
                "#f5a582fa",
                "#e66101fa",
                "#ca0020fa")
tmp <-
  data.frame(lossRate_slope = routes_shp$lossRate_trend,
             recRate_slope = routes_shp$recRate_trend,
             growthRate_slope = routes_shp$growthRate_trend,
             abundance_slope = routes_shp$ab_trend) %>%
  mutate(color = case_when(
    (recRate_slope > 0 & lossRate_slope > 0 & abs(recRate_slope) > abs(lossRate_slope)) ~ "#b3df8afa",
    (recRate_slope > 0 & lossRate_slope > 0 & abs(recRate_slope) < abs(lossRate_slope)) ~ "#f5a582fa",
    (recRate_slope < 0 & lossRate_slope > 0 & abs(recRate_slope) < abs(lossRate_slope)) ~ "#ca0020fa",
    (recRate_slope < 0 & lossRate_slope > 0 & abs(recRate_slope) > abs(lossRate_slope)) ~ "#e66101fa",
    (recRate_slope < 0 & lossRate_slope < 0 & abs(recRate_slope) > abs(lossRate_slope)) ~ "#feb963fa",
    (recRate_slope < 0 & lossRate_slope < 0 & abs(recRate_slope) < abs(lossRate_slope)) ~ "#a6cee3fa",
    (recRate_slope > 0 & lossRate_slope < 0 & abs(recRate_slope) < abs(lossRate_slope)) ~ "#1f78b5fa",
    (recRate_slope > 0 & lossRate_slope < 0 & abs(recRate_slope) > abs(lossRate_slope)) ~ "#33a02cfa"
  )) %>%
  mutate(nSign = ifelse(abundance_slope < 0, "decrease", "increase"))
tmp %>%
  arrange(color) %>%
  ggplot()+
  geom_point(aes(lossRate_slope, recRate_slope, color = color),
             shape = ifelse(tmp$nSign == "increase", "\u2191", "\u2193"),
             size = 4, show.legend = F)+
  scale_color_manual(values = unique(sort(tmp$color)))+
  xlim(-max(abs(routes_shp$lossRate_trend)), max(abs(routes_shp$lossRate_trend)))+
  ylim(-max(abs(routes_shp$recRate_trend)), max(abs(routes_shp$recRate_trend)))+
  scale_x_continuous(labels = scales::scientific)+
  geom_hline(yintercept = 0, linetype="dashed", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", size = 1)+
  geom_abline(slope = 1, linetype = "longdash")+
  geom_abline(slope = -1, linetype = "longdash")+
  theme_bw()
dev.off()

svg("figures/rec_vs_loss_smooth.svg",
    height = 4.13, width = 5.83)
tmp <-
  data.frame(lossRate_slope_gam = routes_shp$lossRate_trend_gam,
             recRate_slope_gam = routes_shp$recRate_trend_gam,
             growthRate_slope_gam = routes_shp$growthRate_trend_gam,
             abundance_slope_gam = routes_shp$ab_trend_gam) %>%
  mutate(color = case_when(
    (recRate_slope_gam > 0 & lossRate_slope_gam > 0 & abs(recRate_slope_gam) > abs(lossRate_slope_gam)) ~ "#b3df8afa",
    (recRate_slope_gam > 0 & lossRate_slope_gam > 0 & abs(recRate_slope_gam) < abs(lossRate_slope_gam)) ~ "#f5a582fa",
    (recRate_slope_gam < 0 & lossRate_slope_gam > 0 & abs(recRate_slope_gam) < abs(lossRate_slope_gam)) ~ "#ca0020fa",
    (recRate_slope_gam < 0 & lossRate_slope_gam > 0 & abs(recRate_slope_gam) > abs(lossRate_slope_gam)) ~ "#e66101fa",
    (recRate_slope_gam < 0 & lossRate_slope_gam < 0 & abs(recRate_slope_gam) > abs(lossRate_slope_gam)) ~ "#feb963fa",
    (recRate_slope_gam < 0 & lossRate_slope_gam < 0 & abs(recRate_slope_gam) < abs(lossRate_slope_gam)) ~ "#a6cee3fa",
    (recRate_slope_gam > 0 & lossRate_slope_gam < 0 & abs(recRate_slope_gam) < abs(lossRate_slope_gam)) ~ "#1f78b5fa",
    (recRate_slope_gam > 0 & lossRate_slope_gam < 0 & abs(recRate_slope_gam) > abs(lossRate_slope_gam)) ~ "#33a02cfa"
  )) %>%
  mutate(nSign = ifelse(abundance_slope_gam < 0, "decrease", "increase"))
tmp %>%
  arrange(color) %>%
  ggplot()+
  geom_point(aes(lossRate_slope_gam, recRate_slope_gam, color = color),
             shape = ifelse(tmp$nSign == "increase", "\u2191", "\u2193"),
             size = 4, show.legend = F)+
  scale_color_manual(values = unique(sort(tmp$color)))+
  xlim(-max(abs(routes_shp$lossRate_trend_gam)), max(abs(routes_shp$lossRate_trend_gam)))+
  ylim(-max(abs(routes_shp$recRate_trend_gam)), max(abs(routes_shp$recRate_trend_gam)))+
  geom_hline(yintercept = 0, linetype="dashed", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", size = 1)+
  geom_abline(slope = 1, linetype = "longdash")+
  geom_abline(slope = -1, linetype = "longdash")+
  theme_bw()
dev.off()

### per species change at the North-American level (i.e. all the roads together)
result_list <- list()
for(sp in names(d_list)){
  print(sp)
  for(m in names(d_list[[sp]])){
    # print(m)
    if(m == "N"){
      # result <- c(sp, m, (lm(colSums(d_list[[sp]][[m]][,-1]) ~ c(1988:2021)))$coefficients[2])
      result <- c(sp, m, (lm(colSums(d_list[[sp]][[m]]) ~ c(1987:2021)))$coefficients[2])
    }else{
      result <- c(sp, m, (lm(colSums(d_list[[sp]][[m]]) ~ c(1988:2021)))$coefficients[2])
    }
    result_list[[length(result_list) + 1]] <- result
  }
}

# Combine the results into a data frame
per_sp_overall <- do.call(rbind, result_list)
colnames(per_sp_overall) <- c("species", "metric", "slope")

## Compute per species Rates at the North American level (all the roads together)
## For the per species rates, we discard the first values. This is because some species have N close to 0 at time 1, leading to extreme values of growth, rec and loss rates
result_list <- list()
for(sp in names(d_list)){
  N0 <- colSums(d_list[[sp]][["N"]][,-ncol(d_list[[sp]][["N"]])])
  N1 <- colSums(d_list[[sp]][["N"]][,-1])
  R <- colSums(d_list[[sp]][["R_raw"]])
  L <- colSums(d_list[[sp]][["L_raw"]])

  N_rate <- (N1 - N0)/N0
  N_rate[is.nan(N_rate)] <- 0
  N_rate[N_rate == "Inf"] <- 0
  R_rate <- R/N0
  R_rate[is.nan(R_rate)] <- 0
  R_rate[R_rate == "Inf"] <- 0
  L_rate <- L/N0
  L_rate[is.nan(L_rate)] <- 0
  L_rate[L_rate == "Inf"] <- 0

  N_rate <- N_rate[-1]
  R_rate <- R_rate[-1]
  L_rate <- L_rate[-1]

  for(m in c("N_rate","R_rate","L_rate")){
    result_list[[length(result_list) + 1]] <- c(sp, m, (lm(get(m) ~ c(1989:2021)))$coefficients[2])
  }
}

# Combine the results into a data frame
per_sp_rates <- do.call(rbind, result_list)
colnames(per_sp_rates) <- c("species", "metric", "slope")
per_sp_overall <- rbind(per_sp_overall, per_sp_rates) %>% as.data.frame() %>% mutate(slope = as.numeric(slope))
rm(per_sp_rates)

# Create the dataset used to plot the per species time series
per_sp_timeseries <- data.frame(species = as.character(),
                                metric = as.character(),
                                year = as.numeric(),
                                value = as.numeric())
## Add the habitat
hab <- readRDS("data/selected_species.rds")
## Numbers
pdf("figures/perSpecies.pdf",
    height = 4.13, width = 5.83)
for(m in c("N","R_raw","L_raw",
           "N_rate", "R_rate", "L_rate")){
  print(
    per_sp_overall %>%
      left_join(hab %>% rename(species = COM_NAME_updated), by = "species") %>%
      filter(metric == m) %>%
      arrange(slope) %>%
      filter(row_number() <= 20 | row_number() > n() - 20) %>%
      ggplot()+
      geom_col(aes(reorder(species, -slope), slope))+
      scale_y_continuous(labels = scales::scientific)+
      ggtitle(m)+
      coord_flip()+
      theme_minimal()+
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.y = element_blank())
  )
}

## Compute the slopes
for(sp in names(d_list)){
  for(metric in c("N", "R_raw", "L_raw")){
    per_sp_timeseries <-
      rbind(per_sp_timeseries,
            data.frame(species = sp,
                       metric = metric,
                       year = if(metric == "N") as.numeric(1987:2021) else(as.numeric(1988:2021)),
                       value = as.numeric(colSums(d_list[[sp]][[metric]]))
            ))
  }
  per_sp_timeseries <-
    rbind(per_sp_timeseries,
          data.frame(species = sp,
                     metric = "N_rate",
                     year = 1989:2021,
                     value = (colSums(d_list[[sp]][["N"]][,-1])[-1] - colSums(d_list[[sp]][["N"]][,-1])[-length(colSums(d_list[[sp]][["N"]][,-1]))])/colSums(d_list[[sp]][["N"]][,-1])[-length(colSums(d_list[[sp]][["N"]][,-1]))]
          ))
  per_sp_timeseries <-
    rbind(per_sp_timeseries,
          data.frame(species = sp,
                     metric = "R_rate",
                     year = 1989:2021,
                     value = colSums(d_list[[sp]][["R_raw"]][,-1])/colSums(d_list[[sp]][["N"]][,-1])[-length(colSums(d_list[[sp]][["N"]][,-1]))]
          ))
  per_sp_timeseries <-
    rbind(per_sp_timeseries,
          data.frame(species = sp,
                     metric = "L_rate",
                     year = 1989:2021,
                     value = colSums(d_list[[sp]][["L_raw"]][,-1])/colSums(d_list[[sp]][["N"]][,-1])[-length(colSums(d_list[[sp]][["N"]][,-1]))]
          ))

}
## Make the per species time series
for(m in unique(per_sp_timeseries$metric)){
  print(
    per_sp_timeseries %>%
      filter(metric == m) %>%
      ggplot()+
      geom_line(aes(year, value, group = species), color = "grey")+
      stat_summary(aes(year, value), fun.y = median, geom = "line", color = "blue",
                   linewidth = 1, linetype = 2)+
      ylab(m)+
      scale_y_continuous(trans= ssqrt_trans)+
      geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
      theme_light()
  )
}
## Make the per_species histograms
for(m in unique(per_sp_overall$metric)){
  print(
    per_sp_overall %>%
      filter(metric == m) %>%
      ggplot()+geom_histogram(aes(slope), bins = 100)+
      geom_vline(aes(xintercept = mean(slope)), color="red")+
      ylab(m)+
      geom_vline(aes(xintercept = 0), linetype = "dashed")+
      scale_x_continuous(trans = ssqrt_trans)+
      theme_bw()
  )
}
dev.off()

svg("figures/rec_vs_loss_species.svg",
    height = 4.13, width = 5.83)
tmp <-
  per_sp_overall %>%
  filter(metric %in% c("L_rate", "R_rate", "N")) %>%
  pivot_wider(names_from = "metric", values_from = "slope") %>%
  mutate(color = case_when(
    (R_rate > 0 & L_rate > 0 & abs(R_rate) > abs(L_rate)) ~ "#b3df8afa",
    (R_rate > 0 & L_rate > 0 & abs(R_rate) < abs(L_rate)) ~ "#f5a582fa",
    (R_rate < 0 & L_rate > 0 & abs(R_rate) < abs(L_rate)) ~ "#ca0020fa",
    (R_rate < 0 & L_rate > 0 & abs(R_rate) > abs(L_rate)) ~ "#e66101fa",
    (R_rate < 0 & L_rate < 0 & abs(R_rate) > abs(L_rate)) ~ "#feb963fa",
    (R_rate < 0 & L_rate < 0 & abs(R_rate) < abs(L_rate)) ~ "#a6cee3fa",
    (R_rate > 0 & L_rate < 0 & abs(R_rate) < abs(L_rate)) ~ "#1f78b5fa",
    (R_rate > 0 & L_rate < 0 & abs(R_rate) > abs(L_rate)) ~ "#33a02cfa"
  )) %>%
  arrange(color) %>%
  mutate(nSign = ifelse(N < 0, "decrease", "increase"))

print(
  ## recRate vs. lossRate per species
  tmp %>%
    ggplot()+
    geom_point(aes(L_rate, R_rate, color = color),
               shape = ifelse(tmp$nSign == "increase", "\u2191", "\u2193"),
               show.legend = F, size = 4)+
    geom_abline(slope = 1, linetype = "longdash", color = "#595959ff")+
    geom_abline(slope = -1, linetype = "longdash", color = "#595959ff")+
    geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
    geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
    geom_text_repel(aes(L_rate, R_rate, label = species), check_overlap = T,
                    color = "black",
                    min.segment.length = unit(0, 'lines'),
                    segment.size = 1,
                    fontface = "bold",
                    size = 3) +
    xlim(-max(abs(per_sp_overall$slope[per_sp_overall$metric == "L_rate"])), max(abs(per_sp_overall$slope[per_sp_overall$metric == "L_rate"])))+
    ylim(-max(abs(per_sp_overall$slope[per_sp_overall$metric == "R_rate"])), max(abs(per_sp_overall$slope[per_sp_overall$metric == "R_rate"])))+
    scale_y_continuous(trans= ssqrt_trans, labels = scales::scientific)+
    scale_x_continuous(trans= ssqrt_trans, labels = scales::scientific)+
    scale_color_manual(values = unique(sort(tmp$color)))+
    theme_bw()
)

dev.off()

pdf("figures/boxplot_perspecies.pdf",
    height = 2, width = 5.83)
tmp %>%
  mutate(order = case_when(
    color == "#e66101fa" ~ 1,
    color == "#feb963fa" ~ 2,
    color == "#ca0020fa" ~ 3,
    color == "#f5a582fa" ~ 4,
    color == "#33a02cfa" ~ 5,
    color == "#b3df8afa" ~ 6,
    color == "#1f78b5fa" ~ 7,
    color == "#a6cee3fa" ~ 8
  )) %>%
  ggplot()+
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  geom_violin(aes(reorder(color, order), N, fill = color), show.legend = F)+
  ggforce::geom_sina(aes(reorder(color, order), N), show.legend = F, alpha = .5, size = .5)+
  geom_boxplot(aes(reorder(color, order), N), fill = NA, show.legend = F, width = 0.1, color = "darkgrey", outlier.shape = NA)+
  scale_fill_manual(values = unique(sort(tmp$color)))+
  scale_y_continuous(trans= ggallin::pseudolog10_trans, breaks = c(-4000,-2000,-10,0,10,2000))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title  = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11))
dev.off()


