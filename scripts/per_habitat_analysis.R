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
            'beta.wind[7]',
            'mean.phi',
            'mean.gamma',
            'mean.p',
            'mean.lam')

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

### If you use the models with the covariates, take off the species with high values of rhats ####
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
    # floor() %>%
    ## 34 because no data for year 1
    matrix(nrow = 1033, ncol = 34)
  ## extinction from survival and abundance
  L_raw <- N[,-ncol(N)] - S
  # Abundance
  N <- cbind(N[,1], S+R_raw)
  # print(sp)
  d_list[[which(species == sp)]] <- list("N"=N,
                                         "S"=S,"R_raw"=R_raw,"L_raw"=L_raw)
  names(d_list)[which(species == sp)] <- sp
}
## Keep only the species for which rhat are decent
d_list <- d_list[names(d_list) %in% ok_sp]
## Cleaning
rm(list = ls()[ls() %in% species])
gc()

## Create the lists of each habitat ##########
hab <- readRDS("data/selected_species.rds")
per_habitat <- data.frame(routes = as.character(),
                          year = as.numeric(),
                          value = as.numeric(),
                          metric = as.character(),
                          habitat = as.character())

pdf("figures/per_habitat_spatial.pdf",
    width=11, height=8.5)
## Here, for each habitat type, I sum all the species together at each route
for(habitat in unique(hab$Cornell_Habitat)){
  tmp_list <- d_list[names(d_list) %in% hab$COM_NAME_updated[hab$Cornell_Habitat == habitat]]
  ## Now compute the overall abundance, survival, recruitment for each year and route
  tmp_ab <- matrix(nrow = 1033, ncol = 35)
  tmp_surv <- matrix(nrow = 1033, ncol = 34)
  tmp_lossRaw <- matrix(nrow = 1033, ncol = 34)
  tmp_recRaw <- matrix(nrow = 1033, ncol = 34)
  for(i in 1:length(tmp_list)){
    tmp_ab <- ifelse(is.na(tmp_ab), 0, tmp_ab)
    tmp_ab <- tmp_ab + tmp_list[[i]][["N"]]
    tmp_surv <- ifelse(is.na(tmp_surv), 0, tmp_surv)
    tmp_surv <- tmp_surv + tmp_list[[i]][["S"]]
    tmp_lossRaw <- ifelse(is.na(tmp_lossRaw), 0, tmp_lossRaw)
    tmp_lossRaw <- tmp_lossRaw + tmp_list[[i]][["L_raw"]]
    tmp_recRaw <- ifelse(is.na(tmp_recRaw), 0, tmp_recRaw)
    tmp_recRaw <- tmp_recRaw + tmp_list[[i]][["R_raw"]]
  }
  ## Add the year as colname and the road
  colnames(tmp_ab) <- 1987:2021
  tmp_ab <- cbind.data.frame(routes, tmp_ab)
  colnames(tmp_surv) <- 1988:2021
  tmp_surv <- cbind.data.frame(routes, tmp_surv)
  colnames(tmp_lossRaw) <- 1988:2021
  tmp_lossRaw <- cbind.data.frame(routes, tmp_lossRaw)
  colnames(tmp_recRaw) <- 1988:2021
  tmp_recRaw <- cbind.data.frame(routes, tmp_recRaw)
  ## Plot all the road trends together for each habitat
  print(
    tmp_ab %>%
      pivot_longer(cols = -routes, names_to = "year", values_to = "ab") %>%
      mutate(year = as.numeric(year), ab = as.numeric(ab)) %>%
      ggplot()+
      ggtitle(habitat)+
      geom_line(aes(year, ab, group = routes), color = "grey")+
      stat_summary(aes(year, ab), fun = median, geom = "line", color = "blue", size = 1, linetype = 2)+
      theme_light()
  )
  ## Boxplots with the ext and rec number
  # Compute the lower and upper limits for the y-axis
  y_lower <- quantile(unlist(c(tmp_lossRaw[,-1], tmp_recRaw[,-1])), 0.25) - 1.5 * IQR(unlist(c(tmp_lossRaw[,-1], tmp_recRaw[,-1])))
  y_upper <- quantile(unlist(c(tmp_lossRaw[,-1], tmp_recRaw[,-1])), 0.75) + 1.5 * IQR(unlist(c(tmp_lossRaw[,-1], tmp_recRaw[,-1])))

  print(
    tmp_recRaw %>%
      pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>%
      mutate(year = as.numeric(year), value = as.numeric(value), metric = "rec_number") %>%
      rbind(tmp_lossRaw %>%
              pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>%
              mutate(year = year, value = as.numeric(value), metric = "loss_number")) %>%
      ggplot()+
      ggtitle(habitat)+
      geom_boxplot(aes(year, value, color = metric), outlier.shape = NA)+
      stat_summary(aes(year, value, group = metric), fun = median,
                   geom = "line", linetype = "dashed", size = 1)+
      coord_cartesian(ylim = c(0, y_upper))+
      theme_bw()
  )

  ## Compute the trend
  #Abundance
  ab_trend <-
    tmp_ab %>%
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

  # Number of Recruitment
  recRaw_trend <-
    tmp_recRaw %>%
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
    tmp_lossRaw %>%
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

  # Growth rate, Extinction and Recruitment rate
  tmp_abPercent <- matrix(nrow = nrow(tmp_ab), ncol = ncol(tmp_ab))
  tmp_lossPercent <- matrix(nrow = nrow(tmp_lossRaw), ncol = ncol(tmp_lossRaw))
  tmp_recPercent <- matrix(nrow = nrow(tmp_recRaw), ncol = ncol(tmp_recRaw))
  for(i in 1:nrow(tmp_lossRaw)){
    for(j in 2:ncol(tmp_lossRaw)){
      tmp_abPercent[i,j+1] <- ((tmp_ab[i,j+1] - tmp_ab[i,j])/tmp_ab[i,j])
      tmp_lossPercent[i,j] <- (tmp_lossRaw[i,j]/tmp_ab[i,j])
      tmp_recPercent[i,j] <- (tmp_recRaw[i,j]/tmp_ab[i,j])
    }
  }
  colnames(tmp_abPercent) <- c("route", 1987:2021)
  tmp_abPercent <- cbind(routes, tmp_abPercent[,-1])
  tmp_abPercent <- as.data.frame(tmp_abPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  # tmp_abPercent[is.na(tmp_abPercent)] <- 0
  colnames(tmp_lossPercent) <- c("route", 1988:2021)
  tmp_lossPercent <- cbind(routes, tmp_lossPercent[,-1])
  tmp_lossPercent <- as.data.frame(tmp_lossPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  colnames(tmp_recPercent) <- c("route", 1988:2021)
  tmp_recPercent <- cbind(routes, tmp_recPercent[,-1])
  tmp_recPercent <- as.data.frame(tmp_recPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))

  ## Plot the growth rate, extinction and recruitment ratio
  print(
    tmp_abPercent %>%
      pivot_longer(cols = -routes, names_to = "year", values_to = "growth_rate") %>%
      mutate(year = as.numeric(year), growth_rate = as.numeric(growth_rate)) %>%
      ggplot()+
      ggtitle(habitat)+
      geom_line(aes(year, growth_rate, group = routes), color = "grey")+
      stat_summary(aes(year, growth_rate), fun = median, geom = "line", color = "blue", size = .5, linetype = 2)+
      theme_light()
  )
  ## Boxplots with the ext and rec number
  # Compute the lower and upper limits for the y-axis
  # compute lower and upper whiskers
  y_lower <- quantile(unlist(c(tmp_lossPercent[,-1], tmp_recPercent[,-1])), 0.25, na.rm = TRUE) - 1.5 * IQR(unlist(c(tmp_lossPercent[,-1], tmp_recPercent[,-1])), na.rm = T)
  y_upper <- quantile(unlist(c(tmp_lossPercent[,-1], tmp_recPercent[,-1])), 0.75, na.rm = TRUE) + 1.5 * IQR(unlist(c(tmp_lossPercent[,-1], tmp_recPercent[,-1])), na.rm = T)

  print(
    tmp_recPercent %>%
      pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>%
      mutate(year = as.numeric(year), value = as.numeric(value), metric = "rec_rate") %>%
      rbind(tmp_lossPercent %>%
              pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>%
              mutate(year = year, value = as.numeric(value), metric = "loss_rate")) %>%
      ggplot()+
      ggtitle(habitat)+
      geom_boxplot(aes(year, value, color = metric), outlier.shape = NA)+
      stat_summary(aes(year, value, group = metric), fun = median,
                   geom = "line", linetype = "dashed", size = 1)+
      coord_cartesian(ylim = c(0, y_upper))+
      theme_bw()
  )
  ## Compute the trends for rates
  #Abundance
  abPercent_trend <-
    as.data.frame(tmp_abPercent) %>%
    select(-"1987") %>%
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
    tmp_recPercent %>%
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
    tmp_lossPercent %>%
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


  #### Maps #####
  routes_shp <- st_read("data/routes_shp/routes_selected_lines.shp")
  ## Join the slopes of ab, rec, ext, surv to the shapefile
  routes_shp <- routes_shp %>%
    left_join(ab_trend %>% select(routes,slope) %>% rename(ab_trend = slope), by = "routes") %>%
    left_join(lossRaw_trend %>% select(routes,slope) %>% rename(loss_trend = slope), by = "routes") %>%
    left_join(recRaw_trend %>% select(routes,slope) %>% rename(rec_trend = slope), by = "routes")

  # Smoothing maps using GAMs
  ## Add lat and long of route centroid
  routes_shp <- routes_shp %>%
    mutate(lon = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
           lat = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(Y))

  ## k is the number of basis functions used the build the spline. This is somehow the degree of freedom
  ab_gp <- gam(ab_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
               data = routes_shp)
  rec_gp <- gam(rec_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                data = routes_shp)
  loss_gp <- gam(loss_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                 data = routes_shp)
  ## add the smoothed values to the shp
  routes_shp$ab_trend_gam <- ab_gp$fitted.values
  routes_shp$rec_trend_gam <- rec_gp$fitted.values
  routes_shp$loss_trend_gam <- loss_gp$fitted.values
  ## smoothing for rates
  routes_shp <-
    routes_shp %>%
    left_join(abPercent_trend %>% select(routes,slope) %>% rename(growthRate_trend = slope), by = "routes") %>%
    left_join(lossPercent_trend %>% select(routes,slope) %>% rename(lossRate_trend = slope), by = "routes") %>%
    left_join(recPercent_trend %>% select(routes,slope) %>% rename(recRate_trend = slope), by = "routes")

  ## k is the number of basis functions used the build the spline. This is somehow the degree of freedom
  growthRate_gp <- gam(growthRate_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                       data = routes_shp)
  recRate_gp <- gam(recRate_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                    data = routes_shp)
  lossRate_gp <- gam(lossRate_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                     data = routes_shp)
  ## add the smoothed values to the shp
  routes_shp$growthRate_trend_gam <- growthRate_gp$fitted.values
  routes_shp$recRate_trend_gam <- recRate_gp$fitted.values
  routes_shp$lossRate_trend_gam <- lossRate_gp$fitted.values
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

  for(metric in c("ab_trend",
                  "rec_trend",
                  "loss_trend",
                  "growthRate_trend",
                  "recRate_trend",
                  "lossRate_trend")){
    if(grepl("ab_trend", metric)){
      print(
        routes_shp %>%
          st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
          st_centroid() %>%
          rename(variable = metric) %>%
          mutate(sign = sign(variable),
                 variable = log(abs(variable))) %>%
          mutate(variable = variable*sign) %>%
          ggplot()+
          geom_sf(data = states_shp)+
          geom_sf(aes(color = variable), size = 5, show.legend = F)+
          scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
          geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
          ggtitle(paste0(habitat ," - ","log(",metric,")"))+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5),
                text = element_text(size = 18))
      )
    }else{
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
          ggtitle(paste0(habitat ," - ",metric))+
          theme_bw()+
          # labs(colour = paste0("log(",metric,")"))+
          theme(plot.title = element_text(hjust = 0.5),
                text = element_text(size = 18))
      )
    }
  }

  for(metric in c("ab_trend_gam","rec_trend_gam","loss_trend_gam",
                  "growthRate_trend_gam","recRate_trend_gam","lossRate_trend_gam")){
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
        ggtitle(paste0(habitat ," - ",metric))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size = 18))
    )
  }

  ## Create the dataset to compare all the habitats together
  per_habitat <- rbind(per_habitat,
                       tmp_ab %>% pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>% mutate(metric = "abundance") %>%
                         rbind(tmp_recRaw %>% pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>% mutate(metric = "recruitment")) %>%
                         rbind(tmp_lossRaw %>% pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>% mutate(metric = "loss")) %>%
                         rbind(tmp_abPercent %>% pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>% mutate(metric = "growth rate")) %>%
                         rbind(tmp_recPercent %>% pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>% mutate(metric = "recruitment rate")) %>%
                         rbind(tmp_lossPercent %>% pivot_longer(cols = -routes, names_to = "year", values_to = "value") %>% mutate(metric = "loss rate")) %>%
                         mutate(habitat = paste0(habitat, " (s = ", length(tmp_list), ")")))
}
############

### plots the numbers for all the habitats together
per_habitat %>%
  filter(!grepl("rate$", .$metric)) %>%
  filter(!grepl("abundance", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment","loss","abundance"))) %>%
  ggplot()+
  stat_summary(aes(year, value, group = habitat, color = habitat), fun = median,
               geom = "line",
               # linetype = "dashed",
               size = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_y_continuous(trans=pseudolog10_trans)+
  facet_wrap(~metric, scale="fixed")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

per_habitat %>%
  filter(!grepl("rate$", .$metric)) %>%
  filter(grepl("abundance", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment","loss","abundance"))) %>%
  ggplot()+
  stat_summary(aes(year, value, group = habitat, color = habitat), fun = median,
               geom = "line",
               # linetype = "dashed",
               size = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_y_continuous(trans=pseudolog10_trans)+
  facet_wrap(~metric, scale="fixed")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()


### plots the rates for all the habitats together
per_habitat %>%
  filter(grepl("rate$", .$metric)) %>%
  filter(!grepl("growth rate", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment rate","loss rate","growth rate"))) %>%
  ggplot()+
  stat_summary(aes(year, value, group = habitat, color = habitat), fun = median,
               geom = "line",
               # linetype = "dashed",
               size = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_y_continuous(trans=pseudolog10_trans)+
  facet_wrap(~metric, scale="free")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

per_habitat %>%
  filter(grepl("rate$", .$metric)) %>%
  filter(grepl("growth rate", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment rate","loss rate","growth rate"))) %>%
  ggplot()+
  stat_summary(aes(year, value, group = habitat, color = habitat), fun = median,
               geom = "line",
               # linetype = "dashed",
               size = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_y_continuous(trans=pseudolog10_trans)+
  facet_wrap(~metric, scale="free")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

## Compute trends for each habitat and each route
trends <-
  per_habitat %>%
  group_by(metric, habitat, routes) %>%
  mutate(year = as.numeric(year)) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(value ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })
dev.off()

pdf("figures/perhabitat_allRoad.pdf",
    width = 7, height = 4.13)

per_habitat <- data.frame(year = as.numeric(),
                          value = as.numeric(),
                          metric = as.character(),
                          habitat = as.character())

## per habitat, all routes
for(habitat in unique(hab$Cornell_Habitat)){
  tmp_list <- d_list[names(d_list) %in% hab$COM_NAME_updated[hab$Cornell_Habitat == habitat]]
  ## Now compute the overall abundance, survival, recruitment for each year and route
  tmp_ab <- matrix(nrow = 1033, ncol = 35)
  tmp_lossRaw <- matrix(nrow = 1033, ncol = 34)
  tmp_recRaw <- matrix(nrow = 1033, ncol = 34)
  for(i in 1:length(tmp_list)){
    tmp_ab <- ifelse(is.na(tmp_ab), 0, tmp_ab)
    tmp_ab <- tmp_ab + tmp_list[[i]][["N"]]
    tmp_lossRaw <- ifelse(is.na(tmp_lossRaw), 0, tmp_lossRaw)
    tmp_lossRaw <- tmp_lossRaw + tmp_list[[i]][["L_raw"]]
    tmp_recRaw <- ifelse(is.na(tmp_recRaw), 0, tmp_recRaw)
    tmp_recRaw <- tmp_recRaw + tmp_list[[i]][["R_raw"]]
  }
  ## Add the year as colname and the road
  tmp_ab <- t(colSums(tmp_ab))
  colnames(tmp_ab) <- 1987:2021
  tmp_lossRaw <- t(colSums(tmp_lossRaw))
  colnames(tmp_lossRaw) <- 1988:2021
  tmp_recRaw <- t(colSums(tmp_recRaw))
  colnames(tmp_recRaw) <- 1988:2021

  # Growth rate, Extinction and Recruitment rate
  tmp_abPercent <- matrix(nrow = nrow(tmp_ab), ncol = ncol(tmp_ab))
  tmp_lossPercent <- matrix(nrow = nrow(tmp_lossRaw), ncol = ncol(tmp_lossRaw))
  tmp_recPercent <- matrix(nrow = nrow(tmp_recRaw), ncol = ncol(tmp_recRaw))
  for(i in 1:nrow(tmp_lossRaw)){
    for(j in 1:ncol(tmp_lossRaw)){
      tmp_abPercent[i,j+1] <- ((tmp_ab[i,j+1] - tmp_ab[i,j])/tmp_ab[i,j])
      tmp_lossPercent[i,j] <- (tmp_lossRaw[i,j]/tmp_ab[i,j])
      tmp_recPercent[i,j] <- (tmp_recRaw[i,j]/tmp_ab[i,j])
    }
  }
  colnames(tmp_abPercent) <- c(1987:2021)
  tmp_abPercent <- as.data.frame(tmp_abPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  tmp_abPercent[sapply(tmp_abPercent, is.nan)] <- 0
  tmp_abPercent[tmp_abPercent == "Inf"] <- 0

  colnames(tmp_lossPercent) <- c(1988:2021)
  tmp_lossPercent <- as.data.frame(tmp_lossPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  tmp_lossPercent[sapply(tmp_lossPercent, is.nan)] <- 0
  tmp_lossPercent[tmp_lossPercent == "Inf"] <- 0

  colnames(tmp_recPercent) <- c(1988:2021)
  tmp_recPercent <- as.data.frame(tmp_recPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  tmp_recPercent[sapply(tmp_recPercent, is.nan)] <- 0
  tmp_recPercent[tmp_recPercent == "Inf"] <- 0

  ## Create the dataset to compare all the habitats together
  per_habitat <- rbind(per_habitat,
                       as.data.frame(tmp_ab) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "abundance") %>%
                         rbind(as.data.frame(tmp_recRaw) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "recruitment")) %>%
                         rbind(as.data.frame(tmp_lossRaw) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "loss")) %>%
                         rbind(as.data.frame(tmp_abPercent) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "growth rate")) %>%
                         rbind(as.data.frame(tmp_recPercent) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "recruitment rate")) %>%
                         rbind(as.data.frame(tmp_lossPercent) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "loss rate")) %>%
                         mutate(habitat = paste0(habitat, " (s = ", length(tmp_list), ")")))


}

### plots the numbers for all the families together
per_habitat %>%
  filter(!grepl("rate$", .$metric)) %>%
  filter(grepl("abundance", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment","loss","abundance"))) %>%
  ggplot()+
  geom_line(aes(year, value, color = habitat, group = habitat),
            show.legend = F)+
  scale_y_continuous(labels = scales::scientific)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_y_continuous(trans=pseudolog10_trans)+
  facet_wrap(~metric, scale="fixed")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

per_habitat %>%
  filter(!grepl("rate$", .$metric)) %>%
  filter(!grepl("abundance", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment","loss","abundance"))) %>%
  ggplot()+
  geom_line(aes(year, value, color = habitat, group = habitat),
            show.legend = F)+
  scale_y_continuous(labels = scales::scientific)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_y_continuous(trans=pseudolog10_trans)+
  facet_wrap(~metric, scale="fixed")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

per_habitat %>%
  filter(grepl("rate$", .$metric)) %>%
  filter(grepl("growth rate", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment rate","loss rate","growth rate"))) %>%
  filter(year != 1988) %>%
  ggplot()+
  geom_line(aes(year, value, color = habitat, group = habitat),
            show.legend = F)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_y_continuous(trans=ggallin::ssqrt_trans,
                     labels = scales::scientific)+
  facet_wrap(~metric, scale="free")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

per_habitat %>%
  filter(grepl("rate$", .$metric)) %>%
  filter(!grepl("growth rate", .$metric)) %>%
  mutate(across(metric, factor, levels=c("recruitment rate","loss rate","growth rate"))) %>%
  filter(year != 1988) %>%
  ggplot()+
  geom_line(aes(year, value, color = habitat, group = habitat),
            show.legend = F)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_y_continuous(trans=ggallin::ssqrt_trans,
                     labels = scales::scientific)+
  facet_wrap(~metric, scale="free")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

## Compute trends for each family at the continent scale
trends <-
  per_habitat %>%
  mutate(year = as.numeric(year)) %>%
  group_by(metric, habitat) %>%
  arrange(year, .by_group = T) %>%
  do({
    mod <- lm(value ~ year, data = .)
    data.frame(slope = coef(mod)[2],                             #### extract the slope
               intercept = coef(mod)[1],                         #### extract the intercept
               pval = (summary(mod)$coefficients[2,"Pr(>|t|)"]))
  })

## Make barplots
trends %>%
  filter(metric %in% c("abundance", "loss", "recruitment")) %>%
  # filter(!grepl("n = 1" , family))  %>%
  ggplot()+
  geom_col(aes(reorder(habitat, -slope), slope, fill = habitat),
           show.legend = F)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("Slope")+
  scale_y_continuous(labels = scales::scientific)+
  coord_flip()+
  theme_minimal()+
  facet_wrap(~metric, scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())


trends %>%
  filter(metric %in% c("growth rate", "loss rate", "recruitment rate")) %>%
  mutate(metric = factor(metric, levels = c("growth rate", "loss rate", "recruitment rate"))) %>%
  arrange(metric, -slope) %>%
  mutate(habitat = factor(habitat, levels = unique(.$habitat))) %>%
  ###############
  mutate(fill = as.character(habitat)) %>%
  ggplot()+
  geom_col(aes(habitat, slope, fill = fill),
           show.legend = F)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("Slope")+
  scale_y_continuous(labels = scales::scientific)+
  coord_flip()+
  theme_minimal()+
  facet_wrap(~metric, scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

plot <-
  trends %>%
  filter(metric %in% c("growth rate", "loss rate", "recruitment rate")) %>%
  mutate(metric = factor(metric, levels = c("growth rate", "loss rate", "recruitment rate"))) %>%
  arrange(metric, -slope) %>%
  mutate(habitat = factor(habitat, levels = unique(.$habitat))) %>%
  ###############
  mutate(fill = as.character(habitat)) %>%
  ggplot()+
  geom_col(aes(habitat, slope, fill = fill),
           show.legend = T)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("Slope")+
  scale_y_continuous(labels = scales::scientific)+
  coord_flip()+
  theme_minimal()+
  facet_wrap(~metric, scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

  grid::grid.newpage()
  grid::grid.draw(cowplot::get_legend(plot))

dev.off()

svg("figures/rec_vs_loss_habitat.svg",
    height = 4.13, width = 5.83)

tmp <-
  trends %>%
  filter(metric %in% c("loss rate", "recruitment rate", "abundance")) %>%
  select(metric, habitat, slope) %>%
  pivot_wider(names_from = "metric", values_from = "slope") %>%
  mutate(color = case_when(
    (`recruitment rate` > 0 & `loss rate` > 0 & abs(`recruitment rate`) > abs(`loss rate`)) ~ "#b3df8afa",
    (`recruitment rate` > 0 & `loss rate` > 0 & abs(`recruitment rate`) < abs(`loss rate`)) ~ "#f5a582fa",
    (`recruitment rate` < 0 & `loss rate` > 0 & abs(`recruitment rate`) < abs(`loss rate`)) ~ "#ca0020fa",
    (`recruitment rate` < 0 & `loss rate` > 0 & abs(`recruitment rate`) > abs(`loss rate`)) ~ "#e66101fa",
    (`recruitment rate` < 0 & `loss rate` < 0 & abs(`recruitment rate`) > abs(`loss rate`)) ~ "#feb963fa",
    (`recruitment rate` < 0 & `loss rate` < 0 & abs(`recruitment rate`) < abs(`loss rate`)) ~ "#a6cee3fa",
    (`recruitment rate` > 0 & `loss rate` < 0 & abs(`recruitment rate`) < abs(`loss rate`)) ~ "#1f78b5fa",
    (`recruitment rate` > 0 & `loss rate` < 0 & abs(`recruitment rate`) > abs(`loss rate`)) ~ "#33a02cfa"
  )) %>%
  arrange(color) %>%
  mutate(nSign = ifelse(abundance < 0, "decrease", "increase"))

tmp %>%
  ggplot()+
  geom_point(aes(`loss rate`,`recruitment rate`, color = color),
             shape = ifelse(tmp$nSign == "increase", "\u2191", "\u2193"),
             show.legend = F, size = 6)+
  geom_abline(slope = 1, linetype = "longdash", color = "#595959ff")+
  geom_abline(slope = -1, linetype = "longdash", color = "#595959ff")+
  geom_hline(yintercept = 0, linetype="dashed", size = 1, color = "#595959ff")+
  geom_vline(xintercept = 0, linetype="dashed", size = 1, color = "#595959ff")+
  ggrepel::geom_text_repel(aes(`loss rate`, `recruitment rate`, label = habitat), check_overlap = T,
                           color = "black",
                           min.segment.length = unit(0, 'lines'),
                           force = 1,
                           direction = c("both"),
                           segment.size = 1,
                           fontface = "bold",
                           size = 3)+
  # vjust = "inward", hjust = "inward")+
  xlim(-max(abs(trends$slope[trends$metric == "loss rate"])), max(abs(trends$slope[trends$metric == "loss rate"])))+
  ylim(-max(abs(trends$slope[trends$metric == "recruitment rate"])), max(abs(trends$slope[trends$metric == "recruitment rate"])))+
  scale_y_continuous(trans= ssqrt_trans, labels = scales::scientific)+
  scale_x_continuous(trans= ssqrt_trans, labels = scales::scientific)+
  scale_color_manual(values = unique(sort(tmp$color)))+
  theme_bw()

dev.off()

