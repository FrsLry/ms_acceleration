## Figures
library(sf)
library(mgcv)
library(dplyr)
library(rnaturalearth)
library(ggplot2)
library(tidyr)
library(tibble)
library(ggforce)
library(ggrepel)
## Load the models' output
for(file in list.files("mixed_model/outputs/", full.names = T)){
  assign(gsub("^mixed_model/outputs/|\\.rds$","",file),
         readRDS(file = file))
}

## Load shapefile
routes_shp <- st_read("data/routes_shp/routes_selected_lines.shp") %>% arrange(routes)

## Add values of beta1 (slope) to the shapefile
routes_shp$ab_trend <- overall_N_output$mean$beta1
routes_shp$growth_rate <- overall_g_output$mean$beta1
routes_shp$rec_rate <- overall_r_output$mean$beta1
routes_shp$loss_rate <- overall_l_output$mean$beta1
routes_shp$gt0 <- overall_gt0_output$mean$beta1
routes_shp$g_lower <- overall_g_output$q2.5$beta1
routes_shp$g_upper <- overall_g_output$q97.5$beta1

## Add lat and long of route centroid
routes_shp <- routes_shp %>%
  mutate(lon = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
         lat = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(Y))

## Perform GAMs
ab_gam <- gam(ab_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
              data = routes_shp)
growth_rate_gam <- gam(growth_rate ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                       data = routes_shp)
rec_rate_gam <- gam(rec_rate ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                    data = routes_shp)
loss_rate_gam <- gam(loss_rate ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                     data = routes_shp)
gt0_gam <- gam(gt0 ~ s(lon, lat, bs = "gp", k = 100, m = 2),
               data = routes_shp)

## Add GAMs values to shapefile
routes_shp$ab_trend_gam <- ab_gam$fitted.values
routes_shp$growth_rate_gam <- growth_rate_gam$fitted.values
routes_shp$rec_rate_gam <- rec_rate_gam$fitted.values
routes_shp$loss_rate_gam <- loss_rate_gam$fitted.values
routes_shp$gt0_gam <- gt0_gam$fitted.values

## get the states polygon
states_shp <-
  ne_states(returnclass = "sf") %>%
  filter(adm0_a3 == "CAN" | adm0_a3 == "USA" | adm0_a3 == "MEX") %>%
  filter(woe_name != "Hawaii") %>%
  st_crop(st_bbox(c(xmin = -165,
                    xmax = -60,
                    ymax = 62,
                    ymin = 16)))

##### Spatial plots ###########
## Make gam maps
pdf("mixed_model/figures/trend_maps_gam.pdf",
    width=11, height=8.5)
for(metric in c("ab_trend_gam",
                "growth_rate_gam",
                "rec_rate_gam",
                "loss_rate_gam",
                "gt0_gam")){
  print(
    routes_shp %>%
      st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
      st_centroid() %>%
      rename(variable = metric) %>%
      ggplot()+
      geom_sf(data = states_shp)+
      geom_sf(aes(color = variable), size = 5, show.legend = FALSE)+
      scale_color_gradient2(low = "#e31a1cff", midpoint = 0, mid = "#ffff99", high = "#1f78b4ff")+
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
    scale_color_gradient2(low = "#e31a1cff", midpoint = 0, mid = "#ffff99", high = "#1f78b4ff",
                          labels = scales::scientific_format())+
    ggtitle(metric)+
    labs(colour = metric)+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18))

  grid::grid.newpage()
  grid::grid.draw(cowplot::get_legend(plot))

}

dev.off()

pdf("mixed_model/figures/trend_maps.pdf",
    width=11, height=8.5)
## Make non gam maps.
# We use the color palette of the GAM in order to show that the patterns infered from the GAM make sense
for(metric in c("ab_trend",
                "growth_rate",
                "rec_rate",
                "loss_rate",
                "gt0")){
    print(
      routes_shp %>%
        st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
        st_centroid() %>%
        rename(variable = metric) %>%
        ggplot()+
        geom_sf(data = states_shp)+
        geom_sf(aes(color = variable), size = 5, show.legend = F)+
        # scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
        scale_color_gradientn(colors = c("#99000D", "#e31a1cff", "#ffffbf", "#1f78b4ff", "#0C4B8E"),
                              values = scales::rescale(c(min(st_drop_geometry(routes_shp[,metric])),
                                                         min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")])),
                                                         0,
                                                         ifelse(metric %in% c("ab_trend", "gt0"), abs(min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))), max(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))),
                                                         max(st_drop_geometry(routes_shp[,metric])))))+
        geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
        ggtitle(metric)+
        theme_bw()+
        # labs(colour = paste0("log(",metric,")"))+
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
    # scale_color_gradient2(low = "#a50026", midpoint = 0, mid = "#ffffbf", high = "#313695")+
    scale_color_gradientn(colors = c("#99000D", "#e31a1cff", "#ffffbf", "#1f78b4ff", "#0C4B8E"),
                          values = scales::rescale(c(min(st_drop_geometry(routes_shp[,metric])),
                                                     min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")])),
                                                     0,
                                                     ifelse(metric %in% c("ab_trend", "gt0"), abs(min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))), max(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))),
                                                     max(st_drop_geometry(routes_shp[,metric])))))+
    geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
    ggtitle(metric)+
    theme_bw()+
    labs(colour = metric)+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18))

  grid::grid.newpage()
  grid::grid.draw(cowplot::get_legend(plot))

}

dev.off()

## Histograms of the posteriors of the grand.beta1 of the mixed models
pdf("mixed_model/figures/hist_posterior.pdf",
    width = 5.83, height = 4.13)

for(metric in c("overall_N_output", "overall_g_output",
                "overall_gt0_output", "overall_l_output",
                "overall_r_output")){

  mcmc_samples <- get(metric)$samples[,c("grand.beta1")]
  mcmc_df <- do.call(rbind, lapply(mcmc_samples, as.data.frame))

  print(
    ggplot(mcmc_df, aes(x = var1)) +
      geom_density(alpha = 0.5, fill = "black")+
      geom_vline(aes(xintercept = get(metric)$mean$grand.beta1), color="#e31a1cff") +
      geom_vline(aes(xintercept = 0), linetype = "dashed") +
      geom_vline(aes(xintercept = get(metric)$q2.5$grand.beta1), color="#1f78b4ff", linetype = "dashed", size = .5)+
      geom_vline(aes(xintercept = get(metric)$q97.5$grand.beta1), color="#1f78b4ff", linetype = "dashed", size = .5)+
      ylab(metric) +
      theme_bw()
  )

}

dev.off()

## Trends plot
pdf("mixed_model/figures/trends_numbers.pdf",
    width = 5.83, height = 4.13)
for(metric in c("N", "gt0", "g", "rr", "lr")){

  d <- readRDS(paste0("mixed_model/save_samples/mean_",metric,"_overall.rds"))

  d <- d %>% as.data.frame()
  colnames(d) <- c(1:ncol(d))

  print(d %>% rownames_to_column(var = "route") %>%
          pivot_longer(cols = -route, names_to = "year", values_to = "value") %>%
          mutate(year = as.numeric(year)) %>%
          ggplot()+
          geom_line(aes(year, value, group = route), color = "grey")+
          geom_hline(yintercept = ifelse(metric %in% c("gt0", "g"), 0, NA), linetype="dashed", color = "black")+
          stat_summary(aes(year, value), fun = median, geom = "line", color = "blue",
                       linewidth = 1, linetype = 2)+
          geom_hline(yintercept = ifelse(metric %in% c("N"), median(d[,1]), NA), linetype="dashed", color = "black")+
          scale_y_continuous(trans= ggallin::ssqrt_trans)+
          ylab(metric)+
          theme_light())

}

dev.off()

## Rec rate vs. Loss rate
## Add the color palette
routes_shp <- routes_shp %>% mutate(recRate_vs_lossRate = case_when(
  (rec_rate > 0 & loss_rate > 0 & abs(rec_rate) > abs(loss_rate)) ~ "#B2DF8A",
  (rec_rate > 0 & loss_rate > 0 & abs(rec_rate) < abs(loss_rate)) ~ "#FB9A99",
  (rec_rate < 0 & loss_rate > 0 & abs(rec_rate) < abs(loss_rate)) ~ "#E31A1C",
  (rec_rate < 0 & loss_rate > 0 & abs(rec_rate) > abs(loss_rate)) ~ "#FF7F00",
  (rec_rate < 0 & loss_rate < 0 & abs(rec_rate) > abs(loss_rate)) ~ "#FDBF6F",
  (rec_rate < 0 & loss_rate < 0 & abs(rec_rate) < abs(loss_rate)) ~ "#A6CEE3",
  (rec_rate > 0 & loss_rate < 0 & abs(rec_rate) < abs(loss_rate)) ~ "#1F78B4",
  (rec_rate > 0 & loss_rate < 0 & abs(rec_rate) > abs(loss_rate)) ~ "#33A02C"),
  recRate_vs_lossRate_gam = case_when(
    (rec_rate_gam > 0 & loss_rate_gam > 0 & abs(rec_rate_gam) > abs(loss_rate_gam)) ~ "#B2DF8A",
    (rec_rate_gam > 0 & loss_rate_gam > 0 & abs(rec_rate_gam) < abs(loss_rate_gam)) ~ "#FB9A99",
    (rec_rate_gam < 0 & loss_rate_gam > 0 & abs(rec_rate_gam) < abs(loss_rate_gam)) ~ "#E31A1C",
    (rec_rate_gam < 0 & loss_rate_gam > 0 & abs(rec_rate_gam) > abs(loss_rate_gam)) ~ "#FF7F00",
    (rec_rate_gam < 0 & loss_rate_gam < 0 & abs(rec_rate_gam) > abs(loss_rate_gam)) ~ "#FDBF6F",
    (rec_rate_gam < 0 & loss_rate_gam < 0 & abs(rec_rate_gam) < abs(loss_rate_gam)) ~ "#A6CEE3",
    (rec_rate_gam > 0 & loss_rate_gam < 0 & abs(rec_rate_gam) < abs(loss_rate_gam)) ~ "#1F78B4",
    (rec_rate_gam > 0 & loss_rate_gam < 0 & abs(rec_rate_gam) > abs(loss_rate_gam)) ~ "#33A02C"))

## Fig rec rate vs loss rate
svg("mixed_model/figures/r_vs_l.svg",
    height = 4.13, width = 5.83)

df <-
  routes_shp %>% select(ab_trend, rec_rate, loss_rate, recRate_vs_lossRate) %>% st_drop_geometry() %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase")) %>%
  ## Add the sd
  add_column(upper_r = overall_r_output$q97.5$beta1,
             lower_r = overall_r_output$q2.5$beta1,
             upper_l = overall_l_output$q97.5$beta1,
             lower_l = overall_l_output$q2.5$beta1) %>%
  arrange(recRate_vs_lossRate)

ggplot(df)+
  geom_errorbar(aes(x = loss_rate, y = rec_rate,
                    ymin = lower_r, ymax = upper_r,
                    colour = recRate_vs_lossRate),
                alpha = .3, size = .2,
                show.legend = F)+
  geom_errorbar(aes(x = loss_rate, y = rec_rate,
                    xmin = lower_l, xmax = upper_l,
                    colour = recRate_vs_lossRate),
                alpha = .3, size = .2,
                show.legend = F)+
  scale_colour_manual(values = unique(sort(df$recRate_vs_lossRate)))+
  geom_point(aes(x = loss_rate, y = rec_rate, color = recRate_vs_lossRate),
             shape = ifelse(df$nSign == "increase", "\u2191", "\u2193"),
             size = 4, show.legend = F)+
  scale_color_manual(values = unique(sort(df$recRate_vs_lossRate)))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)+
  geom_hline(yintercept = 0, linetype="dashed", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", size = 1)+
  geom_abline(slope = 1, linetype = "longdash")+
  geom_abline(slope = -1, linetype = "longdash")+
  theme_bw()

dev.off()

## Fig rec rate vs loss rate gam
svg("mixed_model/figures/r_vs_l_gam.svg",
    height = 4.13, width = 5.83)

df <-
  routes_shp %>% select(ab_trend_gam, rec_rate_gam, loss_rate_gam, recRate_vs_lossRate_gam) %>% st_drop_geometry() %>%
  arrange(recRate_vs_lossRate_gam) %>%
  mutate(nSign = ifelse(ab_trend_gam < 0, "decrease", "increase"))

ggplot(df)+
  geom_point(aes(loss_rate_gam, rec_rate_gam, color = recRate_vs_lossRate_gam),
             shape = ifelse(df$nSign == "increase", "\u2191", "\u2193"),
             size = 4, show.legend = F)+
  scale_color_manual(values = unique(sort(df$recRate_vs_lossRate_gam)))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)+
  geom_hline(yintercept = 0, linetype="dashed", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", size = 1)+
  geom_abline(slope = 1, linetype = "longdash")+
  geom_abline(slope = -1, linetype = "longdash")+
  theme_bw()

dev.off()

pdf("mixed_model/figures/maps_rec_vs_loss.pdf",
    width=11, height=8.5)
## Rec vs loss maps
routes_shp %>%
  st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
  st_centroid() %>%
  ggplot()+
  geom_sf(data = states_shp)+
  geom_sf(aes(color = recRate_vs_lossRate), size = 5, show.legend = F)+
  scale_color_manual(values = unique(sort(routes_shp$recRate_vs_lossRate)))+
  geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18))

routes_shp %>%
  st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
  st_centroid() %>%
  ggplot()+
  geom_sf(data = states_shp)+
  geom_sf(aes(color = recRate_vs_lossRate_gam), size = 5, show.legend = F)+
  scale_color_manual(values = unique(sort(routes_shp$recRate_vs_lossRate_gam)))+
  geom_sf(fill = NA, data = states_shp, color = alpha("grey",0.5))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18))

dev.off()

## Boxplot per route not gam
pdf("mixed_model/figures/boxplot_route.pdf",
    height = 2, width = 5.83)
routes_shp %>%
  mutate(order = case_when(
    recRate_vs_lossRate == "#FF7F00" ~ 1,
    recRate_vs_lossRate == "#FDBF6F" ~ 2,
    recRate_vs_lossRate == "#E31A1C" ~ 3,
    recRate_vs_lossRate == "#FB9A99" ~ 4,
    recRate_vs_lossRate == "#33A02C" ~ 5,
    recRate_vs_lossRate == "#B2DF8A" ~ 6,
    recRate_vs_lossRate == "#1F78B4" ~ 7,
    recRate_vs_lossRate == "#A6CEE3" ~ 8
  )) %>%
  ggplot()+
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  geom_boxplot(aes(reorder(recRate_vs_lossRate, order), ab_trend, colour = recRate_vs_lossRate),
               show.legend = F, width = 0.5, outlier.shape = NA)+
  # geom_violin(aes(reorder(rec_vs_loss, order), ab_trend, fill = rec_vs_loss), show.legend = F, width = 1)+
  ggforce::geom_sina(aes(reorder(recRate_vs_lossRate, order), ab_trend), show.legend = F, alpha = .5, size = .5)+
  scale_fill_manual(values = unique(sort(routes_shp$recRate_vs_lossRate)))+
  scale_color_manual(values = unique(sort(routes_shp$recRate_vs_lossRate)))+
  scale_y_continuous(trans= ggallin::pseudolog10_trans)+
                     # breaks = c(-4000,-2000,-10,0,10,2000))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title  = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11))

dev.off()

## Boxplot per route gam
pdf("mixed_model/figures/boxplot_route_gam.pdf",
    height = 2, width = 5.83)
routes_shp %>%
  mutate(order = case_when(
    recRate_vs_lossRate_gam == "#FF7F00" ~ 1,
    recRate_vs_lossRate_gam == "#FDBF6F" ~ 2,
    recRate_vs_lossRate_gam == "#E31A1C" ~ 3,
    recRate_vs_lossRate_gam == "#FB9A99" ~ 4,
    recRate_vs_lossRate_gam == "#33A02C" ~ 5,
    recRate_vs_lossRate_gam == "#B2DF8A" ~ 6,
    recRate_vs_lossRate_gam == "#1F78B4" ~ 7,
    recRate_vs_lossRate_gam == "#A6CEE3" ~ 8
  )) %>%
  ggplot()+
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  geom_boxplot(aes(reorder(recRate_vs_lossRate_gam, order), ab_trend_gam, colour = recRate_vs_lossRate_gam),
               show.legend = F, width = 0.5, outlier.shape = NA)+
  # geom_violin(aes(reorder(rec_vs_loss, order), ab_trend, fill = rec_vs_loss), show.legend = F, width = 1)+
  ggforce::geom_sina(aes(reorder(recRate_vs_lossRate_gam, order), ab_trend_gam), show.legend = F, alpha = .5, size = .5)+
  scale_fill_manual(values = unique(sort(routes_shp$recRate_vs_lossRate_gam)))+
  scale_color_manual(values = unique(sort(routes_shp$recRate_vs_lossRate_gam)))+
  scale_y_continuous(trans= ggallin::pseudolog10_trans)+
                     # breaks = c(-4000,-2000,-10,0,10,2000))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title  = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11))

dev.off()


#### Per species plots ########
## Load the models
perspecies_N_output <- readRDS("mixed_model/outputs/perspecies_N_output.rds")
perspecies_r_output <- readRDS("mixed_model/outputs/perspecies_r_output.rds")
perspecies_l_output <- readRDS("mixed_model/outputs/perspecies_l_output.rds")

## Get back the species names based of median Rhat < 1.1 #####
file.list <-list.files(path = "summary_models/", full.names = T)
species <- gsub("^summary_|\\.rds$","",list.files(path = "summary_models/"))
for(i in 1:length(file.list)){
  assign(species[i],
         readRDS(file.list[i]))
}
## Let's check the Rhats
rhats <- data.frame(species = as.character(), param = as.character(), Rhat = as.numeric())
for(sp in species){
  d <- get(sp)
  rhats <- rbind(rhats,
                 c(sp, "avg_rhat", median(d$Rhat, na.rm = T)))
}
colnames(rhats) <- c("species", "param", "rhat")
rhats$rhat <- as.numeric(rhats$rhat)

## Filter the species with median Rhat <= 1.1
species <-
  rhats %>% filter(!grepl("mean.", .$param)) %>%
  group_by(species) %>%
  filter(all(rhat <= 1.1)) %>%
  distinct(species) %>% pull()

## Remove the per species models
rm(list = ls()[ls() %in% gsub("^summary_|\\.rds$","",list.files(path = "summary_models/"))])
################

species_data <- data.frame(species = as.character(),
                           ab_trend = as.numeric(),
                           r = as.numeric(),
                           r_upper = as.numeric(),
                           r_lower = as.numeric(),
                           l = as.numeric(),
                           l_upper = as.numeric(),
                           l_lower = as.numeric(),
                           g = as.numeric(),
                           g_upper = as.numeric(),
                           g_lower = as.numeric())

for(i in 1:length(perspecies_N_output$mean$beta1)){

  s = species[i]
  a = perspecies_N_output$mean$beta1[i]
  r = perspecies_r_output$mean$beta1[i]
  r_upper = perspecies_r_output$q97.5$beta1[i]
  r_lower = perspecies_r_output$q2.5$beta1[i]
  l = perspecies_l_output$mean$beta1[i]
  l_upper = perspecies_l_output$q97.5$beta1[i]
  l_lower = perspecies_l_output$q2.5$beta1[i]
  g = perspecies_g_output$mean$beta1[i]
  g_upper = perspecies_g_output$q97.5$beta1[i]
  g_lower = perspecies_g_output$q2.5$beta1[i]

  species_data <- rbind(species_data,
                        data.frame(species = s,
                                   ab_trend = a,
                                   r = r,
                                   r_upper = r_upper,
                                   r_lower = r_lower,
                                   l = l,
                                   l_upper = l_upper,
                                   l_lower = l_lower,
                                   g = g,
                                   g_upper = g_upper,
                                   g_lower = g_lower))

}

## Add the color palette
tmp <-
species_data %>%
  mutate(rec_vs_loss = case_when(
    (r > 0 & l > 0 & abs(r) > abs(l)) ~ "#B2DF8A",
    (r > 0 & l > 0 & abs(r) < abs(l)) ~ "#FB9A99",
    (r < 0 & l > 0 & abs(r) < abs(l)) ~ "#E31A1C",
    (r < 0 & l > 0 & abs(r) > abs(l)) ~ "#FF7F00",
    (r < 0 & l < 0 & abs(r) > abs(l)) ~ "#FDBF6F",
    (r < 0 & l < 0 & abs(r) < abs(l)) ~ "#A6CEE3",
    (r > 0 & l < 0 & abs(r) < abs(l)) ~ "#1F78B4",
    (r > 0 & l < 0 & abs(r) > abs(l)) ~ "#33A02C"
  )) %>%
  arrange(rec_vs_loss) %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase")) %>%
  ## remove one species that increase the scale of the plot
  filter(species != "Eurasian Collared-Dove")

## Rec vs loss per species
svg("mixed_model/figures/rec_vs_loss_species.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = l, y = r,
                    ymin = r_lower, ymax = r_upper,
                    colour = rec_vs_loss),
                alpha = .3, size = .2,
                show.legend = F)+
  geom_errorbar(aes(x = l, y = r,
                    xmin = l_lower, xmax = l_upper,
                    colour = rec_vs_loss),
                alpha = .3, size = .2,
                show.legend = F)+
  scale_colour_manual(values = unique(sort(tmp$rec_vs_loss)))+
  geom_point(aes(l, r, color = rec_vs_loss),
             shape = ifelse(tmp$nSign == "increase", "\u2191", "\u2193"),
             show.legend = F, size = 4)+
  geom_abline(slope = 1, linetype = "longdash", color = "#595959ff")+
  geom_abline(slope = -1, linetype = "longdash", color = "#595959ff")+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(l, r, label = species), check_overlap = T,
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10,  ## too control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_x_continuous(labels = scales::scientific)+
  scale_y_continuous(labels = scales::scientific)+
  theme_bw()

dev.off()

## Boxplot per species
pdf("mixed_model/figures/boxplot_perspecies.pdf",
    height = 2, width = 5.83)

tmp <-
  species_data %>%
  mutate(rec_vs_loss = case_when(
    (r > 0 & l > 0 & abs(r) > abs(l)) ~ "#B2DF8A",
    (r > 0 & l > 0 & abs(r) < abs(l)) ~ "#FB9A99",
    (r < 0 & l > 0 & abs(r) < abs(l)) ~ "#E31A1C",
    (r < 0 & l > 0 & abs(r) > abs(l)) ~ "#FF7F00",
    (r < 0 & l < 0 & abs(r) > abs(l)) ~ "#FDBF6F",
    (r < 0 & l < 0 & abs(r) < abs(l)) ~ "#A6CEE3",
    (r > 0 & l < 0 & abs(r) < abs(l)) ~ "#1F78B4",
    (r > 0 & l < 0 & abs(r) > abs(l)) ~ "#33A02C"
  )) %>%
  arrange(rec_vs_loss) %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))

tmp %>%
  mutate(order = case_when(
    rec_vs_loss == "#FF7F00" ~ 1,
    rec_vs_loss == "#FDBF6F" ~ 2,
    rec_vs_loss == "#E31A1C" ~ 3,
    rec_vs_loss == "#FB9A99" ~ 4,
    rec_vs_loss == "#33A02C" ~ 5,
    rec_vs_loss == "#B2DF8A" ~ 6,
    rec_vs_loss == "#1F78B4" ~ 7,
    rec_vs_loss == "#A6CEE3" ~ 8
  )) %>%
  ggplot()+
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  geom_boxplot(aes(reorder(rec_vs_loss, order), ab_trend, colour = rec_vs_loss),
               show.legend = F, width = 0.5, outlier.shape = NA)+
  # geom_violin(aes(reorder(rec_vs_loss, order), ab_trend, fill = rec_vs_loss), show.legend = F, width = 1)+
  ggforce::geom_sina(aes(reorder(rec_vs_loss, order), ab_trend), show.legend = F, alpha = .5, size = .5)+
  scale_fill_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_y_continuous(trans= ggallin::pseudolog10_trans,
                     breaks = c(-4000,-2000,-10,0,10,2000))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title  = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11))

dev.off()

## Abundance and growth rate trends per species
pdf("mixed_model/figures/trends_species.pdf",
    width = 5.83, height = 4.13)
for(metric in c("N", "g")){

  d <- readRDS(paste0("mixed_model/save_samples/mean_sd_",metric,"_sp.rds"))

  d <- d %>% select(-paste0("sd_",metric)) %>% pivot_wider(names_from = year, values_from = paste0("mean_",metric))

  print(d %>%
          as.data.frame() %>%
          pivot_longer(cols = -species, names_to = "year", values_to = "value") %>%
          mutate(year = as.numeric(year)) %>%
          ggplot()+
          geom_line(aes(year, value, group = species), color = "grey")+
          geom_hline(yintercept = ifelse(metric %in% c("gt0", "g"), 0, NA), linetype="dashed", color = "black")+
          stat_summary(aes(year, value), fun = median, geom = "line", color = "blue",
                       linewidth = 1, linetype = 2)+
          geom_hline(yintercept = ifelse(metric %in% c("N"), median(pull(d[,2])), NA), linetype="dashed", color = "black")+
          scale_y_continuous(trans= ggallin::ssqrt_trans)+
          ylab(metric)+
          theme_light())

}
dev.off()

## Per family plots #####
## Load the models
perfamily_N_output <- readRDS("mixed_model/outputs/perfamily_N_output.rds")
perfamily_r_output <- readRDS("mixed_model/outputs/perfamily_r_output.rds")
perfamily_l_output <- readRDS("mixed_model/outputs/perfamily_l_output.rds")

## Get back the family names
family <- readr::read_csv("data/phylo.csv")

## Filter the species with good Rhat
family = family[family$English_Common_Name %in% species,]

## Extract family names
families <- unique(family$Family)

families_data <- data.frame(family = as.character(),
                            ab_trend = as.numeric(),
                            r = as.numeric(),
                            r_upper = as.numeric(),
                            r_lower = as.numeric(),
                            l = as.numeric(),
                            l_upper = as.numeric(),
                            l_lower = as.numeric())

for(i in 1:length(perfamily_N_output$mean$beta1)){

  f = families[i]
  a = perfamily_N_output$mean$beta1[i]
  r = perfamily_r_output$mean$beta1[i]
  r_upper = perfamily_r_output$q97.5$beta1[i]
  r_lower = perfamily_r_output$q2.5$beta1[i]
  l = perfamily_l_output$mean$beta1[i]
  l_upper = perfamily_l_output$q97.5$beta1[i]
  l_lower = perfamily_l_output$q2.5$beta1[i]

  families_data <- rbind(families_data,
                        data.frame(family = f,
                                   ab_trend = a,
                                   r = r,
                                   r_upper = r_upper,
                                   r_lower = r_lower,
                                   l = l,
                                   l_upper = l_upper,
                                   l_lower = l_lower))

}

## Add the color palette
tmp <-
  families_data %>%
  mutate(rec_vs_loss = case_when(
    (r > 0 & l > 0 & abs(r) > abs(l)) ~ "#B2DF8A",
    (r > 0 & l > 0 & abs(r) < abs(l)) ~ "#FB9A99",
    (r < 0 & l > 0 & abs(r) < abs(l)) ~ "#E31A1C",
    (r < 0 & l > 0 & abs(r) > abs(l)) ~ "#FF7F00",
    (r < 0 & l < 0 & abs(r) > abs(l)) ~ "#FDBF6F",
    (r < 0 & l < 0 & abs(r) < abs(l)) ~ "#A6CEE3",
    (r > 0 & l < 0 & abs(r) < abs(l)) ~ "#1F78B4",
    (r > 0 & l < 0 & abs(r) > abs(l)) ~ "#33A02C"
  )) %>%
  arrange(rec_vs_loss) %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))

## Rec vs loss per family
svg("mixed_model/figures/rec_vs_loss_families.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = l, y = r,
                    ymin = r_lower, ymax = r_upper,
                    colour = rec_vs_loss),
                alpha = .3, size = .2,
                show.legend = F)+
  geom_errorbar(aes(x = l, y = r,
                    xmin = l_lower, xmax = l_upper,
                    colour = rec_vs_loss),
                alpha = .3, size = .2,
                show.legend = F)+
  scale_colour_manual(values = unique(sort(tmp$rec_vs_loss)))+
  geom_point(aes(l, r, color = rec_vs_loss),
             shape = ifelse(tmp$nSign == "increase", "\u2191", "\u2193"),
             show.legend = F, size = 4)+
  geom_abline(slope = 1, linetype = "longdash", color = "#595959ff")+
  geom_abline(slope = -1, linetype = "longdash", color = "#595959ff")+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(l, r, label = family), check_overlap = T,
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 4,  ## too control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_y_continuous(labels = scales::scientific)+
  scale_x_continuous(labels = scales::scientific)+
  theme_bw()

dev.off()

## Boxplot per family
pdf("mixed_model/figures/boxplot_perfamily.pdf",
    height = 2, width = 5.83)

tmp %>%
  mutate(order = case_when(
    rec_vs_loss == "#FF7F00" ~ 1,
    rec_vs_loss == "#FDBF6F" ~ 2,
    rec_vs_loss == "#E31A1C" ~ 3,
    rec_vs_loss == "#FB9A99" ~ 4,
    rec_vs_loss == "#33A02C" ~ 5,
    rec_vs_loss == "#B2DF8A" ~ 6,
    rec_vs_loss == "#1F78B4" ~ 7,
    rec_vs_loss == "#A6CEE3" ~ 8
  )) %>%
  mutate(rec_vs_loss = factor(rec_vs_loss, levels = c(
    "#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99",
    "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3"
  ))) %>%
  ggplot()+
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  geom_boxplot(aes(reorder(rec_vs_loss, order), ab_trend, colour = rec_vs_loss),
               show.legend = F, width = 0.5, outlier.shape = NA)+
  # geom_violin(aes(reorder(rec_vs_loss, order), ab_trend, fill = rec_vs_loss), show.legend = F, width = 1)+
  ggforce::geom_sina(aes(reorder(rec_vs_loss, order), ab_trend), show.legend = F, alpha = .5, size = .5)+
  # scale_fill_manual(values = c(
  #   "#e66101fa", "#feb963fa", "#ca0020fa", "#f5a582fa",
  #   "#33a02cfa", "#b3df8afa", "#1f78b5fa", "#a6cee3fa"
  # ))+
  scale_color_manual(values = c(
    "#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99",
    "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3"
  ))+
  scale_y_continuous(trans= ggallin::pseudolog10_trans,
                     breaks = c(-5000,-2000,-10,0,10,1000))+
  scale_x_discrete(drop = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title  = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11))

dev.off()


## Per habitats plots #####
## Load the models
perhabitat_N_output <- readRDS("mixed_model/outputs/perhabitat_N_output.rds")
perhabitat_r_output <- readRDS("mixed_model/outputs/perhabitat_r_output.rds")
perhabitat_l_output <- readRDS("mixed_model/outputs/perhabitat_l_output.rds")

## Get back the family names
mean_N_habitats <- readRDS("mixed_model/save_samples/mean_N_habitats.rds")

## Extrat family names
habitats <- rownames(mean_N_habitats)

habitats_data <- data.frame(habitat = as.character(),
                            ab_trend = as.numeric(),
                            r = as.numeric(),
                            r_upper = as.numeric(),
                            r_lower = as.numeric(),
                            l = as.numeric(),
                            l_upper = as.numeric(),
                            l_lower = as.numeric())

for(i in 1:length(perhabitat_N_output$mean$beta1)){

  h = habitats[i]
  a = perhabitat_N_output$mean$beta1[i]
  r = perhabitat_r_output$mean$beta1[i]
  r_upper = perhabitat_r_output$q97.5$beta1[i]
  r_lower = perhabitat_r_output$q2.5$beta1[i]
  l = perhabitat_l_output$mean$beta1[i]
  l_upper = perhabitat_l_output$q97.5$beta1[i]
  l_lower = perhabitat_l_output$q2.5$beta1[i]

  habitats_data <- rbind(habitats_data,
                         data.frame(habitat = h,
                                    ab_trend = a,
                                    r = r,
                                    r_upper = r_upper,
                                    r_lower = r_lower,
                                    l = l,
                                    l_upper = l_upper,
                                    l_lower = l_lower))

}

## Add the color palette
tmp <-
  habitats_data %>%
  mutate(rec_vs_loss = case_when(
    (r > 0 & l > 0 & abs(r) > abs(l)) ~ "#B2DF8A",
    (r > 0 & l > 0 & abs(r) < abs(l)) ~ "#FB9A99",
    (r < 0 & l > 0 & abs(r) < abs(l)) ~ "#E31A1C",
    (r < 0 & l > 0 & abs(r) > abs(l)) ~ "#FF7F00",
    (r < 0 & l < 0 & abs(r) > abs(l)) ~ "#FDBF6F",
    (r < 0 & l < 0 & abs(r) < abs(l)) ~ "#A6CEE3",
    (r > 0 & l < 0 & abs(r) < abs(l)) ~ "#1F78B4",
    (r > 0 & l < 0 & abs(r) > abs(l)) ~ "#33A02C"
  )) %>%
  arrange(rec_vs_loss) %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))

## Rec vs loss per habitat
svg("mixed_model/figures/rec_vs_loss_habitats.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = l, y = r,
                    ymin = r_lower, ymax = r_upper,
                    colour = rec_vs_loss),
                alpha = .3, size = .2,
                show.legend = F)+
  geom_errorbar(aes(x = l, y = r,
                    xmin = l_lower, xmax = l_upper,
                    colour = rec_vs_loss),
                alpha = .3, size = .2,
                show.legend = F)+
  scale_colour_manual(values = unique(sort(tmp$rec_vs_loss)))+
  geom_point(aes(l, r, color = rec_vs_loss),
             shape = ifelse(tmp$nSign == "increase", "\u2191", "\u2193"),
             show.legend = F, size = 4)+
  geom_abline(slope = 1, linetype = "longdash", color = "#595959ff")+
  geom_abline(slope = -1, linetype = "longdash", color = "#595959ff")+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(l, r, label = habitat), check_overlap = T,
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  # max.overlaps = 4,  ## too control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_y_continuous(labels = scales::scientific)+
  scale_x_continuous(labels = scales::scientific)+
  theme_bw()

dev.off()

## Trends per habitat
pdf("mixed_model/figures/trends_habitat.pdf",
    width = 7.83, height = 4.13)
for(metric in c("N", "g", "rr", "lr")){

  d <- readRDS(paste0("mixed_model/save_samples/mean_",metric,"_habitats.rds"))

  d <- d %>% as.data.frame()
  colnames(d) <- c(1:ncol(d))

  print(d %>% rownames_to_column(var = "habitat") %>%
          pivot_longer(cols = -habitat, names_to = "year", values_to = "value") %>%
          mutate(year = as.numeric(year)) %>%
          ggplot()+
          geom_line(aes(year, value, group = habitat, color = habitat))+
          geom_hline(yintercept = ifelse(metric %in% c("gt0", "g"), 0, NA), linetype="dashed", color = "black")+
          # stat_summary(aes(year, value), fun = median, geom = "line", color = "blue",
          #              linewidth = 1, linetype = 2)+
          # geom_hline(yintercept = ifelse(metric %in% c("N"), median(d[,1]), NA), linetype="dashed", color = "black")+
          scale_y_continuous(trans= ggallin::ssqrt_trans)+
          ylab(metric)+
          theme_light())

}

dev.off()

