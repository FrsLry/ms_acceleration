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
routes_shp$g_lower <- overall_g_output$q2.5$beta1
routes_shp$g_upper <- overall_g_output$q97.5$beta1
routes_shp$absg <- overall_absg_output$mean$beta1

## Add lat and long of route centroid
routes_shp <- routes_shp %>%
  mutate(lon = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
         lat = st_centroid(routes_shp) %>% st_coordinates() %>% as.data.frame() %>% pull(Y))

## Perform GAMs
ab_gam <- gam(ab_trend ~ s(lon, lat, bs = "gp", k = 100, m = 2),
              data = routes_shp)
growth_rate_gam <- gam(growth_rate ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                       data = routes_shp)
absg_gam <- gam(absg ~ s(lon, lat, bs = "gp", k = 100, m = 2),
                data = routes_shp)

## Add GAMs values to shapefile
routes_shp$ab_trend_gam <- ab_gam$fitted.values
routes_shp$growth_rate_gam <- growth_rate_gam$fitted.values
routes_shp$absg_gam <- absg_gam$fitted.values

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
                "absg_gam")){
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
                          labels = scales::comma_format()
                          )+
    ggtitle(metric)+
    labs(colour = metric)+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18))

  grid::grid.newpage()
  grid::grid.draw(cowplot::get_legend(plot))

}

dev.off()

## Make non gam maps.
# We use the color palette of the GAM in order to show that the patterns inferred from the GAM make sense
pdf("mixed_model/figures/trend_maps.pdf",
    width=11, height=8.5)

for(metric in c("ab_trend",
                "growth_rate",
                "absg")){
    print(
      routes_shp %>%
        st_transform(crs = "+proj=aea +lon_0=-108.4570304 +lat_1=34.300327 +lat_2=64.5191945 +lat_0=49.4097608 +datum=WGS84 +units=m +no_defs") %>%
        st_centroid() %>%
        rename(variable = metric) %>%
        ggplot()+
        geom_sf(data = states_shp)+
        geom_sf(aes(color = variable), size = 5, show.legend = F)+
        scale_color_gradientn(colors = c("#99000D", "#e31a1cff", "#ffff99", "#1f78b4ff", "#0C4B8E"),
                              values = scales::rescale(c(min(st_drop_geometry(routes_shp[,metric])),
                                                         min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")])),
                                                         0,
                                                         ifelse(metric %in% c("ab_trend"), abs(min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))), max(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))),
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
    scale_color_gradientn(colors = c("#99000D", "#e31a1cff", "#ffff99", "#1f78b4ff", "#0C4B8E"),
                          values = scales::rescale(c(min(st_drop_geometry(routes_shp[,metric])),
                                                     min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")])),
                                                     0,
                                                     ifelse(metric %in% c("ab_trend"), abs(min(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))), max(st_drop_geometry(routes_shp[,paste0(metric, "_gam")]))),
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

for(metric in c("overall_N_output", "overall_g_output", "overall_absg_output")){

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

for(metric in c("N", "g", "absg")){

  d <- readRDS(paste0("mixed_model/save_samples/mean_",metric,"_overall.rds"))

  d <- d %>% as.data.frame()
  colnames(d) <- c(1:ncol(d))

  print(d %>% rownames_to_column(var = "route") %>%
          pivot_longer(cols = -route, names_to = "year", values_to = "value") %>%
          mutate(year = as.numeric(year)) %>%
          ggplot()+
          geom_line(aes(year, value, group = route), color = "grey")+
          geom_hline(yintercept = ifelse(metric %in% c("g", "absg"), 0, NA_real_), linetype="dashed", color = "black")+
          stat_summary(aes(year, value), fun = median, geom = "line", color = "blue",
                       linewidth = 1, linetype = 2)+
          scale_y_continuous(trans= if (metric %in% c("g","absg")) ggallin::ssqrt_trans else "identity")+
          ylab(metric)+
          theme_light())

}

dev.off()

#### Per species plots ########
## Load the models
perspecies_N_output <- readRDS("mixed_model/outputs/perspecies_N_output.rds")

## Get back the species names for which all parameters Rhat < 1.1 #####
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
                 data.frame(sp, rownames(d[1:20,]), d[1:20,"Rhat"]))
  # print(sp)
}
colnames(rhats) <- c("species", "param", "rhat")
rhats$rhat <- as.numeric(rhats$rhat)

### Take off the species with at least 1 Rhat > 1.1
species <-
  rhats %>%
  group_by(species) %>%
  filter(all(rhat <= 1.1, na.rm = T)) %>%
  distinct(species) %>% pull()

## Remove the per species models
rm(list = ls()[ls() %in% gsub("^summary_|\\.rds$","",list.files(path = "summary_models/"))])
################

## Histogram of Rhat distributions for all parameters for the 261 selected species ###
pdf("mixed_model/figures/rhats.pdf",
    width = 8.27, height = 5.83)
rhats %>%
  group_by(species) %>%
  filter(all(rhat <= 1.1, na.rm = T)) %>%
  mutate(param = stringr::str_replace(param, "mean", "alpha"),
         param = stringr::str_replace(param, "lam", "lambda")) %>%
  ggplot() +
  geom_histogram(aes(rhat), binwidth = 0.01) +
  facet_wrap(vars(param), scales = "fixed")+
  geom_vline(xintercept = 1.1, linetype = "dashed", color = "red")+
  scale_x_continuous(breaks = seq(1, 1.1, by = 0.05)) +
  theme_bw()
dev.off()

species_data <- data.frame(species = as.character(),
                           ab_trend = as.numeric(),
                           g = as.numeric(),
                           g_upper = as.numeric(),
                           g_lower = as.numeric(),
                           absg = as.numeric(),
                           absg_upper = as.numeric(),
                           absg_lower = as.numeric())

for(i in 1:length(perspecies_N_output$mean$beta1)){

  s = species[i]
  a = perspecies_N_output$mean$beta1[i]
  a_upper = perspecies_N_output$q97.5$beta1[i]
  a_lower = perspecies_N_output$q2.5$beta1[i]
  g = perspecies_g_output$mean$beta1[i]
  g_upper = perspecies_g_output$q97.5$beta1[i]
  g_lower = perspecies_g_output$q2.5$beta1[i]
  absg = perspecies_absg_output$mean$beta1[i]
  absg_upper = perspecies_absg_output$q97.5$beta1[i]
  absg_lower = perspecies_absg_output$q2.5$beta1[i]

  species_data <- rbind(species_data,
                        data.frame(species = s,
                                   ab_trend = a,
                                   a_upper = a_upper,
                                   a_lower = a_lower,
                                   g = g,
                                   g_upper = g_upper,
                                   g_lower = g_lower,
                                   absg = absg,
                                   absg_upper = absg_upper,
                                   absg_lower = absg_lower))

}

## Add the color palette
tmp <-
species_data %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))
  ## remove one species that increase the scale of the plot
  # filter(species != "Eurasian Collared-Dove")


## delta N vs delta g per species
tmp <- tmp %>%
  mutate(color = case_when(
    a_upper < 0 & g_upper < 0 ~ "#E31A1C",
    a_upper < 0 & g_lower > 0 ~ "#FF7F00",
    a_lower > 0 & g_lower > 0 ~ "#33A02C",
    a_lower > 0 & g_upper < 0 ~ "#1F78B4"
  ))

svg("mixed_model/figures/deltaN_vs_deltag_sp.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = ab_trend, y = g,
                    ymin = g_lower, ymax = g_upper, colour = color),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = g,
                    xmin = a_lower, xmax = a_upper, colour= color),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_point(aes(ab_trend, g, color = color),
             show.legend = F, size = 1)+
  scale_colour_manual(values = unique(sort(tmp$color)))+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(ab_trend, g, label = species),
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10,  ## to control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  # scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_x_continuous(trans = ggallin::pseudolog10_trans)+
  scale_y_continuous(labels = scales::scientific)+
  theme_bw()

dev.off()

## delta N vs delta absg per species
tmp <- tmp %>%
  mutate(color_absg = case_when(
    a_upper < 0 & absg_upper < 0 ~ "#E31A1C",
    a_upper < 0 & absg_lower > 0 ~ "#FF7F00",
    a_lower > 0 & absg_lower > 0 ~ "#33A02C",
    a_lower > 0 & absg_upper < 0 ~ "#1F78B4"
  ))

svg("mixed_model/figures/deltaN_vs_deltaabsg_sp.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    ymin = absg_lower, ymax = absg_upper, colour = color_absg),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    xmin = a_lower, xmax = a_upper, colour= color_absg),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_point(aes(ab_trend, absg, color = color_absg),
             show.legend = F, size = 1)+
  scale_colour_manual(values = unique(sort(tmp$color_absg)))+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(ab_trend, absg, label = species),
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10,  ## to control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  # scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_x_continuous(trans = ggallin::pseudolog10_trans)+
  scale_y_continuous(trans = ggallin::pseudolog10_trans)+
  theme_bw()

dev.off()

## Abundance and growth rate trends per species
pdf("mixed_model/figures/trends_species.pdf",
    width = 5.83, height = 4.13)
for(metric in c("N", "g", "absg")){

  d <- readRDS(paste0("mixed_model/save_samples/mean_sd_",metric,"_sp.rds"))

  d <- d %>% select(-paste0("sd_",metric)) %>% pivot_wider(names_from = year, values_from = paste0("mean_",metric))

  print(d %>%
          as.data.frame() %>%
          pivot_longer(cols = -species, names_to = "year", values_to = "value") %>%
          mutate(year = as.numeric(year)) %>%
          ggplot()+
          geom_line(aes(year, value, group = species), color = "grey")+
          geom_hline(yintercept = ifelse(metric %in% c("gt0", "g", "absg"), 0, NA_real_), linetype="dashed", color = "black")+
          stat_summary(aes(year, value), fun = median, geom = "line", color = "blue",
                       linewidth = 1, linetype = 2)+
          geom_hline(yintercept = ifelse(metric %in% c("N", "absg"), median(pull(d[,2])), NA), linetype="dashed", color = "black")+
          scale_y_continuous(trans= if (metric %in% c("g","absg")) ggallin::ssqrt_trans else "identity")+
          ylab(metric)+
          theme_light())

}
dev.off()

## Per family plots #####
## Load the models
perfamily_N_output <- readRDS("mixed_model/outputs/perfamily_N_output.rds")
perfamily_g_output <- readRDS("mixed_model/outputs/perfamily_g_output.rds")
perfamily_absg_output <- readRDS("mixed_model/outputs/perfamily_absg_output.rds")
## Get back the family names
family <- readr::read_csv("data/phylo.csv")

## Filter the species with good Rhat
family = family[family$English_Common_Name %in% species,]

## Extract family names
families <- unique(family$Family)

families_data <- data.frame(family = as.character(),
                            ab_trend = as.numeric(),
                            a_upper = as.numeric(),
                            a_lower = as.numeric(),
                            g = as.numeric(),
                            g_upper = as.numeric(),
                            g_lower = as.numeric(),
                            absg = as.numeric(),
                            absg_upper = as.numeric(),
                            absg_lower = as.numeric())

for(i in 1:length(perfamily_N_output$mean$beta1)){

  f = families[i]
  a = perfamily_N_output$mean$beta1[i]
  a_upper = perfamily_N_output$q97.5$beta1[i]
  a_lower = perfamily_N_output$q2.5$beta1[i]
  g = perfamily_g_output$mean$beta1[i]
  g_upper = perfamily_g_output$q97.5$beta1[i]
  g_lower = perfamily_g_output$q2.5$beta1[i]
  absg = perfamily_absg_output$mean$beta1[i]
  absg_upper = perfamily_absg_output$q97.5$beta1[i]
  absg_lower = perfamily_absg_output$q2.5$beta1[i]

  families_data <- rbind(families_data,
                        data.frame(family = f,
                                   ab_trend = a,
                                   a_upper = a_upper,
                                   a_lower = a_lower,
                                   g = g,
                                   g_upper = g_upper,
                                   g_lower = g_lower,
                                   absg = absg,
                                   absg_upper = absg_upper,
                                   absg_lower = absg_lower))

}

## Add the color palette
tmp <-
  families_data %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))

## delta N vs delta g per family
tmp <- tmp %>%
  mutate(color = case_when(
    a_upper < 0 & g_upper < 0 ~ "#E31A1C",
    a_upper < 0 & g_lower > 0 ~ "#FF7F00",
    a_lower > 0 & g_lower > 0 ~ "#33A02C",
    a_lower > 0 & g_upper < 0 ~ "#1F78B4"
  ))

svg("mixed_model/figures/deltaN_vs_deltag_fam.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = ab_trend, y = g,
                    ymin = g_lower, ymax = g_upper, colour = color),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = g,
                    xmin = a_lower, xmax = a_upper, colour= color),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_point(aes(ab_trend, g, color = color),
             show.legend = F, size = 1)+
  scale_colour_manual(values = unique(sort(tmp$color)))+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(ab_trend, g, label = family), check_overlap = T,
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10,  ## too control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  # scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_x_continuous(trans = ggallin::pseudolog10_trans)+
  scale_y_continuous(labels = scales::scientific)+
  theme_bw()

dev.off()

## delta N vs delta absg per family
tmp <- tmp %>%
  mutate(color_absg = case_when(
    a_upper < 0 & absg_upper < 0 ~ "#E31A1C",
    a_upper < 0 & absg_lower > 0 ~ "#FF7F00",
    a_lower > 0 & absg_lower > 0 ~ "#33A02C",
    a_lower > 0 & absg_upper < 0 ~ "#1F78B4"
  ))

svg("mixed_model/figures/deltaN_vs_deltaabsg_fam.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    ymin = absg_lower, ymax = absg_upper, colour = color_absg),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    xmin = a_lower, xmax = a_upper, colour= color_absg),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_point(aes(ab_trend, absg, color = color_absg),
             show.legend = F, size = 1)+
  scale_colour_manual(values = unique(sort(tmp$color_absg)))+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(ab_trend, absg, label = family), check_overlap = T,
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10,  ## too control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  # scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_x_continuous(trans = ggallin::pseudolog10_trans)+
  scale_y_continuous(trans = ggallin::pseudolog10_trans)+
  theme_bw()

dev.off()

for(metric in c("N", "g", "absg")){

  d <- readRDS(paste0("mixed_model/save_samples/mean_",metric,"_families.rds"))
  sd <- readRDS(paste0("mixed_model/save_samples/sd_",metric,"_families.rds"))

  d <- d %>% as.data.frame()
  colnames(d) <- c(((2021-ncol(d))+1):2021)
  sd <- sd %>% as.data.frame()
  colnames(sd) <- c(((2021-ncol(d))+1):2021)

  print(d %>% rownames_to_column(var = "family") %>%
          pivot_longer(cols = -family, names_to = "year", values_to = "values") %>%
          left_join(
            sd %>% rownames_to_column(var = "family") %>%
              pivot_longer(cols = -family, names_to = "year", values_to = "sd")
          ) %>%
          mutate(year = as.numeric(year)) %>%
          ggplot()+
          ## add 95% credible interval
          geom_ribbon(aes(year, values, group = family,
                          ymin = values-sd*1.96,
                          ymax= values+sd*1.96), colour="grey", fill="grey",alpha=0.5)+
          geom_line(aes(year, values, group = family, color = family), show.legend = F)+
          facet_wrap(vars(family), scales = "free")+
          geom_hline(yintercept = if (metric %in% c("g","absg")) 0 else NULL, linetype="dashed", color = "black")+
          ylab(metric)+
          theme_light())

}

## Per habitats plots #####
## Load the models
perhabitat_N_output <- readRDS("mixed_model/outputs/perhabitat_N_output.rds")
perhabitat_g_output <- readRDS("mixed_model/outputs/perhabitat_g_output.rds")
perhabitat_absg_output <- readRDS("mixed_model/outputs/perhabitat_absg_output.rds")
## Get back the family names
mean_N_habitats <- readRDS("mixed_model/save_samples/mean_N_habitats.rds")

## Extract family names
habitats <- rownames(mean_N_habitats)

habitats_data <- data.frame(habitat = as.character(),
                            ab_trend = as.numeric(),
                            a_upper = as.numeric(),
                            a_lower = as.numeric(),
                            g = as.numeric(),
                            g_upper = as.numeric(),
                            g_lower = as.numeric(),
                            absg = as.numeric(),
                            absg_upper = as.numeric(),
                            absg_lower = as.numeric())

for(i in 1:length(perhabitat_N_output$mean$beta1)){

  h = habitats[i]
  a = perhabitat_N_output$mean$beta1[i]
  a_upper = perhabitat_N_output$q97.5$beta1[i]
  a_lower = perhabitat_N_output$q2.5$beta1[i]
  g = perhabitat_g_output$mean$beta1[i]
  g_upper = perhabitat_g_output$q97.5$beta1[i]
  g_lower = perhabitat_g_output$q2.5$beta1[i]
  absg = perhabitat_absg_output$mean$beta1[i]
  absg_upper = perhabitat_absg_output$q97.5$beta1[i]
  absg_lower = perhabitat_absg_output$q2.5$beta1[i]

  habitats_data <- rbind(habitats_data,
                         data.frame(habitat = h,
                                    ab_trend = a,
                                    a_upper = a_upper,
                                    a_lower = a_lower,
                                    g = g,
                                    g_upper = g_upper,
                                    g_lower = g_lower,
                                    absg = absg,
                                    absg_upper = absg_upper,
                                    absg_lower = absg_lower))

}

## Add the color palette
tmp <-
  habitats_data %>%
  mutate(nSign = ifelse(ab_trend < 0, "decrease", "increase"))

## delta N vs delta g per habitat
tmp <- tmp %>%
  mutate(color = case_when(
    a_upper < 0 & g_upper < 0 ~ "#E31A1C",
    a_upper < 0 & g_lower > 0 ~ "#FF7F00",
    a_lower > 0 & g_lower > 0 ~ "#33A02C",
    a_lower > 0 & g_upper < 0 ~ "#1F78B4"
  ))

svg("mixed_model/figures/deltaN_vs_deltag_hab.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = ab_trend, y = g,
                    ymin = g_lower, ymax = g_upper, colour = color),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = g,
                    xmin = a_lower, xmax = a_upper, colour= color),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_point(aes(ab_trend, g, color = color),
             show.legend = F, size = 1)+
  scale_colour_manual(values = unique(sort(tmp$color)))+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(ab_trend, g, label = habitat), check_overlap = T,
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10,  ## too control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  # scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_x_continuous(trans = ggallin::pseudolog10_trans)+
  scale_y_continuous(labels = scales::scientific)+
  theme_bw()
dev.off()

## delta N vs delta absg per habitat
tmp <- tmp %>%
  mutate(color_absg = case_when(
    a_upper < 0 & absg_upper < 0 ~ "#E31A1C",
    a_upper < 0 & absg_lower > 0 ~ "#FF7F00",
    a_lower > 0 & absg_lower > 0 ~ "#33A02C",
    a_lower > 0 & absg_upper < 0 ~ "#1F78B4"
  ))

svg("mixed_model/figures/deltaN_vs_deltaabsg_hab.svg",
    height = 4.13, width = 5.83)

tmp %>%
  ggplot()+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    ymin = absg_lower, ymax = absg_upper, colour = color_absg),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    xmin = a_lower, xmax = a_upper, colour= color_absg),
                alpha = .3, size = .5,
                show.legend = F)+
  geom_point(aes(ab_trend, absg, color = color_absg),
             show.legend = F, size = 1)+
  scale_colour_manual(values = unique(sort(tmp$color_absg)))+
  geom_hline(yintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_vline(xintercept = 0, linetype="dashed", color = "#595959ff", size = 1)+
  geom_text_repel(aes(ab_trend, absg, label = habitat), check_overlap = T,
                  color = "black",
                  min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10,  ## too control the number of labels
                  segment.size = 1,
                  fontface = "bold",
                  size = 3)+
  # scale_color_manual(values = unique(sort(tmp$rec_vs_loss)))+
  scale_x_continuous(trans = ggallin::pseudolog10_trans)+
  scale_y_continuous(trans = ggallin::pseudolog10_trans)+
  theme_bw()
dev.off()


## Trends per habitat
pdf("mixed_model/figures/trends_habitat.pdf",
    width = 8.27, height = 5.83)
for(metric in c("N", "g", "absg")){

  d <- readRDS(paste0("mixed_model/save_samples/mean_",metric,"_habitats.rds"))
  sd <- readRDS(paste0("mixed_model/save_samples/sd_",metric,"_habitats.rds"))

  d <- d %>% as.data.frame()
  colnames(d) <- c(((2021-ncol(d))+1):2021)
  sd <- sd %>% as.data.frame()
  colnames(sd) <- c(((2021-ncol(d))+1):2021)

  print(d %>% rownames_to_column(var = "habitat") %>%
          pivot_longer(cols = -habitat, names_to = "year", values_to = "values") %>%
          left_join(
            sd %>% rownames_to_column(var = "habitat") %>%
              pivot_longer(cols = -habitat, names_to = "year", values_to = "sd")
          ) %>%
          mutate(year = as.numeric(year)) %>%
          ggplot()+
          ## add 95% credible interval
          geom_ribbon(aes(year, values, group = habitat,
                          ymin = values-sd*1.96,
                          ymax= values+sd*1.96), colour="grey", fill="grey",alpha=0.5)+
          geom_line(aes(year, values, group = habitat, color = habitat), show.legend = F)+
          facet_wrap(vars(habitat), scales = "free")+
          geom_hline(yintercept = if (metric %in% c("g","absg")) 0 else NULL, linetype="dashed", color = "black")+
          ylab(metric)+
          theme_light())

}

dev.off()

## Add 95% credible interval for abundance change
routes_shp$N_lower <- overall_N_output$q2.5$beta1
routes_shp$N_upper <- overall_N_output$q97.5$beta1

## Plot of change in abundance vs. change in growth rate
pdf("mixed_model/figures/deltan_vs_deltag.pdf",
    width = 8.27, height = 5.83)

df <-
  routes_shp %>%
  select(ab_trend, growth_rate, g_upper, g_lower, N_upper, N_lower) %>%
  mutate(data = "Not smoothed") %>% st_drop_geometry() %>%
  mutate(color = case_when(
    N_upper < 0 & g_upper < 0 ~ "#E31A1C",
    N_upper < 0 & g_lower > 0 ~ "#FF7F00",
    N_lower > 0 & g_lower > 0 ~ "#33A02C",
    N_lower > 0 & g_upper < 0 ~ "#1F78B4"
  ))

p1 <-
  df %>%
  ggplot()+
  geom_point(aes(ab_trend, growth_rate, colour = color), show.legend = F)+
  scale_colour_manual(values = unique(sort(df$color)))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_errorbar(aes(x = ab_trend, y = growth_rate,
                    ymin = g_lower, ymax = g_upper, colour = color),
                alpha = .3, size = .2,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = growth_rate,
                    xmin = N_lower, xmax = N_upper, colour = color),
                alpha = .3, size = .2,
                show.legend = F)+
  ylab("Growth rate change")+
  xlab("Abundance change")+
  facet_wrap(. ~ data, scales = "free")+
  theme_bw()


df <-
  routes_shp %>% select(ab_trend_gam, growth_rate_gam) %>%
  rename(ab_trend = ab_trend_gam,
         growth_rate = growth_rate_gam) %>%
  mutate(data = "Smoothed") %>% st_drop_geometry() %>%
  mutate(color = case_when(
    ab_trend < 0 & growth_rate < 0 ~ "#E31A1C",
    ab_trend < 0 & growth_rate > 0 ~ "#FF7F00",
    ab_trend > 0 & growth_rate > 0 ~ "#33A02C",
    ab_trend > 0 & growth_rate < 0 ~ "#1F78B4"
  ))


p2 <-
  df %>%
  ggplot()+
  geom_point(aes(ab_trend, growth_rate, colour = color), show.legend = F)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_colour_manual(values = unique(sort(df$color)))+
  facet_wrap(. ~ data, scales = "free")+
  xlab("Abundance change")+
  ylab("Growth rate change")+
  theme_bw()

cowplot::plot_grid(p1, p2)

dev.off()

## Plot of change in abundance vs. change in absolute growth rate
routes_shp$absg_upper <- overall_absg_output$q97.5$beta1
routes_shp$absg_lower <- overall_absg_output$q2.5$beta1

pdf("mixed_model/figures/deltan_vs_deltaabsg.pdf",
    width = 8.27, height = 5.83)

df <-
  routes_shp %>%
  select(ab_trend, absg, absg_upper, absg_lower, N_upper, N_lower) %>%
  mutate(data = "Not smoothed") %>% st_drop_geometry() %>%
  mutate(color_absg = case_when(
    N_upper < 0 & absg_upper < 0 ~ "#E31A1C",
    N_upper < 0 & absg_lower > 0 ~ "#FF7F00",
    N_lower > 0 & absg_lower > 0 ~ "#33A02C",
    N_lower > 0 & absg_upper < 0 ~ "#1F78B4"
  ))

p1 <-
  df %>%
  ggplot()+
  geom_point(aes(ab_trend, absg, colour = color_absg), show.legend = F)+
  scale_colour_manual(values = unique(sort(df$color_absg)))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    ymin = absg_lower, ymax = absg_upper, colour = color_absg),
                alpha = .3, size = .2,
                show.legend = F)+
  geom_errorbar(aes(x = ab_trend, y = absg,
                    xmin = N_lower, xmax = N_upper, colour = color_absg),
                alpha = .3, size = .2,
                show.legend = F)+
  ylab("Absolute growth rate change")+
  xlab("Abundance change")+
  facet_wrap(. ~ data, scales = "free")+
  theme_bw()


df <-
  routes_shp %>% select(ab_trend_gam, absg_gam) %>%
  rename(ab_trend = ab_trend_gam,
         absg = absg_gam) %>%
  mutate(data = "Smoothed") %>% st_drop_geometry() %>%
  mutate(color_absg = case_when(
    ab_trend < 0 & absg < 0 ~ "#E31A1C",
    ab_trend < 0 & absg > 0 ~ "#FF7F00",
    ab_trend > 0 & absg > 0 ~ "#33A02C",
    ab_trend > 0 & absg < 0 ~ "#1F78B4"
  ))


p2 <-
  df %>%
  ggplot()+
  geom_point(aes(ab_trend, absg, colour = color_absg), show.legend = F)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_colour_manual(values = unique(sort(df$color_absg)))+
  facet_wrap(. ~ data, scales = "free")+
  xlab("Abundance change")+
  ylab("Absolute growth rate change")+
  theme_bw()

cowplot::plot_grid(p1, p2)

dev.off()

