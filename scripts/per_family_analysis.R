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
  d_list[[which(species == sp)]] <- list("N"=N,"R_raw"=R_raw,"L_raw"=L_raw)
  names(d_list)[which(species == sp)] <- sp
}
## Keep only the species for which rhat are decent
d_list <- d_list[names(d_list) %in% ok_sp]
## Cleaning
rm(list = ls()[ls() %in% species])
gc()

## Create the lists of each family ##########
family <- readr::read_csv("data/phylo.csv")
hab <- readRDS("data/selected_species.rds")

family <- family %>%
  right_join(hab %>% rename(English_Common_Name = COM_NAME_updated)) %>%
  select(-AOU) %>%
  distinct()

pdf("figures/perfamily_allRoad.pdf",
    height = 8.27, width = 11.69)

per_family <- data.frame(year = as.numeric(),
                         value = as.numeric(),
                         metric = as.character(),
                         family = as.character())

## per family, all routes
for(fam in na.omit(unique(family$Family))){
  tmp_list <- d_list[names(d_list) %in% family$English_Common_Name[family$Family == fam]]
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
  # tmp_abPercent <- cbind(routes, tmp_abPercent[,-1])
  tmp_abPercent <- as.data.frame(tmp_abPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  tmp_abPercent[sapply(tmp_abPercent, is.nan)] <- 0
  tmp_abPercent[tmp_abPercent == "Inf"] <- 0

  colnames(tmp_lossPercent) <- c(1988:2021)
  # tmp_lossPercent <- cbind(routes, tmp_lossPercent[,-1])
  tmp_lossPercent <- as.data.frame(tmp_lossPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  tmp_lossPercent[sapply(tmp_lossPercent, is.nan)] <- 0
  tmp_lossPercent[tmp_lossPercent == "Inf"] <- 0

  colnames(tmp_recPercent) <- c(1988:2021)
  # tmp_recPercent <- cbind(routes, tmp_recPercent[,-1])
  tmp_recPercent <- as.data.frame(tmp_recPercent) %>% mutate_at(vars(matches("[0-9]")), function(x) as.numeric(x))
  tmp_recPercent[sapply(tmp_recPercent, is.nan)] <- 0
  tmp_recPercent[tmp_recPercent == "Inf"] <- 0

  ## Create the dataset to compare all the habitats together
  per_family <- rbind(per_family,
                      as.data.frame(tmp_ab) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "abundance") %>%
                        rbind(as.data.frame(tmp_recRaw) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "recruitment")) %>%
                        rbind(as.data.frame(tmp_lossRaw) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "loss")) %>%
                        rbind(as.data.frame(tmp_abPercent) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "growth rate")) %>%
                        rbind(as.data.frame(tmp_recPercent) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "recruitment rate")) %>%
                        rbind(as.data.frame(tmp_lossPercent) %>% pivot_longer(everything(), names_to = "year", values_to = "value") %>% mutate(metric = "loss rate")) %>%
                        mutate(family = paste0(fam, " (s = ", length(tmp_list), ")")))


}

## Compute trends for each family at the continent scale
trends <-
  per_family %>%
  mutate(year = as.numeric(year)) %>%
  group_by(metric, family) %>%
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
  geom_col(aes(reorder(family, -slope), slope),
           show.legend = F)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("Slope")+
  coord_flip()+
  theme_minimal()+
  facet_wrap(~metric, scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

trends %>%
  filter(metric %in% c("growth rate", "loss rate", "recruitment rate")) %>%
  ### brut forcing orders #####
  mutate(metric = factor(metric, levels = c("growth rate", "loss rate", "recruitment rate"))) %>%
  arrange(metric, -slope) %>%
  mutate(family = factor(family, levels = unique(.$family))) %>%
  ################
  ggplot()+
  geom_col(aes(family, slope),
           show.legend = F)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("Slope")+
  coord_flip()+
  theme_minimal()+
  facet_wrap(~metric, scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())

dev.off()

pdf("figures/per_family_trends.pdf",
    width = 5.83, height = 4.13)

print(
  per_family %>%
    filter(!grepl("rate$", .$metric)) %>%
    filter(grepl("abundance", .$metric)) %>%
    mutate(year = as.numeric(year), value = as.numeric(value)) %>%
    group_by(year, family) %>% summarize(abundance = sum(value)) %>%
    ggplot()+
    geom_line(aes(year, abundance, group = family), color = "grey")+
    stat_summary(aes(year, abundance), fun = median, geom = "line", color = "blue", size = 1, linetype = 2)+
    theme_light()
)

print(
  per_family %>%
    filter(grepl("growth rate", .$metric)) %>%
    mutate(year = as.numeric(year), value = as.numeric(value)) %>%
    group_by(year, family) %>% summarize(growth_rate = sum(value)) %>%
    ggplot()+
    geom_line(aes(year, growth_rate, group = family), color = "grey")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    scale_y_continuous(trans=pseudolog10_trans)+
    stat_summary(aes(year, growth_rate), fun = median, geom = "line", color = "blue", size = 1, linetype = 2)+
    theme_light()
)

dev.off()

svg("figures/rec_vs_loss_family.svg",
    height = 4.13, width = 5.83)

tmp <-
  trends %>%
  filter(metric %in% c("loss rate", "recruitment rate", "abundance")) %>%
  select(metric, family, slope) %>%
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
             shape = ifelse(tmp$nSign == "increase",  "\u2191", "\u2193"),
             show.legend = F, size = 4)+
  geom_abline(slope = 1, linetype = "longdash", color = "#595959ff")+
  geom_abline(slope = -1, linetype = "longdash", color = "#595959ff")+
  geom_hline(yintercept = 0, linetype="dashed", size = 1, color = "#595959ff")+
  geom_vline(xintercept = 0, linetype="dashed", size = 1, color = "#595959ff")+
  ggrepel::geom_text_repel(aes(`loss rate`, `recruitment rate`, label = family), check_overlap = T,
                           color = "black",
                           min.segment.length = unit(0, 'lines'),
                           max.overlaps = 4,  ## too control the number of labels
                           segment.size = 1,
                           fontface = "bold",
                           size = 3)+
  xlim(-max(abs(trends$slope[trends$metric == "loss rate"])), max(abs(trends$slope[trends$metric == "loss rate"])))+
  ylim(-max(abs(trends$slope[trends$metric == "recruitment rate"])), max(abs(trends$slope[trends$metric == "recruitment rate"])))+
  scale_y_continuous(trans= ssqrt_trans, labels = scales::scientific)+
  scale_x_continuous(trans= ssqrt_trans, labels = scales::scientific)+
  scale_color_manual(values = unique(sort(tmp$color)))+
  theme_bw()

dev.off()

pdf("figures/boxplot_perfamily.pdf",
    height = 2, width = 5.83)
tmp %>%
  ungroup() %>%
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
  geom_boxplot(aes(reorder(color, order), abundance, fill = color), show.legend = F, width = 0.1, outlier.shape = NA)+
  geom_violin(aes(reorder(color, order), abundance, fill = color), show.legend = F)+
  ggforce::geom_sina(aes(reorder(color, order), abundance), show.legend = F, alpha = .5, size = 1)+
  geom_boxplot(aes(reorder(color, order), abundance), fill = NA, show.legend = F, width = 0.1,color = "darkgrey", outlier.shape = NA)+
  scale_fill_manual(values = unique(sort(tmp$color)))+
  scale_y_continuous(trans= ggallin::pseudolog10_trans, breaks = c(-10000,-5000,-10,0,10,1000))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title  = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11))
dev.off()
