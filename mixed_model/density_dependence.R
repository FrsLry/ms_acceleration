## Check for density dependence
library(ggplot2)
library(tidyr)
library(dplyr)
## Load the data
sp_N <- readRDS("mixed_model/save_samples/mean_sd_N_sp.rds")
sp_g <- readRDS("mixed_model/save_samples/mean_sd_g_sp.rds")
sp_r <- readRDS("mixed_model/save_samples/mean_sd_rr_sp.rds")

## Create the data frame
df_g <- data.frame(species = as.character(),
                   slope = as.numeric(),
                   pval = as.numeric())

## For each species, check growth rate density dependence by fitting a linear model of growth rate
## as a function of abundance
for(sp in unique(sp_N$species)[unique(sp_N$species) != "Eurasian Collared-Dove"]){

  lm <- lm(sp_g$mean_g[sp_g$sp == sp] ~ sp_N$mean_N[sp_N$species == sp][-35])
  df_g[which(unique(sp_N$species)[unique(sp_N$species) != "Eurasian Collared-Dove"] == sp),] <- c(sp, lm$coefficients[2], summary(lm)$coefficients[2,4])

}

df_g$type <- "Density dependence of growth rate"

## Create the data frame
df_r <- data.frame(species = as.character(),
                   slope = as.numeric(),
                   pval = as.numeric())

## For each species, check recruitment rate density dependence by fitting a linear model of recruitment rate
## as a function of abundance
for(sp in unique(sp_N$species)[unique(sp_N$species) != "Eurasian Collared-Dove"]){

  lm <- lm(sp_r$mean_R[sp_r$sp == sp] ~ sp_N$mean_N[sp_N$species == sp][-35])
  df_r[which(unique(sp_N$species)[unique(sp_N$species) != "Eurasian Collared-Dove"] == sp),] <- c(sp, lm$coefficients[2], summary(lm)$coefficients[2,4])

}

df_r$type <- "Density dependence of recruitment rate"

## Plot the slopes vs. p-values of these lm
pdf("mixed_model/figures/density_dpendence.pdf",
    width = 5.83, height = 4.13)
df_g %>% rbind(df_r) %>%
  mutate(slope = as.numeric(slope),
         pval = as.numeric(pval)) %>%
  mutate(color = ifelse(.$slope < 0 & .$pval < 0.05, "Negative & Significant", "Other")) %>%
  mutate(color = ifelse(.$slope > 0 & .$pval < 0.05, "Positive & Significant", color)) %>%
  ggplot()+
  geom_point(aes(slope, pval, colour = color), show.legend = F)+
  scale_x_continuous(trans= ggallin::ssqrt_trans)+
  scale_y_continuous(trans= ggallin::ssqrt_trans)+
  geom_hline(yintercept = 0.05, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  scale_color_manual(values = c("Negative & Significant" = "red",
                                "Positive & Significant" = "blue",
                                "Other" = "grey"))+
  facet_wrap(vars(type))+
  theme_bw()
dev.off()

