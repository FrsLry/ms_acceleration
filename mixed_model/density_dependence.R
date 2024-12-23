## Check for density dependence
library(ggplot2)
library(tidyr)
library(dplyr)
## Load the data
sp_N <- readRDS("mixed_model/save_samples/mean_sd_N_sp.rds")
sp_g <- readRDS("mixed_model/save_samples/mean_sd_g_sp.rds")

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

## Plot the slopes vs. p-values of these lm
pdf("mixed_model/figures/density_dependence.pdf",
    width = 5.83, height = 4.13)

df_g %>%
  mutate(slope = as.numeric(slope),
         pval = as.numeric(pval)) %>%
  ggplot()+
  geom_density(aes(slope, y=..scaled..), alpha = 0.5, fill = "black")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = median(as.numeric(df_g$slope)), color = "red")+
  theme_bw()

dev.off()

