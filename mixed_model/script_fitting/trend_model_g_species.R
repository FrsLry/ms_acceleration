## Fitting the trend models
library(jagsUI)
library(dplyr)
library(tidyr)
library(tibble)
## Load the data
sp_g <- readRDS("mixed_model/save_samples/mean_sd_g_sp.rds")

## For the very unique case where a value is missing (only 1 year for 1 Eurasian-Collared Dove)
## Use the mean over the time series
for(sp in unique(sp_g$species)){
  if(all(!is.na(sp_g[sp_g$species == sp,]$mean_g))){
    next
  }else{
    sp_g[sp_g$species == sp,]$mean_g[is.na(sp_g[sp_g$species == sp,]$mean_g)] <-
      mean(sp_g[sp_g$species == sp,]$mean_g[!is.na(sp_g[sp_g$species == sp,]$mean_g)])
  }
}

for(sp in unique(sp_g$species)){
  if(all(sp_g[sp_g$species == sp,]$mean_g != "Inf")){
    next
  }else{
    sp_g[sp_g$species == sp,]$mean_g[sp_g[sp_g$species == sp,]$mean_g == "Inf"] <-
      mean(sp_g[sp_g$species == sp,]$mean_g[sp_g[sp_g$species == sp,]$mean_g != "Inf"])
  }
}


for(sp in unique(sp_g$species)){
  if(all(!is.na(sp_g[sp_g$species == sp,]$sd_g))){
    next
  }else{
    sp_g[sp_g$species == sp,]$sd_g[is.na(sp_g[sp_g$species == sp,]$sd_g)] <-
      mean(sp_g[sp_g$species == sp,]$sd_g[!is.na(sp_g[sp_g$species == sp,]$sd_g)])
  }
}

for(sp in unique(sp_g$species)){
  if(all(sp_g[sp_g$species == sp,]$sd_g != "Inf")){
    next
  }else{
    sp_g[sp_g$species == sp,]$sd_g[sp_g[sp_g$species == sp,]$sd_g == "Inf"] <-
      mean(sp_g[sp_g$species == sp,]$sd_g[sp_g[sp_g$species == sp,]$sd_g != "Inf"])
  }
}

## Format the per species data
mean_g_sp <-
  sp_g[,!(colnames(sp_g) == "sd_g")] %>%
    pivot_wider(names_from = "year", values_from = "mean_g") %>%
    column_to_rownames("species")

sd_g_sp <-
  sp_g[,!(colnames(sp_g) == "mean_g")] %>%
  pivot_wider(names_from = "year", values_from = "sd_g") %>%
  column_to_rownames("species")

## JAGS data
data_list <- list(g = as.matrix(mean_g_sp),
          sd = as.matrix(sd_g_sp),
          species = nrow(mean_g_sp),
          years = ncol(mean_g_sp),
          t = 1:ncol(mean_g_sp))

## JAGS model
cat(
  "model{

  # Priors
  for(i in 1:species){
    beta0[i] ~ dnorm(grand.beta0, grand.tau.beta0)
    beta1[i] ~ dnorm(grand.beta1, grand.tau.beta1)
  }

  sigma ~ dunif(0, 100)
  tau = 1/(sigma*sigma)

  grand.beta0 ~ dnorm(0, 1e-10)
  grand.sigma.beta0 ~ dunif(0, 100)
  grand.tau.beta0 = 1/(grand.sigma.beta0*grand.sigma.beta0)

  grand.beta1 ~ dnorm(0, 1e-10)
  grand.sigma.beta1 ~ dunif(0, 100)
  grand.tau.beta1 = 1/(grand.sigma.beta1*grand.sigma.beta1)

  # Likelihood
  for(sp in 1:species){
    for(year in 1:years){

      tau.obs[sp,year] <- 1/(sd[sp,year]*sd[sp,year])
      error[sp,year] ~ dnorm(0, tau.obs[sp,year])

      g[sp,year] ~ dnorm(mu[sp,year], tau)
      mu[sp,year] <- beta0[sp] + beta1[sp] * t[year] + error[sp,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/g_sp.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/g_sp.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/perspecies_g_output.rds")
