## Fitting the trend models
library(jagsUI)
library(dplyr)
library(tidyr)
library(tibble)
## Load the data
sp_r <- readRDS("mixed_model/save_samples/mean_sd_rr_sp.rds")

## For the very unique case where a value is missing (only 1 year for 1 Eurasian-Collared Dove)
## Use the mean over the time series
for(sp in unique(sp_r$species)){
  if(all(!is.na(sp_r[sp_r$species == sp,]$mean_R))){
    next
  }else{
    sp_r[sp_r$species == sp,]$mean_R[is.na(sp_r[sp_r$species == sp,]$mean_R)] <-
      mean(sp_r[sp_r$species == sp,]$mean_R[!is.na(sp_r[sp_r$species == sp,]$mean_R)])
  }
}

for(sp in unique(sp_r$species)){
  if(all(sp_r[sp_r$species == sp,]$mean_R != "Inf")){
    next
  }else{
    sp_r[sp_r$species == sp,]$mean_R[sp_r[sp_r$species == sp,]$mean_R == "Inf"] <-
      mean(sp_r[sp_r$species == sp,]$mean_R[sp_r[sp_r$species == sp,]$mean_R != "Inf"])
  }
}


for(sp in unique(sp_r$species)){
  if(all(!is.na(sp_r[sp_r$species == sp,]$sd_R))){
    next
  }else{
    sp_r[sp_r$species == sp,]$sd_R[is.na(sp_r[sp_r$species == sp,]$sd_R)] <-
      mean(sp_r[sp_r$species == sp,]$sd_R[!is.na(sp_r[sp_r$species == sp,]$sd_R)])
  }
}

for(sp in unique(sp_r$species)){
  if(all(sp_r[sp_r$species == sp,]$sd_R != "Inf")){
    next
  }else{
    sp_r[sp_r$species == sp,]$sd_R[sp_r[sp_r$species == sp,]$sd_R == "Inf"] <-
      mean(sp_r[sp_r$species == sp,]$sd_R[sp_r[sp_r$species == sp,]$sd_R != "Inf"])
  }
}

## Format the per species data
mean_r_sp <-
  sp_r[,!(colnames(sp_r) == "sd_R")] %>%
    pivot_wider(names_from = "year", values_from = "mean_R") %>%
    column_to_rownames("species")

sd_r_sp <-
  sp_r[,!(colnames(sp_r) == "mean_R")] %>%
  pivot_wider(names_from = "year", values_from = "sd_R") %>%
  column_to_rownames("species")

## JAGS data
data_list <- list(rr = as.matrix(mean_r_sp),
          sd = as.matrix(sd_r_sp),
          species = nrow(mean_r_sp),
          years = ncol(mean_r_sp),
          t = 1:ncol(mean_r_sp))

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

      rr[sp,year] ~ dnorm(mu[sp,year], tau)
      mu[sp,year] <- beta0[sp] + beta1[sp] * t[year] + error[sp,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/r_sp.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/r_sp.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/perspecies_r_output.rds")
