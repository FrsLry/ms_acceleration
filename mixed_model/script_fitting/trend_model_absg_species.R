## Fitting the trend models
library(jagsUI)
library(dplyr)
library(tidyr)
library(tibble)
## Load the data
sp_absg <- readRDS("mixed_model/save_samples/mean_sd_absg_sp.rds")

## Format the per species data
mean_absg_sp <-
  sp_absg[,!(colnames(sp_absg) == "sd_absg")] %>%
  pivot_wider(names_from = "year", values_from = "mean_absg") %>%
  column_to_rownames("species")

sd_absg_sp <-
  sp_absg[,!(colnames(sp_absg) == "mean_absg")] %>%
  pivot_wider(names_from = "year", values_from = "sd_absg") %>%
  column_to_rownames("species")


## JAGS data
data_list <- list(absg = as.matrix(mean_absg_sp),
                  sd = as.matrix(sd_absg_sp),
                  species = nrow(mean_absg_sp),
                  years = ncol(mean_absg_sp),
                  t = 1:ncol(mean_absg_sp))

## JAGS model
cat(
  "model{

  # Priors
  grand.beta0 ~ dnorm(0, 1.0E-4)
  grand.beta1 ~ dnorm(0, 1.0E-4)

  # Group SDs: half-t(df=3, scale=1.5)
  grand.sigma.beta0 ~ dt(0, 1/pow(1.5,2), 3) T(0,)
  grand.sigma.beta1 ~ dt(0, 1/pow(1.5,2), 3) T(0,)

  for(i in 1:species){
    z0[i] ~ dnorm(0, 1)
    z1[i] ~ dnorm(0, 1)
    beta0[i] <- grand.beta0 + z0[i] * grand.sigma.beta0
    beta1[i] <- grand.beta1 + z1[i] * grand.sigma.beta1
  }

  # Residual SD (keep as in your script if you want; shown here as Uniform)
  sigma ~ dunif(0, 1000000)
  tau = 1/(sigma*sigma)

  # Likelihood
  for(sp in 1:species){
    for(year in 1:years){

      tau.obs[sp,year] <- 1/(sd[sp,year]*sd[sp,year])
      error[sp,year] ~ dnorm(0, tau.obs[sp,year])

      absg[sp,year] ~ dnorm(mu[sp,year], tau)
      mu[sp,year] <- beta0[sp] + beta1[sp] * t[year] + error[sp,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/absg_sp.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/absg_sp.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/perspecies_absg_output.rds")

