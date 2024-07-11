## Fitting the trend models
library(jagsUI)
library(dplyr)
library(tidyr)
library(tibble)
## Load the data
sp_N <- readRDS("mixed_model/save_samples/mean_sd_N_sp.rds")

## Format the per species data
mean_N_sp <-
  sp_N[,!(colnames(sp_N) == "sd_N")] %>%
    pivot_wider(names_from = "year", values_from = "mean_N") %>%
    column_to_rownames("species")

sd_N_sp <-
  sp_N[,!(colnames(sp_N) == "mean_N")] %>%
  pivot_wider(names_from = "year", values_from = "sd_N") %>%
  column_to_rownames("species")

## JAGS data
data_list <- list(N = as.matrix(mean_N_sp),
          sd = as.matrix(sd_N_sp),
          species = nrow(mean_N_sp),
          years = ncol(mean_N_sp),
          t = 1:ncol(mean_N_sp))

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

      N[sp,year] ~ dnorm(mu[sp,year], tau)
      mu[sp,year] <- beta0[sp] + beta1[sp] * t[year] + error[sp,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/N_sp.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/N_sp.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/perspecies_N_output.rds")
