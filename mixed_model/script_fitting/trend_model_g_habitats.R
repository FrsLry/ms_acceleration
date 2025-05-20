## Fitting the trend models
library(jagsUI)
library(dplyr)
library(tidyr)
library(tibble)
## Load the data
mean_g_habitats <- readRDS("mixed_model/save_samples/mean_g_habitats.rds")
sd_g_habitats <- readRDS("mixed_model/save_samples/sd_g_habitats.rds")

## JAGS data
data_list <- list(g = as.matrix(mean_g_habitats),
                  sd = as.matrix(sd_g_habitats),
                  habitats = nrow(mean_g_habitats),
                  years = ncol(mean_g_habitats),
                  t = 1:ncol(mean_g_habitats))

## JAGS model
cat(
  "model{

  # Priors
  for(i in 1:habitats){
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
  for(hab in 1:habitats){
    for(year in 1:years){

      tau.obs[hab,year] <- 1/(sd[hab,year]*sd[hab,year])
      error[hab,year] ~ dnorm(0, tau.obs[hab,year])

      g[hab,year] ~ dnorm(mu[hab,year], tau)
      mu[hab,year] <- beta0[hab] + beta1[hab] * t[year] + error[hab,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/g_hab.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/g_hab.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/perhabitat_g_output.rds")
