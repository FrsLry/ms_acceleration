## Fitting the trend models
library(jagsUI)
library(dplyr)
library(tidyr)
library(tibble)
## Load the data
mean_N_habitats <- readRDS("mixed_model/save_samples/mean_N_habitats.rds")
sd_N_habitats <- readRDS("mixed_model/save_samples/sd_N_habitats.rds")

## JAGS data
data_list <- list(N = as.matrix(mean_N_habitats),
                  sd = as.matrix(sd_N_habitats),
                  habitats = nrow(mean_N_habitats),
                  years = ncol(mean_N_habitats),
                  t = 1:ncol(mean_N_habitats))

## JAGS model
cat(
  "model{

  # Priors
  for(i in 1:habitats){
    beta0[i] ~ dnorm(grand.beta0, grand.tau.beta0)
    beta1[i] ~ dnorm(grand.beta1, grand.tau.beta1)
  }

  sigma ~ dunif(0, 1000000)
  tau = 1/(sigma*sigma)

  grand.beta0 ~ dnorm(0, 1e-100)
  grand.sigma.beta0 ~ dunif(0, 1000000)
  grand.tau.beta0 = 1/(grand.sigma.beta0*grand.sigma.beta0)

  grand.beta1 ~ dnorm(0, 1e-100)
  grand.sigma.beta1 ~ dunif(0, 1000000)
  grand.tau.beta1 = 1/(grand.sigma.beta1*grand.sigma.beta1)

  # Likelihood
  for(hab in 1:habitats){
    for(year in 1:years){

      tau.obs[hab,year] <- 1/(sd[hab,year]*sd[hab,year])
      error[hab,year] ~ dnorm(0, tau.obs[hab,year])

      N[hab,year] ~ dnorm(mu[hab,year], tau)
      mu[hab,year] <- beta0[hab] + beta1[hab] * t[year] + error[hab,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/N_hab.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/N_hab.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/perhabitat_N_output.rds")
