## Fitting the trend models
library(jagsUI)
## Load the data
mean_g_families <- readRDS("mixed_model/save_samples/mean_g_families.rds")
sd_g_families <- readRDS("mixed_model/save_samples/sd_g_families.rds")

## JAGS data
data_list <- list(g = as.matrix(mean_g_families),
                  sd = as.matrix(sd_g_families),
                  families = nrow(mean_g_families),
                  years = ncol(mean_g_families),
                  t = 1:ncol(mean_g_families))

## JAGS model
cat(
  "model{

  # Priors
  for(i in 1:families){
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
  for(fam in 1:families){
    for(year in 1:years){

      tau.obs[fam,year] <- 1/(sd[fam,year]*sd[fam,year])
      error[fam,year] ~ dnorm(0, tau.obs[fam,year])

      g[fam,year] ~ dnorm(mu[fam,year], tau)
      mu[fam,year] <- beta0[fam] + beta1[fam] * t[year] + error[fam,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/g_fam.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/g_fam.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/perfamily_g_output.rds")
