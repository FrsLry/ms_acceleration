## Fitting the trend models
library(jagsUI)
## Load the data
mean_gt0_overall <- readRDS("mixed_model/save_samples/mean_gt0_overall.rds")
sd_gt0_overall <- readRDS("mixed_model/save_samples/sd_gt0_overall.rds")

## For JAGS, put non zero values (but very small) of sd for time 1
sd_gt0_overall[,1] <- 1e-10

## JAGS data
data_list <- list(gt0 = as.matrix(mean_gt0_overall),
                  sd = as.matrix(sd_gt0_overall),
                  routes = nrow(mean_gt0_overall),
                  years = ncol(mean_gt0_overall),
                  t = 1:ncol(mean_gt0_overall))

## JAGS model
cat(
  "model{

  # Priors (uninformative)
  for(i in 1:routes){
    beta0[i] ~ dnorm(grand.beta0, grand.tau.beta0)
    beta1[i] ~ dnorm(grand.beta1, grand.tau.beta1)
  }

  sigma ~ dunif(0, 100)
  tau = 1 / (sigma * sigma)

  grand.beta0 ~ dnorm(0, 1e-10)
  grand.sigma.beta0 ~ dunif(0, 100)
  grand.tau.beta0 = 1 / (grand.sigma.beta0 * grand.sigma.beta0)

  grand.beta1 ~ dnorm(0, 1e-10)
  grand.sigma.beta1 ~ dunif(0, 100)
  grand.tau.beta1 = 1 / (grand.sigma.beta1 * grand.sigma.beta1)


  # Likelihood
  for(r in 1:routes){
    for(year in 1:years){

      tau.obs[r,year] <- 1/(sd[r,year]*sd[r,year])
      error[r,year] ~ dnorm(0, tau.obs[r,year])

      gt0[r,year] ~ dnorm(mu[r,year], tau)
      mu[r,year] = beta0[r] + beta1[r] * t[year] + error[r,year]

    }
  }

  }", fill = TRUE, file = "mixed_model/models/gt0_overall.txt")

## Hyperparams
params <- c("beta0", "beta1", "grand.beta0", "grand.beta1")
na <- 1000 ; ni <- 10000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

## Fit the model
model <- jags(model.file = "mixed_model/models/gt0_overall.txt",
              parameters.to.save = params,
              data = data_list,
              n.adapt = na,
              n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb)

model

saveRDS(model, "mixed_model/outputs/overall_gt0_output.rds")

