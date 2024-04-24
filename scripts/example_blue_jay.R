# This code is meant to run the dynamic N-mixture model for the Blue Jay.
# The jagsData contains the input data for running the model for 564 species.
# Please note that all covariates were standardized and that some covariates are not used in the models due to convergence issues (namely land cover and elevation).
# Users can select the species they want by modifying the `species` argument.
rm(list=ls())
library(jagsUI)
library(MCMCvis)
jagsData<-readRDS('data/jags_data_par.rds')
inits<-readRDS('data/inits.rds')

species <- "Blue Jay" # change the species name accordingly (species names can be found in the jagsData list)
jagsData <- jagsData[species]
inits <- inits[species] # created following Kéry & Royle (2020)

# Parameters monitored
params <- c('beta.time',
            'beta.temp',
            'beta.sky',
            'beta.wind',
            'mean.phi',
            'mean.gamma',
            'mean.p',
            'mean.lam',
            'S', 'R',
            'N')

# MCMC settings
na <- 1000 ; ni <- 100000 ; nt <- 10 ; nb <- 75000 ; nc <- 3

# Fit the model
# please note that the maximum estimated time is ca. 2 days, as it took this time for running the 564 models on the Ohio Super Computer
# please note that the function jags.basic() was used in order to not overwhelm the supercomputer node RAM
set.seed(333)
jagsOut<-jagsUI::jags.basic(jagsData[[1]], inits[[1]], params, 'model/jagsModel.txt', n.adapt = na, n.chains = nc,
                            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
saveRDS(jagsOut, paste('model_output/', species, '.rds', sep=''))

# Now summarize the output of the jags.basic() function
output <- readRDS(paste('model_output/', species, '.rds', sep=''))
sum <- MCMCsummary(output,
                   params = c('beta.time',
                              'beta.temp',
                              'beta.sky',
                              'beta.wind',
                              'mean.phi',
                              'mean.gamma',
                              'mean.p',
                              'mean.lam',
                              'S', 'R', 'N'))
saveRDS(sum, file = paste0('summary_models/summary_', species, '.rds', sep=''))
