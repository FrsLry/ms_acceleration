# This script is made to demonstrate how the dynamic N-mixture model (or Dail-Madsen, DM, 10.1111/j.1541-0420.2010.01465.x) works on simulated data.
# It is strongly inspired from AHM vol. 2 by Kéry and Royle (ISBN 978-0-12-823768-7). If you use this script, please give THEM credits.
# We simulate data and fit a DM model to reproduce our analysis as much as possible.
# Users can modify the data simulation in order to observe the consequences on fitting the model.

library(jagsUI)

# First, we set the characteristics of the dataset
# We warn that using the same number of sites (1033) and the same length of time series (35) as in our analysis
# will take a very long time to fit
sites = 100 # the number of sites, or routes of the BBS data
surveys = 1 # we used the DM model without repeated count, thus this variable should be set to 1
years = 10 # length of the time series

# Now, set values to the parameters we want to predict
lambda = 4 # abundance at time 1
phi = 0.8 # probability of survival
gamma = 1.5 # recruitment parameter
p = 0.7 # detection probability

### Simulator function
simDM0 <- function(nsites = sites, nsurveys = surveys, nyears = years, lambda = lambda,
                   phi = phi, gamma = gamma, p = p){

  ## No covariates (like in our analysis)
  y <- array(NA, dim = c(nsites, nyears, nsurveys))
  N <- matrix(NA, nsites, nyears)
  S <- R <- matrix(NA, nsites, nyears-1)
  N[,1] <- rpois(nsites, lambda) # Initial state
  for(t in 1:(nyears-1)) { # State dynamics
    S[,t] <- rbinom(nsites, N[,t], phi) # Survival process
    R[,t] <- rpois(nsites, gamma) # Recruitment process
    N[,t+1] <- S[,t] + R[,t]
  }
  for(j in 1:nsurveys){ # Observation process
    y[,,j] <- rbinom(nsites*nyears, N, p)
  }

  # Put observed data into two dimensions
  yy <- array(NA, dim = c(nsites, nsurveys*nyears))
  for(t in 1:nyears){
    yy[,(nsurveys * t-(nsurveys-1)):(nsurveys*t)] <- y[,t,]
  }
  return(list(nsites = nsites, nsurveys = nsurveys, nyears = nyears,
              lambda = lambda, phi = phi, gamma = gamma, p = p, N = N, S = S, R = R,
              y = y, yy = yy))
}

# Execute function
set.seed(2017, kind = "L'Ecuyer")
str(data <- simDM0(nsites = sites, nsurveys = surveys, nyears = years, lambda = lambda,
                   phi = phi, gamma = gamma, p = p))

# Bundle data set
str(bdata <- list(C = data$y, nsites = dim(data$y)[1], nsurveys = dim(data$y)[3],
                  nyears = dim(data$y)[2]))

# Specify model in BUGS language
cat(file = "simulations/DM1.txt","
model {
  # Priors
  lambda ~ dunif(0, 100)   # Initial site-specific abundance
  phi ~ dunif(0, 1)        # Apparent survival
  gamma ~ dunif(0, 5)      # Recruitment rate
  p ~ dunif(0, 1)          # Detection probability

  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t])     # Survival process
      R[i,t+1] ~ dpois(gamma)        # 'absolute' recruitment = 'constant'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,t,j] ~ dbin(p, N[i,t])
      }
    }
  }
}
")

# Initial values that we used in our analysis
Nst <- apply(data$y, c(1,2), max) + 2
Nst[, 2:dim(data$y)[2]] <- NA # Deterministic, N <- S + R
R1 <- apply(data$y, c(1,2), max) # Observed max. counts + 1 as inits
R1[,1] <- NA
inits <- function(){list( lambda = runif(1, 6, 16), phi = runif(1),
                          gamma = runif(1), p = runif(1), N = Nst, R = R1 + 1 )}

# Parameters monitored
params <- c("lambda", "phi", "gamma", "p") # in our analysis, we also monitored N, S, R

# MCMC settings
# In our analysis, the settings were:
na <- 1000 ; ni <- 100000 ; nt <- 10 ; nb <- 75000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "simulations/DM1.txt", n.adapt = na, n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
traceplot(out1, layout=c(2,3))
print(out1, 3)
