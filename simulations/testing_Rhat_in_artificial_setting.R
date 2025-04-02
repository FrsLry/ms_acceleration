
library(coda)


# POISSON CONVERGENCE

# parameters of the initial gamma distribution of lambda parameter of poisson
shape <- 0.005
rate <- 8

# plot the distribution
hist(rgamma(10000, shape = shape, rate = rate), breaks= 100)

# chain length
N <- 2000

# sample lambda parameters from an exponentiated normal (i.e., log link)
lambda.chain1 <- as.mcmc(rgamma(N, shape = shape, rate = rate))
lambda.chain2 <- as.mcmc(rgamma(N, shape = shape, rate = rate))
lambda.chains <- as.mcmc.list(list(lambda.chain1, lambda.chain2))

# sample counts from poisson with the lambdas
count.chain1 <- as.mcmc(rpois(N, lambda = lambda.chain1))
count.chain2 <- as.mcmc(rpois(N, lambda = lambda.chain2))
count.chains <- as.mcmc.list(list(count.chain1, count.chain2))

# Rhat
Rhat.lambda<- gelman.diag(lambda.chains)$psrf[1]
Rhat.count <- gelman.diag(count.chains)$psrf[1]

# check convergence
par(mfrow=c(2,2))
plot(lambda.chains, main = paste("Rhat =", round(Rhat.lambda,3)), auto.layout = FALSE)
plot(count.chains, main = paste("Rhat =", round(Rhat.count,3)), auto.layout = FALSE)

# ------------------------------------------------------------------------------

# BINOMIAL CONVERGENCE

# parameters of the initial beta distribution
a <- 9
b <- 5
size <- 10 # think about this as the true abundance at a site i time t

# chain length
N <- 1000

# plot the distribution
hist(rbeta(10000, a, b), breaks=100, xlim = c(0,1))


# sample P parameters from beta distribution
P.chain1 <- as.mcmc(rbeta(N, a, b))
P.chain2 <- as.mcmc(rbeta(N, a, b))
P.chains <- as.mcmc.list(list(P.chain1, P.chain2))

# sample counts from binomial with the Ps
count.chain1 <- as.mcmc(rbinom(N, size = size, prob = P.chain1))
count.chain2 <- as.mcmc(rbinom(N, size = size, prob = P.chain2))
count.chains <- as.mcmc.list(list(count.chain1, count.chain2))

# Rhat
Rhat.P <- gelman.diag(P.chains)$psrf[1]
Rhat.count <- gelman.diag(count.chains)$psrf[1]

# check convergence
par(mfrow=c(2,2))
plot(P.chains, main = paste("Rhat =", round(Rhat.P,3)), auto.layout = FALSE)
plot(count.chains, main = paste("Rhat =", round(Rhat.count,3)), auto.layout = FALSE)








