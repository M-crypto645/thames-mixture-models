rm(list=ls())
library(mclust)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# get data
data("banknote")
x <- banknote[, -1]
cl_target <- banknote$Status
n <- dim(x)[1]
p <- dim(x)[2]

# hyperparameters
G <- 3
alpha0 <- 1
nu0 <- p
beta0 <- 1#1e-5
m0 <- colMeans(x)
#W0 <- solve(var(x))/(2*nu0)
W0 <- var(x)*(2*nu0)

stan_data <- list(
  p = p, n = n, G = G, x = x,
  alpha0 = alpha0, nu0 = nu0, beta0 = beta0, m0 = m0, W0 = W0
)

# fit model in stan
fit <- stan(
  file = 'gmm_multi.stan', 
  data=stan_data, 
  iter=2000, warmup=1000, chains=4)
fit@sim
print(fit)

# posterior plots
sims <- extract(fit)
#save(sims,file='gmm_multi_stan_G3_Wishart.RData')
#load('gmm_multi_stan_G3_Wishart.RData')

save(sims,file='gmm_multi_stan_G3.RData')
load('gmm_multi_stan_G3.RData')

par(mfrow=c(2,3))
hist(sims$mu[,1,1])
hist(sims$mu[,1,2])
hist(sims$mu[,1,3])
hist(sims$mu[,1,4])
hist(sims$mu[,1,5])
hist(sims$mu[,1,6])

plot(sims$mu[,1,1])
plot(sims$mu[,1,2])
plot(sims$mu[,1,3])
plot(sims$mu[,1,4])
plot(sims$mu[,1,5])
plot(sims$mu[,1,6])

