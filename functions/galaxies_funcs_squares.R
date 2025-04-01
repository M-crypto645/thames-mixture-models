# same as galaxies_funcs, using variances

library(invgamma,LaplacesDemon)

# functions to calculate gaussian mixture model unnormalized log posterior
# (now vectorized)
loglik_gmm <- function(y,theta,G){
  #browser()
  mus <- theta[,1:G]
  sigma_squs <- theta[,(G+1):(2*G)]
  pis <- theta[,(2*G+1):(3*G)]
  
  log_single_y = Vectorize(function(x) 
    log(rowSums(sapply(1:G, 
                       function(g) pis[,g]*dnorm(x,mus[,g],sqrt(sigma_squs[,g]))))
    )
  )
  return(rowSums(log_single_y(y)))
  # return(sum(sapply(y,function(x){
  #   log(sum(pis*dnorm(x,mean=mus,sd=sqrt(sigma_squ))))
  # })))
}

#same as above, but without summing (now vectorized)
loglik_gmm_partial <- function(y,sims,G){
  browser()
  pis <- sims[,1:G,3]
  mus <- sims[,1:G,1]
  sigma_squs <- sims[,1:G,2]
  res = array(dim=c(nrow(pis),length(y),G))
  
  for(g in 1:G){
    log_single_y = Vectorize(function(x) 
      log(pis[,g]) + dnorm(x,mus[,g],sqrt(sigma_squs[,g]),log=TRUE))
    res[,,g] = log_single_y(y)
  }
  return(res)
  #return(sapply(1:G,function(g) log(pis[g])+dnorm(y,mean=mus[g],sd=sqrt(sigma_squs[g]),log=TRUE)))
}

#same as above, but without summing (now vectorized)
loglik_gmm_partial_transform <- function(y,sims,G){
  #browser()
  pis <- sims[,1:G,3]
  mus <- sims[,1:G,1]
  sigma_squs <- exp(sims[,1:G,2])
  res = array(dim=c(nrow(pis),length(y),G))
  
  for(g in 1:G){
    log_single_y = Vectorize(function(x) 
      log(pis[,g]) + dnorm(x,mus[,g],sqrt(sigma_squs[,g]),log=TRUE))
    res[,,g] = log_single_y(y)
  }
  return(res)
  #return(sapply(1:G,function(g) log(pis[g])+dnorm(y,mean=mus[g],sd=sqrt(sigma_squs[g]),log=TRUE)))
}

# loglik_gmm(y,theta,G)

# lp_gmm <- function(y,theta,G,m,R){
#   loglik_gmm(y,theta,G)+logprior_gmm(theta,G,m,R)
# }

# lp_gmm(y,theta,G,m,R)

lp_gmm_marginal <- function(y,theta,G,m,R){
  #browser()
  
  mus <- theta[,1:G]
  sigma_squs <- theta[,(G+1):(2*G)]
  pis <- theta[,(2*G+1):(3*G)]
  
  # set to 0 outside of support
  if(G>2){
    mask = (((pis > 0) & (rowSums(pis[,1:(G-1)])<=1)) & (sigma_squs>0))
  }else{
    mask = (((pis > 0) & (pis[,1]<=1)) & (sigma_squs>0))
  }
  
  l_total = loglik_gmm(y,theta,G)+logprior_gmm_marginal(theta,G,m,R)
  l_total[exp(rowSums(log(mask)))==0] = -Inf
  return(l_total)
}

lp_gmm_marginal_transform <- function(y,theta,G,m,R){
  #browser()
  
  mus <- theta[,1:G]
  
  # apply exp transform
  theta[,(G+1):(2*G)] = exp(theta[,(G+1):(2*G)])
  
  sigma_squs <- theta[,(G+1):(2*G)]
  pis <- theta[,(2*G+1):(3*G)]
  
  # set to 0 outside of support
  if(G>2){
    mask = (((pis > 0) & (rowSums(pis[,1:(G-1)])<=1)) & (sigma_squs>0))
  }else{
    mask = (((pis > 0) & (pis[,1]<=1)) & (sigma_squs>0))
  }
  
  # adjust for log tranform
  jacobian = rowSums(log(sigma_squs))
    
  l_total = loglik_gmm(y,theta,G)+logprior_gmm_marginal(theta,G,m,R)+jacobian
  l_total[exp(rowSums(log(mask)))==0] = -Inf
  # browser()
  return(l_total)
}

# lp_gmm_marginal(y,theta,G,m,R)

# logprior_gmm <- function(theta,G,m,R){
#   browser()
#   pis <- theta[1:G]
#   mus <- theta[(G+1):(2*G)]
#   sigmas <- theta[(2*G+1):(3*G)]  
#   C0 <- theta[length(theta)]
#   
#   lp <- dgamma(x = C0, shape = 0.2, rate = 10/(R^2), log = TRUE) +
#     sum(LaplacesDemon::dinvgamma(x = sigmas^2, shape = 2, scale = C0, log = TRUE)) +
#     # sum(invgamma::dinvgamma(x = sigmas^2, shape = 2, rate = C0, log = TRUE)) +
#     sum(dnorm(x = mus, mean = m, sd = R, log = TRUE)) +
#     LaplacesDemon::ddirichlet(x = pis, alpha = rep(1,G), log=TRUE)
#   return(lp)
# }

# log of the prior of the marginal distribution of pi,mu,sigma_squ (now vectorized)
logprior_gmm_marginal <- function(theta,G,m,R) {
  mus <- theta[,1:G]
  sigma_squs <- theta[,(G+1):(2*G)]
  pis <- theta[,(2*G+1):(3*G)]
  
  l_mus <- rowSums(sapply(1:G, function(g) dnorm(mus[,g], mean = m, sd = R, log = TRUE)))
  l_pis <- LaplacesDemon::ddirichlet(1:G/G, rep(1,G),log=TRUE) # constant wrt pis
  l_sigma_squs <- lgamma(2*G+0.2) - lgamma(0.2) +
    0.2*log(10/R^2) - (2*G+0.2) * log(rowSums(sigma_squs^(-1))+10/R^2) - 3*rowSums(log(sigma_squs))
  return(l_mus + l_pis + l_sigma_squs)
  
  
  
  # l_mus <- sum(function(x) dnorm(mus, mean = m, sd = R, log = TRUE))
  # l_pis <- LaplacesDemon::ddirichlet(pis, rep(1,G),log=TRUE)
  # l_sigma_squs <- lgamma(2*G+0.2) - lgamma(0.2) +
  #   0.2*log(10/R^2) - (2*G+0.2)*log(sum(sigma_squs^(-1))+10/R^2) - 3*sum(log(sigma_squs))
  #log_det_jac = G*log(2)+sum(log(sigmas)) (if we use square roots instead)
}

#' @title Posterior Samples from JAGS
#'
#' @description Simulates samples from the posterior for the galaxies dataset 
#'              using the priors from Richardson & Green (1997)
#'              Notes:  The output gives standard deviations, not variances, 
#'                      for sigma
#' @param n     [int>0]         size of the data
#' @param x     [vector]        the (univariate) data
#' @param R2,m  [constants]     hyperparameters, see Richardson & Green (1997)
#' @param k1    [int>1]         number of groups
#' @param iters [vector]        2*iters=number of simulations from the posterior
#' @param init  [list, optional]initial value          
#' @param seed  [int, optional] seed for the sampler                   
#' @return  results   [array]       tensor with dimensions (iters,G,3) standing
#'                                  for "iterations, group, (mu,sigma,pi)" resp.
#'          alloc_vec [matrix]      (2*iters) x n matrix of allocation vectors
#'          name      [string]      always set to "JAGS"
#'          truth     [int=0]       always set to 0 (there is not truth)
sim_jags = function(n, x, R2, m, k1, iters, init=NULL, seed=NULL){
  # simulate from the posterior using JAGS (instead of STAN)
  # takes as input the parameters of the Gaussian mixture model
  # iters:  number of simulations from the posterior. 
  #         1/2 of those samples, iters, is used as input for the THAMES (split)
  
  if(is.null(seed)){
    seed = 2024
  }
  
  if(!is.null(init)){
    browser()
  }
  model <- BMMmodel(x, k = k1,  priors = list(kind = "independence",
                                              hierarchical = "tau"),
                    initialValues = list(S0 = 2))
  
  #TODO: Find a good intitialization!
  #model$inits$mu=(1:k1)*100#print(round(model$inits$mu/100))
  #model$inits$eta
  #model$inits$tau=rep(100,k1)
  #model$inits$S0
  # control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
  #                        burn.in = 2000, n.iter = 2*iters)
  control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                         burn.in = 2000, n.iter = 2*iters, seed=seed)
  ## add the same prior than Model Selection for Mixture Models-Perspectives and
  #Strategies
  # Gilles Celeux, Sylvia Frühwirth-Schnatter, Christian Robert
  model$data$B0inv <-1/R2
  model$data$b0 <- m
  model$data$nu0Half <- 2
  model$data$g0Half <- 0.2
  model$data$g0G0Half <- 20 / R2

  z <- JAGSrun(x, model = model, control = control)
  #browser()
  results=array(c(z$results[,c((n+k1+1):(n+3*k1),(n+1):(n+k1))]),
                dim=c(2*iters,k1,3))
  
  # apply log transform
  results[,,2] = log(results[,,2])
  return(list(results=results,
              alloc_vec=z$results[,1:n],
              name="JAGS",y=x,truth=0))
}

#' @title Posterior Samples from STAN
#'
#' @description Simulates samples from the posterior for the galaxies dataset 
#'              using the priors from Richardson & Green (1997)
#'              Notes:  STAN provides NO allocation variables
#'                      The output gives standard deviations, not variances, 
#'                      for sigma
#' @param n     [int>0]         size of the data
#' @param y     [vector]        the (univariate) data
#' @param R,m   [constants]     hyperparameters, see Richardson & Green (1997)
#' @param G     [int>1]         number of groups
#' @param iters [vector]        2*iters=number of simulations from the posterior
#' @param init  [list, optional]initial value
#' @param seed  [int, optional] seed for the sampler                               
#' @return  results   [array]       tensor with dimensions (iters,G,3) standing
#'                                  for "iterations, group, (mu,sigma,pi)" resp.
#'          alloc_vec [NULL]        always set to NULL (there are none for STAN)
#'          name      [string]      always set to "STAN"
#'          truth     [int=0]       always set to 0 (there is no truth)
sim_stan = function(n,y,R,m,G,iters, init=NULL, seed=NULL){
  # simulate from the posterior using JAGS (instead of STAN)
  # takes as input the parameters of the Gaussian mixture model
  # iters:  2*iters is the number of simulations from the posterior. 
  #         1/2 of those samples, iters, is used as input for the THAMES (split)
  
  #browser()
  if(is.null(seed)){
    seed=2024
  }
  set.seed(seed)
  
  if(is.null(init)){
    init = runif(1,-2,2)
  }
  
  stan_data <- list(
    n = length(y),
    y = y,
    R = R,
    m = m,
    G = G
  )
  iters=2*iters
  # stan model
  stancode <- 'data {
    int<lower=1> n;
    vector[n] y;
    real<lower=0> R;
    real m;
    int<lower=1> G;
  }
  parameters {
    vector[G] mu; // component means
    vector<lower=0>[G] sigma_squared; // component variances
    simplex[G] pi; // mixture weights
    real<lower=0> C0; // hyperparameter
  }
  transformed parameters {
    vector[G] log_pi = log(pi);  // cache log calculation
    real lp; // log prior
    matrix[n,G] lls; // log-likelihood terms
    for(i in 1:n){
      for(g in 1:G){
        lls[i,g] = log_pi[g] + normal_lpdf(y[i] | mu[g], sqrt(sigma_squared[g]));
      }
    }
    lp = gamma_lpdf(C0 | 0.2, 10/(R^2));
    lp += inv_gamma_lpdf(sigma_squared | 2, C0);
    lp += normal_lpdf(mu | m, R);
    lp += dirichlet_lpdf(pi | rep_vector(1,G));
  }
  model {
    // priors
    target += gamma_lpdf(C0 | 0.2, 10/(R^2));
    target += inv_gamma_lpdf(sigma_squared | 2, C0);
    target += normal_lpdf(mu | m, R);
    target += dirichlet_lpdf(pi | rep_vector(1,G));
    
    // (log) likelihood
    //vector[G] log_pi = log(pi);  // cache log calculation
    //for (i in 1:n) {
      //vector[G] lls = log_pi;
      //for (g in 1:G) {
        //lls[g] += normal_lpdf(y[i] | mu[g], sqrt(sigma_squared[g]));
      //}
      //target += log_sum_exp(lps);
    //}
    for(i in 1:n){
      target += log_sum_exp(lls[i,]);
    }
  } 
  '
  
  # compile model
  stanmodel <- rstan::stan_model(model_code = stancode)
  
  # fit models
  stanfit <- rstan::sampling(stanmodel, data = stan_data, iter = iters/2, seed=seed, init=init)
  
  set.seed(seed)
  sims <- rstan::extract(stanfit,permute=FALSE)
  sims = sims[,sample(1:4,4),]
  reorder_sims = matrix(nrow=iters,ncol=dim(sims)[3])
  for(i in 0:(dim(sims)[2]-1)){
    reorder_sims[1+((iters*i/dim(sims)[2]):(iters*(i+1)/dim(sims)[2]-1)),] = sims[,i+1,]
  }
  sims = reorder_sims
  #browser()
  
  results=array(c(sims[,1:(3*G)]),dim=c(iters,G,3))
  
  # apply log transform
  results[,,2] = log(results[,,2])
  # format for label.switching
  list(results=results,
       alloc_vec=NULL,
       name="STAN",y=y,truth=0)
}
