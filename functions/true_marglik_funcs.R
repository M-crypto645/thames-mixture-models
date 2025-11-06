library(mvtnorm)

# simulates from the prior predictive distribution of the data
rdata = function(n,mustars,taus_tilde,sigma_tilde,seed=2024){
  set.seed(seed)
  mus = numeric(length(mustars))
  G = length(mustars)
  C = apply(rmultinom(n=n,size=1,prob=taus_tilde),2,function(x) which.max(x))
  y=numeric(n)
  for(g in 1:G){
    mus[g] = mustars[g] # rnorm(1,mean=mustars[g],1)
    y[C==g] = rnorm(sum(C==g),mean=mus[g],sd=sigma_tilde)
  }
  return(y)
}

# true marginal likelihood for Gaussian mixture models
# with known proportions and known variance, conditional on the clusters C
log_true_marglik_C = function(y,C_vec,taus_tilde,mustars,sigma_tilde){
  G = length(taus_tilde)
  check_feasible = TRUE
  #browser()
  # not included in the sum if not feasible
  for(g in which(taus_tilde==0)){
    check_feasible = check_feasible & (sum(C_vec==g)==0)
  }
  if(check_feasible){
    log_prior_given_c = sum(sapply((1:G)[taus_tilde!=0], function(g) log(taus_tilde[g])*sum(C_vec==g)))
  
    loglik_y_given_c_partial = numeric(length(y))
    mus = numeric(length(y))
    Sigma = diag(length(y))
    for(g in 1:G){
      mus[C_vec==g] = mustars[g]
      Sigma[C_vec==g,C_vec==g] = 1
      #loglik_y_given_c_partial[C_vec==g] = dnorm(y[C_vec==g],mustars[g],sqrt(sigma_tilde^2+1),log=TRUE)
    }
    #browser()
    loglik_y_given_c_partial = mvtnorm::dmvnorm(y, mean=mus, sigma=Sigma+sigma_tilde^2*diag(n), log = TRUE)
    loglik_y_given_c = sum(loglik_y_given_c_partial)
    #browser()
    return(log_prior_given_c + loglik_y_given_c)
  } else{
    return(-Inf)
  }
}

log_true_marglik = function(y,taus_tilde,mustars,sigma_tilde){
  
  G = length(taus_tilde)
  n <- length(y)
  
  # Generate all possible combinations using expand.grid
  set_values <- 1:G
  all_combinations <- expand.grid(rep(list(set_values), n))
  
  # Convert it to a matrix for easy access
  all_Cs <- as.matrix(all_combinations)
  
  # sum over all clusters using the log trick
  #browser()
  partial_log_marginal_liks  = apply(all_Cs,1,function(C_vec) log_true_marglik_C(y,C_vec,taus_tilde,mustars,sigma_tilde))
  true_marglik = log(sum(exp(partial_log_marginal_liks[!is.infinite(partial_log_marginal_liks)] - min(partial_log_marginal_liks[!is.infinite(partial_log_marginal_liks)])))) + 
    min(partial_log_marginal_liks[!is.infinite(partial_log_marginal_liks)])
  #browser()
  #log(rnorm())
  return(true_marglik)
}

rgibbs_mus = function(y,C_vec,mustars, sigma_tilde){
  G = length(mustars)
  abs_C_vec = sapply(1:G, function(g) sum(C_vec==g))
  sum_C_vec_Y = sapply(1:G, function(g) sum(y[C_vec==g]))
  sd = sqrt(1/(1+abs_C_vec/(sigma_tilde^2)))
  #browser()
  return(rnorm(G,mean=(sd^2)*(mustars + sum_C_vec_Y)/(sigma_tilde^2),sd=sd))
}

rgibbs_C_vec = function(y,mus,sigma_tilde,taus_tilde){
  G = length(taus_tilde)
  n = length(y)
  C_vec = numeric(n)
  #browser()
  for(i in 1:n){
    prob = taus_tilde*exp(-(mus-y[i])^2/(2*(sigma_tilde^2)))
    C_vec[i] = which(rmultinom(1,size=1,prob=prob)==1)
  }
  return(C_vec)
}

gibbs_sampling = function(iters, mus_init, y, sigma_tilde, taus_tilde, mustars, seed=2024){
  set.seed(seed)
  n = length(y)
  G = length(mustars)
  C_mat = matrix(nrow=iters,ncol=n)
  mus_mat = matrix(nrow=iters,ncol=G)
  #browser()
  mus_mat[1,] = mus_init
  for(i in 1:iters){
    C_mat[i,] = rgibbs_C_vec(y,mus_mat[i,],sigma_tilde,taus_tilde)
    if(i!=iters){
      mus_mat[i+1,] = rgibbs_mus(y,C_mat[i,],mustars, sigma_tilde)
    }
  }
  #browser()
  return(list(mus_mat=mus_mat,C_mat=C_mat))
}

# combines the sampling of y with the gibbs sampling of mu
y_mu_sampler = function(n,mus_tilde,sigma_tilde,taus_tilde,
                      mustars,sigmastar,taustars,
                      init,iters,
                      name,burn_in = 2000,seed=2024){
  
  # renaming stuff
  #browser()
  #mus_init = init
  data_taus_tilde = taustars
  prior_mustars = mus_tilde 
  
  G = length(mus_tilde)
  y = rdata(n=n,mustars=mustars,taus_tilde = data_taus_tilde, 
            sigma_tilde=sigmastar, seed=seed)
  #browser()
  mus_init = init(y)
  (log_true_marglik_G2 = log_true_marglik(y=y,taus_tilde=taus_tilde,mustars=prior_mustars,sigma_tilde=sigma_tilde))
  posterior_sample = gibbs_sampling(iters+burn_in, mus_init=mus_init, y, 
                                    sigma_tilde, taus_tilde,
                                    mustars=prior_mustars,seed=seed)
  posterior_sample$mus_mat = posterior_sample$mus_mat[-(1:burn_in),]
  posterior_sample$C_mat = posterior_sample$C_mat[-(1:burn_in),]
  
  return(list(results=array(c(posterior_sample$mus_mat),
                            dim=c(iters,G,1)),
              alloc_vec=posterior_sample$C_mat,
              name=name, y=y, truth=log_true_marglik_G2))
}

logposty = function(sims, y, sigma_tilde, taus_tilde, mus_tilde){
  #browser()
  iters = dim(sims)[1]
  G = dim(sims)[2]
  thetas = array(dim=c(dim(sims)[1:2],3))
  thetas[,,1] = sims[,,1]
  thetas[,,2] = matrix(sigma_tilde,nrow=iters,ncol=G)
  thetas[,,3] = matrix(sapply(1:G,function(g) rep(taus_tilde[g],iters)),ncol=G)
  
  iters = nrow(thetas)
  logliks = loglik_gmm(y,sims=thetas,G=G)
  logpriors = rowSums(sapply(1:G, function(g) dnorm(thetas[,g,1], mean = mus_tilde[g], log = TRUE)))
  return(logliks + logpriors)
}

# returns sample from the prior, evaluated at the loglikelihood
prior_sampler_marglik = function(y, mus_tilde, sigma_tilde, taus_tilde, iters, seed=2024){
  #browser()
  G = length(mus_tilde)
  prior_sample = sapply(1:G,function(g) rnorm(iters,mean=mus_tilde[g],1))
  mc_estim_partial = loglik_gmm(y,theta=cbind(prior_sample,matrix(sigma_tilde,nrow=iters,ncol=G),
                                              matrix(sapply(1:G,function(g) rep(taus_tilde[g],iters)),ncol=G)),G=G)
  return(mc_estim_partial)
}