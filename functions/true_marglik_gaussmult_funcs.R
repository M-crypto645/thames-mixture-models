#source("functions/thames_gmm_multi.R")

### VII model ###

calc_lambda_vii = function(y,nu=2,G=3){
  set.seed(2025)
  res=Mclust(y, G=G,verbose=FALSE)
  #browser()
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  varymclust = matrix(rowSums(simplify2array(mclapply(1:G,function(g) c(res$parameters$variance$sigma[,,g])))/G),
                      nrow=nrow(res$parameters$variance$sigma[,,1]),
                      ncol=ncol(res$parameters$variance$sigma[,,1]))
  return(diag(diag(varymclust))*nu)
}

# combines the sampling of y with the gibbs sampling of mu
y_theta_sampler_gaussmulti_vii = function(n, nu, alpha_0,
                                      mustars,sigmastars,taustars,
                                      init,iters,burn_in = 2000,seed=2024,y=NULL){
  #browser()
  G = length(alpha_0)
  
  if(is.null(y)){
    y = rdata_gaussmulti(n,mustars,taustars,sigmastars,seed=seed)    
  }
  
  # choose beta and lambda via the data
  beta = calc_beta(y)
  lambda = calc_lambda_vii(y,nu,G)
  kappa_0 = calc_kappa_0(y)
  
  R = ncol(y)
  #browser()
  
  # initialize via Mclust
  Cs_init = Mclust(y, G=G)$classification
  # res=Mclust(y, G=G)
  # res$parameters
  
  if(is.null(y)){
    (log_marglik = log_true_marglik_gauss_multi_vii(y, nu, lambda, alpha_0, kappa_0, beta, mustars))
  } else{
    log_marglik = 0 # presume truth is unknown if y is known
  }
  
  tic()
  posterior_sample = gibbs_sampling_gaussmulti_vii(iters+burn_in, Cs_init=Cs_init, y, 
                                               alpha_0, kappa_0, nu, lambda, beta ,seed=seed)
  toc()
  posterior_sample$theta = posterior_sample$theta[-(1:burn_in),,]
  posterior_sample$C_mat = posterior_sample$C_mat[-(1:burn_in),]
  
  return(list(results=posterior_sample$theta,
              alloc_vec=posterior_sample$C_mat,
              name="GIBBS", y=y, truth=log_marglik))
}

# the gibbs sampler
gibbs_sampling_gaussmulti_vii = function(iters, Cs_init=Cs_init, y,
                                     alpha_0, kappa_0, nu, lambda, beta ,seed=seed){
  set.seed(seed)
  #browser()
  G = length(alpha_0)
  R = ncol(y)
  n = nrow(y)
  #browser()
  C_mat = matrix(nrow=iters,ncol=n)
  theta = array(dim=c(iters,G,1+2*R))
  
  C_mat[1,] = Cs_init
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  #browser()
  
  try(load("data/C_mat.RData"))
  # if(iters==11000){ # TODO REMOVE
  #   load("data/theta.RData")
  #   C_mat = C_mat[1:iters,]
  #   theta = theta[1:iters,,]
  # } else{
    # reload if possible (to save time)
  if((max(C_mat[!is.na(C_mat)])!=G)|(ncol(C_mat)!=n)|(nrow(C_mat)!=iters)|
     (which.max(rowSums(is.na(C_mat)))==2)|(mean(C_mat[1,] == Cs_init)<1)){
    C_mat = matrix(nrow=iters,ncol=n)
    C_mat[1,] = Cs_init
  } else{
    load("data/theta.RData")
  }
  #}
  
  # no need to sample if already done
  if(max(rowSums(is.na(C_mat)))>0){
    tic()
    for(i in which.max(rowSums(is.na(C_mat))):iters){
      #tic()
      # if(i==21){
      #browser()
      # }
      #print(i)
      if(round(i/100)==i/100){
        save(C_mat,file="data/C_mat.RData")
        save(theta,file="data/theta.RData")
        print(i)
      }
      C_vec = C_mat[i-1,]
      theta_i = rgibbs_theta_gaussmulti_vii(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=((seed*iters):((seed+1)*(iters)))[i-1])
      #theta_i = rgibbs_theta_gaussmulti(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=((seed*1000):((seed+1)*(1000)))[i-1])
      #browser()
      theta[i-1,,] = t(simplify2array(mclapply(1:G, function(g) c(theta_i[[g]]$mu_g,log(diag(theta_i[[g]]$Sigma_g)),theta_i[[g]]$pi_g),mc.cores = num_use_cores)))
      
      C_vec_i = rgibbs_C_vec_gaussmulti(y,theta_i,sigma_tilde,taus_tilde,seed=i)
      C_mat[i,] = C_vec_i
      #toc()
      gc()
    }
    #browser()
    C_vec = C_mat[iters,]
    theta_i = rgibbs_theta_gaussmulti_vii(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=iters)
    
    theta[iters,,] = t(simplify2array(mclapply(1:G, function(g) c(theta_i[[g]]$mu_g,log(diag(theta_i[[g]]$Sigma_g)),theta_i[[g]]$pi_g),mc.cores = num_use_cores)))
    
    save(C_mat,file="data/C_mat.RData")
    save(theta,file="data/theta.RData")
    #0.21*iters/60/60
    toc()
  }
  
  #browser()
  
  posterior_sample = list(theta=theta,C_mat = C_mat)
  
  return(posterior_sample)
}

# sampling conditioned on theta
rgibbs_theta_gaussmulti_vii = function(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed){
  
  G=length(alpha_0)
  alpha = alpha_0 + sapply(1:G, function(g) sum(C_vec==g))
  pis = rdirichlet(1, alpha) 
  #browser()
  rgibbs_theta_g = function(g){
    set.seed(seed+g)
    # define setting for this particular class
    n_g = sum(C_vec==g)
    y_g = y[C_vec==g,]
    beta_g = beta#s[,g]
    lambda_g = lambda
    nu_g = nu
    alpha_g = alpha_0
    
    # updated parameters on the conjugate prior (see Gelman et al)
    if(n_g==0){
      # set to prior if there is no information
      nu_ng = nu_g
      lambda_ng = lambda_g
      beta_ng = beta_g
      kappa_ng = kappa_0
    } else if (n_g==1){
      ybar_g = y_g
      S_ng = t(t(y_g))%*% t(y_g)
      
      kappa_ng = kappa_0 + n_g
      nu_ng = 2*nu_g + n_g
      nu_ng = nu_ng/2
      beta_ng = (1/(1+n_g/kappa_0))*beta_g+ybar_g*n_g/(kappa_0+n_g)
      lambda_ng = lambda_g*2 + S_ng + (t(t(ybar_g))-t(t(beta_g)))%*%t(t(t(ybar_g))-t(t(beta_g))) *(kappa_0*n_g)/(kappa_0+n_g)
      lambda_ng = lambda_ng/2
    } else{
      ybar_g = colMeans(y_g)
      S_ng = var(y_g)*(n_g-1)
      beta_ng = (1/(1+n_g/kappa_0))*beta_g+ybar_g*n_g/(kappa_0+n_g)
      kappa_ng = kappa_0 + n_g
      nu_ng = 2*nu_g + n_g
      nu_ng = nu_ng/2
      lambda_ng = 2*lambda_g + S_ng + (t(t(ybar_g))-t(t(beta_g)))%*%t(t(t(ybar_g))-t(t(beta_g))) *(kappa_0*n_g)/(kappa_0+n_g)
      lambda_ng = lambda_ng/2
    }
    
    # draw from the posterior given Z
    browser()
    Sigma_g = diag(rinvgamma(R,nu_ng, diag(lambda_ng)))
    mu_g = mvtnorm::rmvnorm(1,mean=beta_ng,
                            sigma=Sigma_g/kappa_ng)
    return(list(pi_g=pis[g],mu_g=mu_g,Sigma_g=Sigma_g))
  }
  #browser()
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  theta_i = mclapply(1:G,function(g) rgibbs_theta_g(g),mc.cores = num_use_cores)
  return(theta_i)
  
}


# calculates the true marginal likelihood given the data
log_true_marglik_gauss_multi_vii = function(y, nu, lambda, alpha_0, kappa_0, beta, mustars){
  
  R = ncol(lambda)
  G = length(alpha_0)
  C_vec = G-rowSums(sapply(1:G, function(g) rowSums(y) <= (R*mustars[1,]+100)[g]))+1
  #C_vec = G-(y[,1] <= beta) 
  # we can do this because the clusters are already extremely well seperated
  
  log_marglik = 0
  for(g in 1:G){
    #browser()
    # define setting for this particular class
    n_g = sum(C_vec==g)
    y_g = y[C_vec==g,]
    beta_g = beta#s[,g]
    lambda_g = lambda
    nu_g = nu
    alpha_g = alpha_0
    
    # updated parameters on the conjugate prior (see Gelman et al)
    ybar_g = colMeans(y_g)
    beta_ng = (1/(1+n_g/kappa_0))*beta_g+ybar_g*n_g/(kappa_0+n_g)
    kappa_ng = kappa_0 + n_g
    nu_ng = 2*nu_g + n_g
    
    S_ng = var(y_g)*(n_g-1)
    lambda_ng = 2*lambda_g + S_ng + (t(t(ybar_g))-t(t(beta_g)))%*%t(t(t(ybar_g))-t(t(beta_g))) *(kappa_0*n_g)/(kappa_0+n_g)
    lambda_ng = lambda_ng/2
    nu_ng = nu_ng/2
    #browser()
    # compute marginal likelihood via Bayes' rule
    mu_g = mustars[,g]
    sigma_g = diag(R) # result should not change with sigma_g,mu_g
    log_prior_g = sum(dinvgamma(diag(sigma_g),nu_g,diag(lambda_g),log=TRUE)) + mvtnorm::dmvnorm(mu_g,beta_g,sigma_g/kappa_0,log=TRUE)
    #log_prior_g = dinvwishart(sigma_g,nu_g,lambda_g,log=TRUE) + mvtnorm::dmvnorm(mu_g,beta_g,sigma_g/kappa_0,log=TRUE)
    loglik_g = sum(mvtnorm::dmvnorm(y_g,mu_g,sigma_g,log=TRUE))  
    logpost_g = sum(dinvgamma(diag(sigma_g),nu_ng,diag(lambda_ng),log=TRUE)) + mvtnorm::dmvnorm(mu_g,beta_ng,sigma_g/kappa_ng,log=TRUE)
    (log_marglik_g = log_prior_g + loglik_g - logpost_g)
    #browser()
    # (log_marglik_g = (nu_g/2)*sum(log(eigen(lambda_g)$values))-(nu_g*R/2)*log(2)-lmvgamma(nu_g/2,R)+
    #     (nu_ng*R/2)*log(2)+lmvgamma(nu_ng/2,R)+(-nu_ng/2)*sum(log(eigen(lambda_ng)$values))-
    #     (n_g*R)*log(sqrt(2*pi))-log(kappa_ng)*R/2+log(kappa_0)*R/2)
    
    log_marglik = log_marglik + log_marglik_g
  }
  
  # log_prior_given_c = sum(sapply((1:G)[taus_g!=0], function(g) log(taus_g[g])*sum(C_vec==g)))
  ns = sapply(1:G,function(g) sum(C_vec==g))
  lprior_z = lgamma(sum(alpha_0)) - lgamma(sum(alpha_0)+n) + sum(lgamma(ns+alpha_0)-lgamma(alpha_0))
  
  log_marglik = log_marglik +lprior_z + lfactorial(G) # only if there is a symmetric prior
  return(log_marglik)
}

logprior_gmm_gaussmulti_vii <- function(theta, G, nu, lambda, beta, alpha_0, kappa_0) {
  
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  logpriors = mclapply(1:iters,
                       function(i) logprior_gmm_i_gaussmulti_vii(theta[i,], G, nu, lambda, beta, alpha_0, kappa_0),
                       mc.cores = num_use_cores)
  return(simplify2array(logpriors))
  
}

logprior_gmm_i_gaussmulti_vii <- function(theta, G, nu, lambda, beta, alpha_0, kappa_0) {
  browser()
  R = length(beta)
  
  param_list = transform_to_params_vii(theta,G,R)
  
  pis = param_list$pis
  mus = param_list$mus
  sigmas = param_list$sigmas
  
  log_pis = log((sum(pis[-1])<=1)&(mean(pis>0)==1))
  if(!is.infinite(log_pis)){
    log_pis = log((sum(pis[-1])<=1)&mean(pis>0))+ddirichlet(x=pis,alpha=alpha_0,log=TRUE)
  }
  
  # sigmas are added due to the chain rule of the log transform
  log_sigmas = try(sum(sapply(1:G, function(g) log(diag(sigmas[[g]]))+dinvgamma(diag(sigmas[[g]]),nu,diag(lambda),log=TRUE))))
  #log_sigmas = try(sum(dInvWishart(simplify2array(sigmas),df=nu,Sigma=lambda,log=TRUE)))
  if(is.character(log_sigmas)|is.na(log_sigmas)){
    log_sigmas = -Inf
  }
  
  log_mus = try(sum(sapply(1:G, function(g) mvtnorm::dmvnorm(mus[g,],mean=beta,sigmas[[g]]/kappa_0,log=TRUE))))
  if(is.character(log_mus)){
    log_mus = -Inf
  }
  
  return(log_pis+log_sigmas+log_mus)
  
}

p0hat_gaussmulti_vii <- function(y, theta, G) {
  
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  p0hats = simplify2array(mclapply(1:iters,
                                   function(i) p0hat_i_vii(y, theta[i,], G),
                                   mc.cores = num_use_cores))
  #browser()
  return(mean(p0hats))
  
}

p0hat_i_vii  <- function(y, theta, G) {
  browser()
  R = dim(y)[2]
  param_list = transform_to_params_vii(theta,G,R)
  
  sigmas = param_list$sigmas
  browser()
  postprob_mat = simplify2array(lapply(1:G,function(g) log(param_list$pis[g]) + 
                                         mvtnorm::dmvnorm(y, param_list$mus[g,], 
                                                          sigmas[[g]],log=TRUE)))
  maxlogrows = do.call(pmax, c(as.data.frame(postprob_mat)))
  postprob_mat_normalized = exp(postprob_mat - maxlogrows) / rowSums(exp(postprob_mat - maxlogrows))
  p0estim = mean(sapply(1:G, function(g) prod(1-postprob_mat_normalized[,g])))
  return(p0estim)
}

loglik_gmm_gaussmulti_vii <- function(y, theta, G){
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  logliks = mclapply(1:iters,
                     function(i) loglik_gmm_i_gaussmulti_vii(y, theta[i,], G),
                     mc.cores = num_use_cores)
  return(simplify2array(logliks ))
}

transform_to_matrix_vii = function(uppertri,R){
  A = matrix(0, R, R)
  A[upper.tri(A, diag = TRUE)] = uppertri
  A = A + t(A)
  diag(A) = diag(A) / 2
  return(A)
}

transform_to_params_vii = function(theta,G,R){
  
  index = 0
  mus = t(matrix(theta[(index+1):(index+R*G)], R, G))
  index = index+R*G
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  
  sigmas = mclapply(1:G, 
                    function(g) diag(exp(theta[(index+1+(g-1)*R):(index+g*R)])),
                    mc.cores = num_use_cores )
  
  index = index + R*G
  
  pis = theta[(length(theta)-G+1):length(theta)]
  
  return(list(mus=mus,sigmas=sigmas,pis=pis))
}

# ASSUMING SIGMA IS THE ACTUAL COVARIANCE MATRIX (NOT ITS INVERSE)
loglik_gmm_i_gaussmulti_vii <- function(y, theta, G) {
  
  #print("likelihood function: start")
  browser()
  n = dim(y)[1]
  R = dim(y)[2]
  
  param_list = transform_to_params_vii(theta,G,R)
  
  pis = param_list$pis
  mus = param_list$mus
  sigmas = param_list$sigmas
  
  isinfinite = !((sum(pis[-1])<=1)&(mean(pis>0)==1)&(!is.character(try(dInvWishart(simplify2array(sigmas),df=R,Sigma=diag(R),log=TRUE)))))
  if(isinfinite){
    log_sum_lik = -Inf # check for invalid values (can happen with MC)
  } else{
    log_sum_lik = mvnfast::dmixn(X=y,mus,sigmas,pis,log=TRUE)
  }
  return(sum(log_sum_lik))
}


logposty_gaussmulti_vii = function(thetas, G, y, nu, alpha_0, kappa_0){
  #browser()
  beta = calc_beta(y)
  lambda = calc_lambda_vii(y,nu,G)
  kappa_0 = calc_kappa_0(y)
  
  iters = nrow(thetas)
  
  logliks  = loglik_gmm_gaussmulti_vii(y,thetas,G)
  logpriors = logprior_gmm_gaussmulti_vii(thetas, G, nu, lambda, beta, alpha_0, kappa_0) 
  return(logliks + logpriors)
}


### VII model ###


### VVV model (no transform) ###

calc_lambda = function(y,nu,G) {
  #browser()
  set.seed(2025)
  res=Mclust(y, G=G,verbose=FALSE)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  varymclust = matrix(rowSums(simplify2array(mclapply(1:G,function(g) c(res$parameters$variance$sigma[,,g])))/G),
                      nrow=nrow(res$parameters$variance$sigma[,,1]),
                      ncol=ncol(res$parameters$variance$sigma[,,1]))
  #diag(var(y)/(G^(R/2)))
  #browser()
  return(varymclust*(nu+ncol(y)+1))
}

calc_kappa_0 = function(y){
  kappa_0 = 0.00001#.001
  #kappa_0=1/diff(range(y))^2
  return(kappa_0)
}

# calc_lambda = function(y,nu) {
#   return(diag(ncol(y))*(nu+ncol(y)+1))
# }

# combines the sampling of y with the gibbs sampling of mu
y_theta_sampler_gaussmulti = function(n, nu, alpha_0,
                                      mustars,sigmastars,taustars,
                                      init,iters,burn_in = 2000,seed=2024){
  #browser()
  G = length(alpha_0)
  y = rdata_gaussmulti(n,mustars,taustars,sigmastars,seed=seed)
  
  # choose beta and lambda via the data
  beta = calc_beta(y)
  lambda = calc_lambda(y,nu,G)
  kappa_0 = calc_kappa_0(y)
  
  R = ncol(y)
  browser()
  
  # initialize via Mclust
  Cs_init = Mclust(y, G=G)$classification
  # res=Mclust(y, G=G)
  # res$parameters
  
  (log_marglik = log_true_marglik_gauss_multi(y, nu, lambda, alpha_0, kappa_0, beta, mustars))
  tic()
  posterior_sample = gibbs_sampling_gaussmulti(iters+burn_in, Cs_init=Cs_init, y, 
                                               alpha_0, kappa_0, nu, lambda, beta ,seed=seed)
  toc()
  posterior_sample$theta = posterior_sample$theta[-(1:burn_in),,]
  posterior_sample$C_mat = posterior_sample$C_mat[-(1:burn_in),]

  return(list(results=posterior_sample$theta,
              alloc_vec=posterior_sample$C_mat,
              name="GIBBS", y=y, truth=log_marglik))
}

# # combines the sampling of y with the gibbs sampling of mu
# y_theta_sampler_gaussmulti = function(n, nu, alpha_0,
#                                       mustars,sigmastars,taustars,
#                                       init,iters,burn_in = 2000,seed=2024){
#   #browser()
#   G = length(alpha_0)
#   y = rdata_gaussmulti(n,mustars,taustars,sigmastars,seed=seed)
#   
#   # choose beta and lambda via the data
#   beta = calc_beta(y)
#   lambda = calc_lambda(y,nu)
#   kappa_0 = calc_kappa_0(y)
#   
#   R = ncol(y)
#   #browser()
#   
#   # initialize via Mclust
#   Cs_init = Mclust(y, G=G)$classification
#   # res=Mclust(y, G=G)
#   # res$parameters
#   
#   (log_marglik = log_true_marglik_gauss_multi(y, nu, lambda, alpha_0, kappa_0, beta, mustars))
#   tic()
#   posterior_sample = gibbs_sampling_gaussmulti(iters+burn_in, Cs_init=Cs_init, y, 
#                                                alpha_0, kappa_0, nu, lambda, beta ,seed=seed)
#   toc()
#   posterior_sample$theta = posterior_sample$theta[-(1:burn_in),,]
#   posterior_sample$C_mat = posterior_sample$C_mat[-(1:burn_in),]
#   
#   return(list(results=posterior_sample$theta,
#               alloc_vec=posterior_sample$C_mat,
#               name="GIBBS", y=y, truth=log_marglik))
# }

# # combines the sampling of y with the gibbs sampling of mu
# y_theta_sampler_gaussmulti = function(n, nu, alpha_0,
#                                       mustars,sigmastars,taustars,
#                                       init,iters,burn_in = 2000,seed=2024){
#   #browser()
#   G = length(alpha_0)
#   y = rdata_gaussmulti(n,mustars,taustars,sigmastars,seed=seed)
#   
#   # choose beta and lambda via the data
#   beta = calc_beta(y)
#   lambda = calc_lambda(y,nu)
#   kappa_0 = calc_kappa_0(y)
#   
#   R = ncol(y)
#   #browser()
#   
#   # initialize via Mclust
#   Cs_init = Mclust(y, G=G)$classification
#   # res=Mclust(y, G=G)
#   # res$parameters
#   
#   (log_marglik = log_true_marglik_gauss_multi(y, nu, lambda, alpha_0, kappa_0, beta, mustars))
#   tic()
#   posterior_sample = gibbs_sampling_gaussmulti(iters+burn_in, Cs_init=Cs_init, y, 
#                                                alpha_0, kappa_0, nu, lambda, beta ,seed=seed)
#   toc()
#   posterior_sample$theta = posterior_sample$theta[-(1:burn_in),,]
#   posterior_sample$C_mat = posterior_sample$C_mat[-(1:burn_in),]
#   
#   return(list(results=posterior_sample$theta,
#               alloc_vec=posterior_sample$C_mat,
#               name="GIBBS", y=y, truth=log_marglik))
# }

# rgibbs_C_vec_gaussmulti = function(y,theta_i,sigma_tilde,taus_tilde){
#   
#   G=length(theta_i)
#   R=ncol(y)
#   n=nrow(y)
#   
#   num_cores <- detectCores()
#   num_use_cores = min(c(num_cores-2,9))
#   postprob_mat = simplify2array(mclapply(1:G,function(g) log(theta_i[[g]]$pi_g) + 
#                                            mvtnorm::dmvnorm(y, theta_i[[g]]$mu_g, 
#                                                             theta_i[[g]]$Sigma_g,log=TRUE),mc.cores = num_use_cores))
#   
#   
#   # tic()
#   # for(i in 1:20000){
#   #   maxlogrows1 = apply(postprob_mat,1,max)
#   # }
#   # toc()
#   # 
#   # 
#   # tic()
#   # for(i in 1:20000){
#   #  maxlogrows2 = do.call(pmax, c(as.data.frame(postprob_mat)))
#   # }
#   # toc()
#   # this one is faster
#   
#   # normalize to deal with potential numeric issues
#   maxlogrows = do.call(pmax, c(as.data.frame(postprob_mat)))
#   postprob_mat_normalized = exp(postprob_mat - maxlogrows) / rowSums(exp(postprob_mat - maxlogrows))
#   
#   # postprobmat_test = t(matrix(c(1/4,2/4,1/4,0/4,3/4,1/4,0/4,0/4),nrow=G))
#   # linearized sampler from the posterior (works if there are no ties)
#   
#   # matrix of cumulative probabilities
#   postprob_mat_cumul = postprob_mat_normalized + 
#     matrix(rowSums(simplify2array(
#       mclapply(1:G,
#                function(g) c(rep(0,n*g),
#                              rep(c(postprob_mat_normalized)[(n*(g-1)+1):(n*g)],G-g)), 
#                mc.cores = num_use_cores))),nrow=n)
#   
#   # random draw from the cumulative probabilities
#   postrob_mat_randomized = matrix(c(runif(n) <= postprob_mat_cumul + 0)*rep(G:1,each=n),nrow=n)
#   
#   # taking the argmax (linearized)
#   C_vec_i = G-do.call(pmax, c(as.data.frame(postrob_mat_randomized)))+1
#   browser()
#   #postprobvec = c(2/8,1/8,5/8)
#   #table(sapply(1:10000000,function(s) which.max(runif(1) <= cumsum(postprobvec))))/10000000
#   return(C_vec_i)
# }

# used to determine the hyperparameter beta of the prior
calc_beta = function(y){
  return(colMeans(y))
}

# the gibbs sampler
gibbs_sampling_gaussmulti = function(iters, Cs_init=Cs_init, y,
                                     alpha_0, kappa_0, nu, lambda, beta ,seed=seed){
  set.seed(seed)
  #browser()
  G = length(alpha_0)
  R = ncol(y)
  n = nrow(y)
  #browser()
  C_mat = matrix(nrow=iters,ncol=n)
  theta = array(dim=c(iters,G,1+R+(R*(R-1)/2+R))) 
  
  C_mat[1,] = Cs_init
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  
  tic()
  for(i in 2:iters){
    #tic()
    # if(i==21){
    #browser()
    # }
    #print(i)
    if(round(i/100)==i/100){
      print(i)
    }
    C_vec = C_mat[i-1,]
    theta_i = rgibbs_theta_gaussmulti(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=((seed*iters):((seed+1)*(iters)))[i-1])
    #theta_i = rgibbs_theta_gaussmulti(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=((seed*1000):((seed+1)*(1000)))[i-1])
    
    theta[i-1,,] = t(simplify2array(mclapply(1:G, function(g) c(theta_i[[g]]$mu_g,theta_i[[g]]$Sigma_g[upper.tri(theta_i[[g]]$Sigma_g,diag=TRUE)],theta_i[[g]]$pi_g),mc.cores = num_use_cores)))

    C_vec_i = rgibbs_C_vec_gaussmulti(y,theta_i,sigma_tilde,taus_tilde,seed=i)
    C_mat[i,] = C_vec_i
    #toc()
  }
  #browser()
  C_vec = C_mat[iters,]
  theta_i = rgibbs_theta_gaussmulti(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=iters)
  
  theta[iters,,] = t(simplify2array(mclapply(1:G, function(g) c(theta_i[[g]]$mu_g,theta_i[[g]]$Sigma_g[upper.tri(theta_i[[g]]$Sigma_g,diag=TRUE)],theta_i[[g]]$pi_g),mc.cores = num_use_cores)))
  
  
  #0.21*iters/60/60
  toc()
  #browser()
  
  posterior_sample = list(theta=theta,C_mat = C_mat)
  
  return(posterior_sample)
}

p0hat_gaussmulti <- function(y, theta, G) {
  
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  p0hats = simplify2array(mclapply(1:iters,
                       function(i) p0hat_i(y, theta[i,], G),
                       mc.cores = num_use_cores))
  #browser()
  return(mean(p0hats))
  
}

p0hat_i  <- function(y, theta, G) {
  browser()
  R = dim(y)[2]
  param_list = transform_to_params(theta,G,R)
  logcholsigmas = param_list$sigmas
  
  cholsigmas = lapply(1:G, function(g) logcholsigmas[[g]] - diag(diag(logcholsigmas[[g]])) + diag(diag(exp(logcholsigmas[[g]]))))
  for(g in 1:G){
    cholsigmas[[g]][upper.tri(cholsigmas[[g]])] = 0
  }
  
  sigmas = lapply(1:G, function(g) cholsigmas[[g]] %*% t(cholsigmas[[g]]))
  browser()
  postprob_mat = simplify2array(lapply(1:G,function(g) log(param_list$pis[g]) + 
                                           mvtnorm::dmvnorm(y, param_list$mus[g,], 
                                                            sigmas[[g]],log=TRUE)))
  maxlogrows = do.call(pmax, c(as.data.frame(postprob_mat)))
  postprob_mat_normalized = exp(postprob_mat - maxlogrows) / rowSums(exp(postprob_mat - maxlogrows))
  p0estim = mean(sapply(1:G, function(g) prod(1-postprob_mat_normalized[,g])))
  return(p0estim)
}

# sampling conditioned on C
rgibbs_C_vec_gaussmulti = function(y,theta_i,sigma_tilde,taus_tilde,seed){
  set.seed(seed)
  G=length(theta_i)
  R=ncol(y)
  n=nrow(y)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  #browser()
  postprob_mat = simplify2array(mclapply(1:G,function(g) log(theta_i[[g]]$pi_g) + 
                                           mvtnorm::dmvnorm(y, theta_i[[g]]$mu_g, 
                                                            theta_i[[g]]$Sigma_g,log=TRUE),mc.cores = num_use_cores))
  
  
  # tic()
  # for(i in 1:20000){
  #   maxlogrows1 = apply(postprob_mat,1,max)
  # }
  # toc()
  # 
  # 
  # tic()
  # for(i in 1:20000){
  #  maxlogrows2 = do.call(pmax, c(as.data.frame(postprob_mat)))
  # }
  # toc()
  # this one is faster
  
  # normalize to deal with potential numeric issues
  maxlogrows = do.call(pmax, c(as.data.frame(postprob_mat)))
  postprob_mat_normalized = exp(postprob_mat - maxlogrows) / rowSums(exp(postprob_mat - maxlogrows))
  p0estim = mean(sapply(1:G, function(g) prod(1-postprob_mat_normalized[,g])))
  # postprobmat_test = t(matrix(c(1/4,2/4,1/4,0/4,3/4,1/4,0/4,0/4),nrow=G))
  # linearized sampler from the posterior (works if there are no ties)
  
  # matrix of cumulative probabilities
  postprob_mat_cumul = postprob_mat_normalized + 
    matrix(rowSums(simplify2array(
      mclapply(1:G,
               function(g) c(rep(0,n*g),
                             rep(c(postprob_mat_normalized)[(n*(g-1)+1):(n*g)],G-g)), 
               mc.cores = num_use_cores))),nrow=n)
  
  # random draw from the cumulative probabilities
  postrob_mat_randomized = matrix(c(runif(n) <= postprob_mat_cumul + 0)*rep(G:1,each=n),nrow=n)
  
  # taking the argmax (linearized)
  C_vec_i = G-do.call(pmax, c(as.data.frame(postrob_mat_randomized)))+1
  #browser()
  #postprobvec = c(2/8,1/8,5/8)
  #table(sapply(1:10000000,function(s) which.max(runif(1) <= cumsum(postprobvec))))/10000000
  return(C_vec_i)
}

# sampling conditioned on theta
rgibbs_theta_gaussmulti = function(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed){
  
  G=length(alpha_0)
  alpha = alpha_0 + sapply(1:G, function(g) sum(C_vec==g))
  pis = rdirichlet(1, alpha) 
  #browser()
  rgibbs_theta_g = function(g){
    set.seed(seed+g)
    # define setting for this particular class
    n_g = sum(C_vec==g)
    y_g = y[C_vec==g,]
    beta_g = beta#s[,g]
    lambda_g = lambda
    nu_g = nu
    alpha_g = alpha_0
    
    # updated parameters on the conjugate prior (see Gelman et al)
    
    if(n_g==0){
      # set to prior if there is no information
      nu_ng = nu_g
      lambda_ng = lambda_g
      beta_ng = beta_g
      kappa_ng = kappa_0
    } else if (n_g==1){
      ybar_g = y_g
      S_ng = t(y_g)%*% t(t(y_g))
      
      kappa_ng = kappa_0 + n_g
      nu_ng = nu_g + n_g
      beta_ng = (1/(1+n_g/kappa_0))*beta_g+ybar_g*n_g/(kappa_0+n_g)
      lambda_ng = lambda_g + S_ng + t(t(t(ybar_g))-t(beta_g))%*%(t(t(ybar_g))-t(beta_g)) *(kappa_0*n_g)/(kappa_0+n_g)
      beta_ng = t(beta_ng)
      
    } else{
      ybar_g = colMeans(y_g)
      S_ng = var(y_g)*(n_g-1)
      beta_ng = (1/(1+n_g/kappa_0))*beta_g+ybar_g*n_g/(kappa_0+n_g)
      kappa_ng = kappa_0 + n_g
      nu_ng = nu_g + n_g
      lambda_ng = lambda_g + S_ng + (t(t(ybar_g))-t(t(beta_g)))%*%t(t(t(ybar_g))-t(t(beta_g))) *(kappa_0*n_g)/(kappa_0+n_g)
    }
    
    # draw from the posterior given Z
    Sigma_g = CholWishart::rInvWishart(1, nu_ng, lambda_ng)[,,1]
    mu_g = mvtnorm::rmvnorm(1,mean=beta_ng,
                            sigma=Sigma_g/kappa_ng)
    return(list(pi_g=pis[g],mu_g=mu_g,Sigma_g=Sigma_g))
  }
  #browser()
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  theta_i = mclapply(1:G,function(g) rgibbs_theta_g(g),mc.cores = num_use_cores)
  return(theta_i)
  
}


# samples from the prior and returns the loglikelihood
# the prior mean beta needs to be set within the function, since it depends on the data

prior_sampler_marglik_gaussmulti = function(y, nu, iters, alpha_0, seed=2024){
  #browser()
  
  # specification of the hyperparameters that determine the prior
  R = ncol(y)
  G = length(alpha_0)
  
  beta = calc_beta(y)
  lambda = calc_lambda(y, nu,G)
  kappa_0 = calc_kappa_0(y)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  
  # sample from the prior
  prior_sample_taus = t(simplify2array(mclapply(1:iters,function(i) rdirichlet(1, alpha_0)[1,],mc.cores = num_use_cores)))
  uppertri_length = R*(R-1)/2+R
  prior_sample_sigmas = t(simplify2array(mclapply(1:iters, function(i) c(sapply(1:G, function(g) CholWishart::rInvWishart(1,nu,lambda)[,,1][upper.tri(matrix(ncol=R,nrow=R), diag=TRUE)])))))
  prior_sample_means = t(simplify2array(mclapply(1:iters, 
                                                 function(i) c(sapply(1:G, function(g) mvtnorm::rmvnorm(1,mean=beta,
                                                                                                        sigma=transform_to_matrix(prior_sample_sigmas[i,((g-1)*uppertri_length+1):(g*uppertri_length)],R)/kappa_0))), 
                                                 mc.cores = num_use_cores)))
  theta = cbind(prior_sample_means,prior_sample_sigmas, prior_sample_taus)
  # tic()
  # logliks = loglik_gmm_i(y,theta[1,],G,1)
  # toc()
  mc_estim_partial = loglik_gmm_gaussmulti(y,theta,G)
  
  #par(mfrow=c(1,1))
  #plot(100000:iters,(log(cumsum(exp(logliks-min(logliks)))/(1:iters)))[100000:iters]+min(logliks))
  
  
  
  
  #abline(h=log_marglik,col="red")
  
  # TODO define mc_estim_partial ?
  return(mc_estim_partial)
}

# calculates the true marginal likelihood given the data
log_true_marglik_gauss_multi = function(y, nu, lambda, alpha_0, kappa_0, beta, mustars){
  
  R = ncol(lambda)
  G = length(alpha_0)
  C_vec = G-rowSums(sapply(1:G, function(g) rowSums(y) <= (R*mustars[1,]+100)[g]))+1
  #C_vec = G-(y[,1] <= beta) 
  # we can do this because the clusters are already extremely well seperated
  
  log_marglik = 0
  for(g in 1:G){
    #browser()
    # define setting for this particular class
    n_g = sum(C_vec==g)
    y_g = y[C_vec==g,]
    beta_g = beta#s[,g]
    lambda_g = lambda
    nu_g = nu
    alpha_g = alpha_0
    
    # updated parameters on the conjugate prior (see Gelman et al)
    ybar_g = colMeans(y_g)
    beta_ng = (1/(1+n_g/kappa_0))*beta_g+ybar_g*n_g/(kappa_0+n_g)
    kappa_ng = kappa_0 + n_g
    nu_ng = nu_g + n_g
    
    S_ng = var(y_g)*(n_g-1)
    lambda_ng = lambda_g + S_ng + (t(t(ybar_g))-t(t(beta_g)))%*%t(t(t(ybar_g))-t(t(beta_g))) *(kappa_0*n_g)/(kappa_0+n_g)
    
    # compute marginal likelihood via Bayes' rule
    mu_g = mustars[,g]
    sigma_g = diag(R) # result should not change with sigma_g,mu_g
    #log_prior_g = dinvwishart(sigma_g,nu_g,lambda_g,log=TRUE) + mvtnorm::dmvnorm(mu_g,beta_g,sigma_g/kappa_0,log=TRUE)
    #loglik_g = sum(mvtnorm::dmvnorm(y_g,mu_g,sigma_g,log=TRUE))  
    #logpost_g = dinvwishart(sigma_g,nu_ng,lambda_ng,log=TRUE) + mvtnorm::dmvnorm(mu_g,beta_ng,sigma_g/kappa_ng,log=TRUE)
    #(log_marglik_g = log_prior_g + loglik_g - logpost_g)
    #browser()
    (log_marglik_g = (nu_g/2)*sum(log(eigen(lambda_g)$values))-(nu_g*R/2)*log(2)-lmvgamma(nu_g/2,R)+
        (nu_ng*R/2)*log(2)+lmvgamma(nu_ng/2,R)+(-nu_ng/2)*sum(log(eigen(lambda_ng)$values))-
        (n_g*R)*log(sqrt(2*pi))-log(kappa_ng)*R/2+log(kappa_0)*R/2)
    
    log_marglik = log_marglik + log_marglik_g
  }
  
  # log_prior_given_c = sum(sapply((1:G)[taus_g!=0], function(g) log(taus_g[g])*sum(C_vec==g)))
  ns = sapply(1:G,function(g) sum(C_vec==g))
  lprior_z = lgamma(sum(alpha_0)) - lgamma(sum(alpha_0)+n) + sum(lgamma(ns+alpha_0)-lgamma(alpha_0))
  
  log_marglik = log_marglik +lprior_z + lfactorial(G) # only if there is a symmetric prior
  return(log_marglik)
}

# simulates the data given the parameters
rdata_gaussmulti = function(n,mustars,taus_stars,sigmas_stars,seed=2024){
  set.seed(seed)
  #browser()
  G = ncol(mustars)
  R = nrow(mustars)
  
  C = apply(rmultinom(n=n,size=1,prob=taus_stars),2,function(x) which.max(x))
  y = matrix(nrow=n, ncol=R)
  for(g in 1:G){
    y[C==g,] = mvtnorm::rmvnorm(sum(C==g),mean=mustars[,g],sigma=sigmas_stars[,,g])
  }
  #browser()
  return(y)
}


logprior_gmm_gaussmulti <- function(theta, G, nu, lambda, beta, alpha_0, kappa_0) {
  
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  logpriors = mclapply(1:iters,
                       function(i) logprior_gmm_i_gaussmulti(theta[i,], G, nu, lambda, beta, alpha_0, kappa_0),
                       mc.cores = num_use_cores)
  return(simplify2array(logpriors))
  
}

logprior_gmm_i_gaussmulti <- function(theta, G, nu, lambda, beta, alpha_0, kappa_0) {
  browser()
  R = length(beta)
  
  param_list = transform_to_params(theta,G,R)
  
  pis = param_list$pis
  mus = param_list$mus
  sigmas = param_list$sigmas
  
  log_pis = log((sum(pis[-1])<=1)&(mean(pis>0)==1))
  if(!is.infinite(log_pis)){
    log_pis = log((sum(pis[-1])<=1)&mean(pis>0))+ddirichlet(x=pis,alpha=alpha_0,log=TRUE)
  }
  
  log_sigmas = try(sum(dInvWishart(simplify2array(sigmas),df=nu,Sigma=lambda,log=TRUE)))
  if(is.character(log_sigmas)){
    log_sigmas = -Inf
  }
  
  log_mus = try(sum(sapply(1:G, function(g) mvtnorm::dmvnorm(mus[g,],mean=beta,sigmas[[g]]/kappa_0,log=TRUE))))
  if(is.character(log_mus)){
    log_mus = -Inf
  }
  
  return(log_pis+log_sigmas+log_mus)
  
}

loglik_gmm_gaussmulti <- function(y, theta, G){
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  logliks = mclapply(1:iters,
                     function(i) loglik_gmm_i_gaussmulti(y, theta[i,], G),
                     mc.cores = num_use_cores)
  return(simplify2array(logliks ))
}

transform_to_matrix = function(uppertri,R){
  A = matrix(0, R, R)
  A[upper.tri(A, diag = TRUE)] = uppertri
  A = A + t(A)
  diag(A) = diag(A) / 2
  return(A)
}

transform_to_params = function(theta,G,R){
  
  index = 0
  mus = t(matrix(theta[(index+1):(index+R*G)], R, G))
  index = index+R*G
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  
  sigmas = mclapply(1:G, 
                    function(g) transform_to_matrix(
                      theta[(index+1+(g-1)*(R*(R-1)/2+R)):(index+g*(R*(R-1)/2+R))],R),
                    mc.cores = num_use_cores )
  
  index = index + (R*(R-1)+R)*G
  
  pis = theta[(length(theta)-G+1):length(theta)]
  
  return(list(mus=mus,sigmas=sigmas,pis=pis))
}

# ASSUMING SIGMA IS THE ACTUAL COVARIANCE MATRIX (NOT ITS INVERSE)
loglik_gmm_i_gaussmulti <- function(y, theta, G) {
  
  #print("likelihood function: start")
  browser()
  n = dim(y)[1]
  R = dim(y)[2]
  
  param_list = transform_to_params(theta,G,R)
  
  pis = param_list$pis
  mus = param_list$mus
  sigmas = param_list$sigmas
  
  isinfinite = !((sum(pis[-1])<=1)&(mean(pis>0)==1)&(!is.character(try(dInvWishart(simplify2array(sigmas),df=R,Sigma=diag(R),log=TRUE)))))
  if(isinfinite){
    log_sum_lik = -Inf # check for invalid values (can happen with MC)
  } else{
    log_sum_lik = mvnfast::dmixn(X=y,mus,sigmas,pis,log=TRUE)
  }
  return(sum(log_sum_lik))
}


logposty_gaussmulti = function(thetas, G, y, nu, alpha_0, kappa_0){
  browser()
  beta = calc_beta(y)
  lambda = calc_lambda(y,nu,G)
  kappa_0 = calc_kappa_0(y)
  
  iters = nrow(thetas)
  
  logliks  = loglik_gmm_gaussmulti(y,thetas,G)
  logpriors = logprior_gmm_gaussmulti(thetas, G, nu, lambda, beta, alpha_0, kappa_0) 
  return(logliks + logpriors)
}


### VVV model (no transform) ###

### VVV model (transform) ###

# combines the sampling of y with the gibbs sampling of mu
y_theta_sampler_gaussmulti_transformed = function(n, nu, alpha_0,
                                      mustars,sigmastars,taustars,
                                      init,iters,burn_in = 2000,seed=2024,y=NULL){
  #browser()
  G = length(alpha_0)
  if(is.null(y)){
    y = rdata_gaussmulti(n,mustars,taustars,sigmastars,seed=seed)    
  }
  
  # choose beta and lambda via the data
  beta = calc_beta(y)
  lambda = calc_lambda(y,nu,G)
  kappa_0 = calc_kappa_0(y)
  
  R = ncol(y)
  #browser()
  
  # initialize via Mclust
  Cs_init = Mclust(y, G=G)$classification
  # res=Mclust(y, G=G)
  # res$parameters
  if(is.null(y)){
    (log_marglik = log_true_marglik_gauss_multi(y, nu, lambda, alpha_0, kappa_0, beta, mustars))
  } else{
    log_marglik = 0 # presume truth is unknown if y is known
  }

  tic()
  posterior_sample = gibbs_sampling_gaussmulti_transformed(iters+burn_in, Cs_init=Cs_init, y, 
                                               alpha_0, kappa_0, nu, lambda, beta ,seed=seed)
  toc()
  posterior_sample$theta = posterior_sample$theta[-(1:burn_in),,]
  posterior_sample$C_mat = posterior_sample$C_mat[-(1:burn_in),]
  
  return(list(results=posterior_sample$theta,
              alloc_vec=posterior_sample$C_mat,
              name="GIBBS", y=y, truth=log_marglik))
}

# the gibbs sampler
gibbs_sampling_gaussmulti_transformed = function(iters, Cs_init=Cs_init, y,
                                     alpha_0, kappa_0, nu, lambda, beta ,seed=seed){
  set.seed(seed)
  #browser()
  G = length(alpha_0)
  R = ncol(y)
  n = nrow(y)
  #browser()
  C_mat = matrix(nrow=iters,ncol=n)
  theta = array(dim=c(iters,G,1+R+(R*(R-1)/2+R))) 
  
  #C_mat[1,] = Cs_init
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  
  try(load("data/C_matG5R6.RData"))
  # if(iters==11000){ # TODO REMOVE
  #   load("data/theta.RData")
  #   C_mat = C_mat[1:iters,]
  #   theta = theta[1:iters,,]
  # } else{
  # reload if possible (to save time)
  #browser()
  if((max(C_mat[!is.na(C_mat)])!=G)|(ncol(C_mat)!=n)|(nrow(C_mat)!=iters)|
     (which.max(rowSums(is.na(C_mat)))==2)|(mean(C_mat[1,] == Cs_init)<1)){
    C_mat = matrix(nrow=iters,ncol=n)
    C_mat[1,] = Cs_init
  } else{
    load("data/thetaG5R6.RData")
  }
  #}
  

  tic()
  if(max(rowSums(is.na(C_mat)))>0){
    for(i in which.max(rowSums(is.na(C_mat))):iters){
      #tic()
      # if(i==21){
      #browser()
      # }
      #print(i)

      if(round(i/100)==i/100){
        save(C_mat,file="data/C_matG5R6.RData")
        save(theta,file="data/thetaG5R6.RData")
        print(i)
      }
      C_vec = C_mat[i-1,]
      theta_i = rgibbs_theta_gaussmulti(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=((seed*iters):((seed+1)*(iters)))[i-1])
      #theta_i = rgibbs_theta_gaussmulti(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=((seed*1000):((seed+1)*(1000)))[i-1])
      #browser()
      theta[i-1,,] = t(simplify2array(mclapply(1:G, function(g) c(theta_i[[g]]$mu_g,map_to_R(theta_i[[g]]$Sigma_g)[upper.tri(theta_i[[g]]$Sigma_g,diag=TRUE)],theta_i[[g]]$pi_g),mc.cores = num_use_cores)))
      #browser()
      C_vec_i = rgibbs_C_vec_gaussmulti(y,theta_i,sigma_tilde,taus_tilde,seed=i)
      C_mat[i,] = C_vec_i
      #toc()
    }
    #browser()
    C_vec = C_mat[iters,]
    
    theta_i = rgibbs_theta_gaussmulti(y, C_vec, alpha_0, kappa_0, nu, lambda, beta, seed=iters)
    
    theta[iters,,] = t(simplify2array(mclapply(1:G, function(g) c(theta_i[[g]]$mu_g,map_to_R(theta_i[[g]]$Sigma_g)[upper.tri(theta_i[[g]]$Sigma_g,diag=TRUE)],theta_i[[g]]$pi_g),mc.cores = num_use_cores)))
  
  save(C_mat,file="data/C_matG5R6.RData")
  save(theta,file="data/thetaG5R6.RData")
  }
  #0.21*iters/60/60
  toc()
  browser()
  
  posterior_sample = list(theta=theta,C_mat = C_mat)
  
  return(posterior_sample)
}

logprior_gmm_gaussmulti_transformed <- function(theta, G, nu, lambda, beta, alpha_0, kappa_0) {
  
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  #browser()
  logpriors = mclapply(1:iters,
                       function(i) logprior_gmm_i_gaussmulti_transformed(theta[i,], G, nu, lambda, beta, alpha_0, kappa_0),
                       mc.cores = num_use_cores)
  return(simplify2array(logpriors))
  
}

logprior_gmm_i_gaussmulti_transformed <- function(theta, G, nu, lambda, beta, alpha_0, kappa_0) {
  browser()
  R = length(beta)
  
  param_list = transform_to_params(theta,G,R)
  
  pis = param_list$pis
  mus = param_list$mus
  # sigmas = param_list$sigmas
  
  logcholsigmas = param_list$sigmas
  
  cholsigmas = lapply(1:G, function(g) logcholsigmas[[g]] - diag(diag(logcholsigmas[[g]])) + diag(diag(exp(logcholsigmas[[g]]))))
  for(g in 1:G){
    cholsigmas[[g]][upper.tri(cholsigmas[[g]])] = 0
  }
  
  sigmas = lapply(1:G, function(g) cholsigmas[[g]] %*% t(cholsigmas[[g]]))
  browser() 
  # cholsigmas = lapply(1:G,function (g) t(chol(sigmas[[g]])))
  # logcholsigmas = lapply(1:G,function (g) t(chol(sigmas[[g]]))-diag(diag(chol(sigmas[[g]])))+diag(log(diag(chol(sigmas[[g]])))))
  logdetcholsigmas = sum(sapply(1:G, function(g) R*log(2)+sum(((R+1):2)*log(diag(cholsigmas[[g]])))))
  
  log_pis = log((sum(pis[-1])<=1)&(mean(pis>0)==1))
  if(!is.infinite(log_pis)){
    log_pis = log((sum(pis[-1])<=1)&mean(pis>0))+ddirichlet(x=pis,alpha=alpha_0,log=TRUE)
  }
  
  log_sigmas = try(sum(dInvWishart(simplify2array(sigmas),df=nu,Sigma=lambda,log=TRUE)))
  if(is.character(log_sigmas)){
    log_sigmas = -Inf
  }
  
  log_mus = try(sum(sapply(1:G, function(g) mvtnorm::dmvnorm(mus[g,],mean=beta,sigmas[[g]]/kappa_0,log=TRUE))))
  if(is.character(log_mus)){
    log_mus = -Inf
  }
  
  return(log_pis+log_sigmas+log_mus+logdetcholsigmas)
  
}

loglik_gmm_gaussmulti_transformed <- function(y, theta, G){
  #browser()
  iters = nrow(theta)
  
  num_cores <- detectCores()
  num_use_cores = min(c(num_cores-2,9))
  logliks = mclapply(1:iters,
                     function(i) loglik_gmm_i_gaussmulti_transformed(y, theta[i,], G),
                     mc.cores = num_use_cores)
  return(simplify2array(logliks ))
}

loglik_gmm_i_gaussmulti_transformed <- function(y, theta, G) {
  
  #print("likelihood function: start")
  browser()
  n = dim(y)[1]
  R = dim(y)[2]
  
  param_list = transform_to_params(theta,G,R)
  
  pis = param_list$pis
  mus = param_list$mus
  logcholsigmas = param_list$sigmas
  
  cholsigmas = lapply(1:G, function(g) logcholsigmas[[g]] - diag(diag(logcholsigmas[[g]])) + diag(diag(exp(logcholsigmas[[g]]))))
  for(g in 1:G){
    cholsigmas[[g]][upper.tri(cholsigmas[[g]])] = 0
  }
  
  sigmas = lapply(1:G, function(g) cholsigmas[[g]] %*% t(cholsigmas[[g]]))
  #browser()
  isinfinite = !((sum(pis[-1])<=1)&(mean(pis>0)==1)&(!is.character(try(dInvWishart(simplify2array(sigmas),df=R,Sigma=diag(R),log=TRUE)))))
  if(isinfinite){
    log_sum_lik = -Inf # check for invalid values (can happen with MC)
  } else{
    log_sum_lik = mvnfast::dmixn(X=as.matrix(y),mus,sigmas,pis,log=TRUE)
  }
  return(sum(log_sum_lik))
}


logposty_gaussmulti_transformed = function(thetas, G, y, nu, alpha_0, kappa_0){
  #browser()
  beta = calc_beta(y)
  lambda = calc_lambda(y,nu,G)
  kappa_0 = calc_kappa_0(y)
  
  iters = nrow(thetas)
  
  logliks  = loglik_gmm_gaussmulti_transformed(y,thetas,G)
  logpriors = logprior_gmm_gaussmulti_transformed(thetas, G, nu, lambda, beta, alpha_0, kappa_0) 
  return(logliks + logpriors)
}


map_to_R = function(covmat){
  L = chol(covmat)
  diag(L) = log(diag(L))
  return(L)
}

### VVV model (transform) ###