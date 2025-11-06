#source("functions/thames_gmm.R")
#source("functions/topsort_funcs.R")
library(label.switching)
library(OptimalDesign)  # to compute the min-vol covering ellipse
library(bayesmix)
library(tidyverse)
library(reshape2)
library(combinat)
library(cvCovEst)
library(glasso)
library(sparsediscrim)
library(igraph)
library(thamesmix)

#needed for relabel
complete.normal.loglikelihood<-function(x,z,pars){
  #x: denotes the n data points
  #z: denotes an allocation vector (size=n)
  #pars: K\times 3 vector of means,variance, weights
  # pars[k,1]: corresponds to the mean of component k
  # pars[k,2]: corresponds to the variance of component k
  # pars[k,3]: corresponds to the weight of component k
  g <- dim(pars)[1]
  n <- length(x)
  logl<- rep(0, n)
  logpi <- log(pars[,3])
  mean <- pars[,1]
  sigma <- sqrt(pars[,2])
  logl<-logpi[z] + dnorm(x,mean = mean[z],sd = sigma[z],log = TRUE)
  return(sum(logl))
}

# relabel = function(algs,sims,n,y,R,m,G,iters,
#                    z_dummy = matrix(rep(1,n*iters),ncol=n)){
relabel = function(lps, loglik_partial, algs, sims, G, iters, z_dummy, y, seed=2024){
  # computes a list of labels for different relabeling algorithms
  # z_dummy:  Posterior allocation vectors;
  #           Only used for ECR-type algorithms, not necessary for Steven's alg
  # algs:     One of "STEPHENS","ECR","PRA", "AIC","STEPHENS","SJW" and others
  
  #browser()
  set.seed(seed)
  prapivot = array(dim=dim(sims)[-1])
  prapivot[,] = sims[which.max(lps),,]
  # needed for stevens (only if loglik_partial is defined)
  if(is.null(dim(y))){
    n=length(y)
  }else{
    n=nrow(y)
  }
  p=array(1/G,dim=c(iters,n,G))
  if(!is.null(loglik_partial)){
    lls_partial = loglik_partial(sims,G)
    p = exp(lls_partial)
    sum = p
    for(g in 1:G){
      for(k in (1:G)[-g]){
        sum[,,g] = sum[,,g] + p[,,k]
      }
    }
    p = p / sum
  } 
  #browser()
  ls<-label.switching(method=algs,
                      zpivot=z_dummy[which.max(lps),],z = z_dummy,K = G, data = y,
                      prapivot = prapivot,mcmc = sims,constraint=1,p=p,
                      complete = complete.normal.loglikelihood)
  #browser()
  return(ls)
}

# melts the sims into a matrix
melt_sims = function(sims, G, num_R){
  #browser()
  params = matrix(nrow=dim(sims)[1],ncol=1)
  for(index in 1:G){
    params = cbind(params,sims[,index,])
  }
  #browser()
  params = params[,-1]
  if(num_R>1){
    index=0
    # the parameters themselves are multidimensional in this case
    #browser()
    if(dim(params)[2]==(G+2*R*G)){
      #browser()
      mus = params[,rep(seq(1,G*dim(sims)[3]-(dim(sims)[3]-1),dim(sims)[3]),each=num_R)+rep(0:(num_R-1),G)]
      index=num_R
      sigmas = params[,rep(seq(1+index,G*dim(sims)[3]-(dim(sims)[3]-1)+index,dim(sims)[3]),each=num_R)+rep(0:(num_R-1),G)]
      index=dim(sims)[3] -1
      pis = params[,rep(seq(1+index,G*dim(sims)[3]-(dim(sims)[3]-1)+index,dim(sims)[3]),each=1)+rep(0:(0),G)]
      params = cbind(mus,sigmas,pis)
    } else{
      mus = params[,rep(seq(1,G*dim(sims)[3]-(dim(sims)[3]-1),dim(sims)[3]),each=num_R)+rep(0:(num_R-1),G)]
      index=num_R
      sigmas = params[,rep(seq(1+index,G*dim(sims)[3]-(dim(sims)[3]-1)+index,dim(sims)[3]),each=num_R*(num_R-1)/2+num_R)+rep(0:(num_R*(num_R-1)/2+num_R-1),G)]
      index = dim(sims)[3] -1
      pis = params[,rep(seq(1+index,G*dim(sims)[3]-(dim(sims)[3]-1)+index,dim(sims)[3]),each=1)+rep(0:(0),G)]
      params = cbind(mus,sigmas,pis)
    }
  } else{
    params = params[,rep(seq(1,G*dim(sims)[3]-(dim(sims)[3]-1),dim(sims)[3]),dim(sims)[3])+
                      rep(0:(dim(sims)[3]-1),each=G)]
  }
  
  #browser()
  return(params)
}

#relabels the parameters for one relabeling algorithm (type)
relabel_params = function(sims,new_labels,G,type,num_R){
  #browser()
  orders = new_labels$permutations[[type, exact = FALSE]]
  orders = t(matrix(c(t(orders)) + rep(0:(dim(sims)[1]-1),each=G)*G,nrow=G))
  for(d in 1:dim(sims)[3]){
    parm = sims[,,d]
    sims[,,d] = matrix(c(t(parm))[orders],ncol=G)
  }
  # params = melt_sims(sims,G,num_R)
  return(list(sims=sims))
}

compute_thames = function(ellipse,params,lps,G,iters,type,logpost,num_R,seed=2024,prior_sampler=NULL,y=NULL,limit=NULL,num_var_g=3,sims=NULL){
  # Computes the THAMES
  # type: "simple" is the efficient version of the mixture THAMES
  #       "standard" is the vanilla THAMES
  #       "permutations" is the mixtures of G! THAMES, each concentrated on 
  #       a different posterior mode
  #browser()
  
  if(type=="mc"){
    browser()
    naive_mc_likelihood = prior_sampler(y=y, G=G, iters=2*iters)
    naive_mc_point_est = log(mean(exp(naive_mc_likelihood - max(naive_mc_likelihood)))) + max(naive_mc_likelihood)
    naive_mc_log_se = log(sd(exp(naive_mc_likelihood - max(naive_mc_likelihood)))/sqrt(2*iters)) + max(naive_mc_likelihood)
    naive_mc_cv = exp(naive_mc_log_se - naive_mc_point_est)
    
    trunc_quantile = function(p,ratio){
      alpha = - 1/ratio
      (qnorm(p=pnorm(alpha)+p*(1-pnorm(alpha))) * (ratio))+1
    } # Calculates quantiles of the truncated normal
    
    naive_mc_U = naive_mc_point_est + log(trunc_quantile(0.025,naive_mc_cv))
    naive_mc_mid = naive_mc_point_est
    naive_mc_L = naive_mc_point_est + log(trunc_quantile(0.975,naive_mc_cv))
    return(list(log_zhat_inv_L=-naive_mc_L,
                log_zhat_inv=-naive_mc_mid,
                log_zhat_inv_U=-naive_mc_U,
                log_cor=NULL,
                len_perms=NULL,
                alpha=NULL,
                c_opt=NULL,
                etas=params))
  } else if(type=="bridge"){
    #browser()
    if(dim(sims)[3]==1){
      logpost_bridge = function(param.row,data){
        tinyparam = rbind(param.row,param.row)
        tinysims = array(c(tinyparam),dim=c(2,dim(sims)[2:3]))
        return(logpost(tinysims)[1])
      }
    } else{
      logpost_bridge = function(param.row,data) logpost(rbind(param.row,param.row),G)[1]
    }
    library(bridgesampling)
    ub = rep(Inf,ncol(params))
    lb = rep(-Inf,ncol(params))
    param = params
    names(lb) <- names (ub) <- colnames(param) <- names(as.data.frame(param))
    bridge=bridge_sampler(samples=param, 
                          log_posterior=logpost_bridge,
                          lb=lb,ub=ub,data=NULL,silent=TRUE)
    
    trunc_quantile = function(p,ratio){
      alpha = - 1/ratio
      (qnorm(p=pnorm(alpha)+p*(1-pnorm(alpha))) * (ratio))+1
    } #Calculates quantiles of the truncated normal
    
    bridge_logml = bridge$logml
    bridge_cv = error_measures(bridge)$cv # bridge$cv
    
    bridge_L = bridge_logml + log(trunc_quantile(0.025,bridge_cv))
    bridge_mid = bridge_logml
    bridge_U = bridge_logml + log(trunc_quantile(0.975,bridge_cv))
    #browser()
    return(list(log_zhat_inv_L=-bridge_L,
                  log_zhat_inv=-bridge$logml,
                  log_zhat_inv_U=-bridge_U,
                  log_cor=NULL,
                  len_perms=NULL,
                  alpha=NULL,
                  c_opt=NULL,
                  etas=params))
  } else{
    
    #browser()
    
    lml_thames_marginal <- thamesmix::thames_mixtures(logpost=logpost,
      sims=sims,
      type=type,
      seed=seed
    )
    
    if((dim(sims)[3]>1)&(num_R==1)&(type=="simple")){
      sims_copy = sims
      dummy_sims = array(dim=c(floor(sqrt(dim(sims)[1]))^2,dim(sims)[2],dim(sims)[3]))
      for(g in 1:G){
        # grid = expand.grid(mu_g = seq(min(sims[,g,1])-0.1,max(sims[,g,1])+0.1,length.out=floor(sqrt(dim(sims)[1]))),
        #                    sigmasqu_g = seq(min(exp(sims[,g,2])),max(exp(sims[,g,2])),length.out=floor(sqrt(dim(sims)[1]))))
        grid = expand.grid(mu_g = seq(min(sims[,,1])-0.1,max(sims[,,1])+0.1,length.out=floor(sqrt(dim(sims)[1]))),
                           sigmasqu_g = seq(min(exp(sims[,,2])),max(exp(sims[,,2])),length.out=floor(sqrt(dim(sims)[1]))))
        dummy_sims[,g,1] = grid$mu_g
        dummy_sims[,g,2] = log(grid$sigmasqu_g)
      }
      # Convert the matrix into a data frame for ggplot
      dummy_params_f_transform = matrix(thamesmix:::reorder_by_lda(lml_thames_marginal$scaling,G,dummy_sims)$W,ncol=G)
      dummy_mat = cbind(dummy_sims[,1:G,1],exp(dummy_sims[,1:G,2]),dummy_params_f_transform)
      #browser()
      params_f_transform = cbind(params[,1:G],exp(params[,(G+1):(2*G)]),params[,(2*G+1):(3*G-1)], 
                                 matrix(thamesmix:::reorder_by_lda(lml_thames_marginal$scaling,G,sims)$W,ncol=G))
      params_f_transform = cbind(params_f_transform[1:nrow(dummy_mat),],dummy_mat)
    } else{
      params_f_transform = NULL
    }
    
    return(c(lml_thames_marginal,
             list(params_f_transform =  params_f_transform)))
  }
}
