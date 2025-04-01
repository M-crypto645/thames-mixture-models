source("functions/thames_gmm.R")
source("functions/topsort_funcs.R")
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
  #p_list_unnorm = apply(sims,1,function(theta) loglik_gmm_partial(y,c(theta[,3],theta[,1],sqrt(theta[,2])),G),simplify = 'FALSE')
  # p_list <- sapply(p_list_unnorm,function(x){exp(x)/rowSums(exp(x))},simplify=FALSE)
  # p <- aperm(sapply(p_list,function(x){x},simplify='array'),c(3,1,2))
  
  ### Plotting stuff ###
  # round_y = sort(unique(round(y)))
  # break_vec=c(sort(c(sort(unique(round(y))),
  #                      sort(unique(round(y)))-5,
  #                      sort(unique(round(y)))+5)))
  # hist(y,freq=TRUE,axes = FALSE,cex.main=2,cex.lab=1.5,cex.sub=1.5,xlab=NULL,cex.axis=1.5,breaks=break_vec)
  # zpivot=z_dummy[which.max(lps_marginal),]
  # axis(2,cex.axis=1.5)
  # axis(1,pos=-6.3,cex.axis=1.5,at=c(100,break_vec[seq(5,length(break_vec)-3,3)],2500))
  # rug(y[y<150],col="green",ticksize = .1,pos=-6,lwd=1.5)
  # colors=viridis::inferno(25)
  # library(RColorBrewer)
  # n <- 60
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # colors=col_vector[c(1:5,19:25,44:60)]
  # #colors=Polychrome::createPalette(25,seedcolors = c("#ff0000","#0000ff","#00ff00"),M=100000)
  # library(Polychrome)
  # test=function(){
  #   for(i in (2:24)){
  #     rug(y[(y>(round_y[i-1]+50))&(y<(round_y[i]+50))],col=colors[i],ticksize = .1,pos=-6,lwd=1.5)
  #   }
  # }
  # test()
  
  #hist(y,axes = FALSE,cex.main=2,cex.lab=1.5,cex.sub=1.5,xlab=NULL,cex.axis=1.5)
  #zpivot=z_dummy[which.max(lps_marginal),]
  # axis(1,pos=-5.125,cex.axis=1.5)
  # axis(2,cex.axis=1.5)
  # rug(y[zpivot==1],col="green",ticksize = .1,pos=-5,lwd=1.5)
  # rug(y[zpivot==2],col="red",ticksize = .1,pos=-5,lwd=1.5)
  # rug(y[zpivot==3],col="blue",ticksize = .1,pos=-5,lwd=1.5)
  
  ### Plotting stuff ###
  
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
  params = melt_sims(sims,G,num_R)
  return(list(params=params,sims=sims))
}

compute_ellipse = function(params,type,iters,lps,limit){
  # computes the hyper-ellipsoid used for the THAMES;
  # returns the mu, M, c, where all points theta in the hyper-ellipsoid fulfill 
  # (theta-mu)^TM(theta-mu)<c^2
  #browser()
  props=NULL # Only used in multimodal setting
  if(type=="standard"){
    #compute sample mean and covariance, take c_opt=sqrt(d+1)
    theta_hat <- colMeans(params[1:iters,])
    sigma_hat <- cov(params[1:iters,])
    c_opt = sqrt(length(theta_hat)+1)
  } else if(type=="glasso"){
    #browser()
    theta_hat <- colMeans(params[1:iters,])
    sigma_hat <- cov(params[1:iters,])
    sigma_hat <- glasso(sigma_hat,rho=.02)$w
    c_opt = sqrt(length(theta_hat)+1)
  } else if(type=="LW"){
    theta_hat <- colMeans(params[1:iters,])
    sigma_hat <- linearShrinkLWEst(params[1:iters,])
    c_opt = sqrt(length(theta_hat)+1)
  }else if(type=="robust"){
    #theta_hat <- params[which.max(lps[1:iters]),]
    
    theta_hat <- params[which.max(lps),] 
    # assures that the THAMES always terminates
    
    sigma_hat = 0
    for (i in 1:(2*iters)){
      sigma_hat <- sigma_hat + t(t(params[i,]-theta_hat)) %*% t(params[i,]-theta_hat)
    }
    sigma_hat = (1/(2*iters))*sigma_hat
    c_opt = sqrt(length(theta_hat)+1)
  } else if(type=="min_vol_60%HPD"){
    #browser()
    #HPDs = (params[1:iters,])[(lps[1:iters] > quantile(lps[1:iters],.4)),] TODO PUT BACK
    HPDs = (params)[(lps > quantile(lps,.5)),]
    #dim(HPDs)
    #HPDs = sapply(1:ncol(HPDs),function(s) (HPDs[,s]-mean(HPDs[,s])))
    means_of_HPDs = colMeans(HPDs)
    sds_of_HPDs = sqrt(diag(cov(HPDs)))
    HPDs = sapply(1:ncol(HPDs),function(s) (HPDs[,s]-mean(HPDs[,s]))/sd(HPDs[,s]))
    min_ellipse = try(mvee_REX(HPDs,alg.AA="MUL"))
    min_ellipse$H = min_ellipse$H# /1.01 # to account for inaccuracies
    params_centered = HPDs - t(matrix(rep(min_ellipse$z,nrow(HPDs)),ncol=nrow(HPDs)))
    num_inA =(rowSums((params_centered %*% (min_ellipse$H)) * params_centered))<=1#*cor3 TODO PUT BACK?
    print(num_inA)
    loop_counter = 0
    while(is.character(min_ellipse[1]) && loop_counter < 100){
      min_ellipse = try(mvee_REX(HPDs))
      loop_counter = loop_counter + 1
    } #Sometimes, the algorithm fails to converge
    theta_hat = c(min_ellipse$z%*%diag(sds_of_HPDs) + means_of_HPDs)
    sigma_hat = diag(sds_of_HPDs)%*%solve(min_ellipse$H)%*%diag(sds_of_HPDs)
    #params_centered = params - t(matrix(rep(theta_hat,2*iters),ncol=2*iters))
    #num_inA =(rowSums((params_centered %*% solve(sigma_hat)) * params_centered))<=1#*cor3 TODO PUT BACK?
    c_opt = 1 # it overestimates a bit
  } else if(type=="max_mode_mclust"){
    #browser()
    ICL <- mclustICL(params[1:iters,][(lps>limit)[1:iters],],modelNames="VVV")
    plot(ICL)
    fitted_K_mclust = which.max(ICL)
    #fitted_K_mclust = which.min(ICL-c(ICL[-1],-Inf))+1
    mod4 <- densityMclust(params[1:iters,][(lps>limit)[1:iters],],G=fitted_K_mclust,modelNames="VVV",plot=FALSE)
    
    counter=0
    while(is.null(mod4)&(counter<10)){
      mod4 <- densityMclust(params[1:iters,][(lps>limit)[1:iters],],G=fitted_K_mclust,modelNames="VVV",plot=FALSE)
      counter = counter + 1
    }
    theta_hat = mod4$parameters$mean
    #browser()
    closest_param = function(s){
      len = dim(params[(iters+1):(2*iters),][(lps>limit)[(iters+1):(2*iters)],])[1]
      centered_params = params[(iters+1):(2*iters),][(lps>limit)[(iters+1):(2*iters)],] - t(matrix(rep(mod4$parameters$mean[,s],len),nrow=ncol(params)))
      return(params[(iters+1):(2*iters),][(lps>limit)[(iters+1):(2*iters)],][which.min(diag(centered_params %*% t(centered_params))),])
    }
    theta_hat = sapply(1:fitted_K_mclust, closest_param)
    
    sigma_hat = mod4$parameters$variance$sigma
    props = mod4$parameters$pro
    #fitted_K_mclust = which.min(ICL-c(ICL[-1],-Inf))+1
    # mod4 <- densityMclust(params[1:iters,],G=fitted_K_mclust,modelNames="VVV",plot=FALSE)
    # max_prop_mode = which.max(mod4$parameters$pro)
    # theta_hat = mod4$parameters$mean[,max_prop_mode]
    # sigma_hat = mod4$parameters$variance$sigma[,,max_prop_mode]
    c_opt = sqrt(nrow(theta_hat)+1)
  }
  #browser()
  ellipse = list(theta_hat=theta_hat,
                 sigma_hat=sigma_hat,
                 c_opt=c_opt,
                 props=props)
  return(ellipse)
}

# computes alpha based on the Kolmogorov-Smirnov test statistic
chisq_find_limit = function(lps,d_par){
  ### add the chi-squared correction ###
  #browser()
  #hist(lps_marginal)
  #lps = lps[(iters+1):(2*iters)]
  # compute minimum Kolmogorov distance
  mean_chisq = d_par#length(lps)
  sd_chisq = sqrt(2*length(lps))
  neg_lps = -lps#[(iters+1):(2*iters)]
  #browser()
  # trunc_i = function(i,plot_res=FALSE){
  #   #browser()
  #   limit_i = sort(neg_lps,decreasing=TRUE)[i]
  #   neg_lps_trunc = neg_lps[neg_lps<=limit_i]
  #   neg_lps_trunc_mu_hat = mean(neg_lps_trunc)
  #   neg_lps_trunc_sigma_hat = sd(neg_lps_trunc)
  #   
  #   neg_lps_trunc_mu_tilde = neg_lps_trunc_mu_hat - mean_chisq
  #   neg_lps_trunc_sigma_tilde = neg_lps_trunc_sigma_hat / sd_chisq
  #   
  #   neg_lps_norm_chisq = ((neg_lps_trunc-neg_lps_trunc_mu_hat)/neg_lps_trunc_sigma_tilde)+mean_chisq
  #   neg_lps_norm_chisq = neg_lps_norm_chisq + runif(1,min=-1,max=1)
  #   #cdf_chisq = Vectorize(function(t) pchisq(t/neg_lps_trunc_sigma_tilde-neg_lps_trunc_mu_tilde,df=length(neg_lps)))
  #   cdf_chisq = Vectorize(function(t) pchisq(t, df=mean_chisq)/pchisq(max(neg_lps_norm_chisq), df=mean_chisq))
  #   ecdf_func = Vectorize(function(x) sum(neg_lps_norm_chisq<=x)/length(neg_lps_norm_chisq))
  #   if(plot_res){
  #     percent_non_zero = round(length(neg_lps_norm_chisq)/length(lps),3)
  #     plot(ecdf_func(sort(neg_lps_norm_chisq)),type="l",col="purple",
  #          xlab="index of tau",ylab="CDF",main=sprintf("Truncation: %f",percent_non_zero))
  #     
  #     lines(cdf_chisq(sort(neg_lps_norm_chisq)),type="l",col="green3")
  #     legend(1, 1, legend=c("pchisq", "ECDF"),
  #            col=c("green3", "purple"), lty=1:2, cex=0.8)
  #   }
  #   #print(i)
  #   max(abs(ecdf_func(neg_lps_norm_chisq) - cdf_chisq(neg_lps_norm_chisq)))
  # }
  # thinned_sequence = seq(2,length(neg_lps)-1,by=100) # thin to save time
  # #Vectorize(trunc_i)(thinned_sequence)
  # kolm_dist = sapply(thinned_sequence,function(i) trunc_i(i))
  # kolm_dist = sapply(thinned_sequence,function(i) trunc_i(i))
  #browser()
  thinned_sequence = seq(2,floor((length(neg_lps)-1)*0.8),by=100)#seq(2,length(neg_lps)-1,by=100)#3:(length(neg_lps)-1) # thin to save time
  # thinned_sequence = seq(2,length(neg_lps)-1,by=100) # thin to save time
  limit_i = sort(neg_lps,decreasing=TRUE)[thinned_sequence]
  neg_lps_trunc = t(matrix(rep(neg_lps,length(limit_i)),ncol=length(limit_i)))
  
  # needed to compute the rowSums
  sum_dummy = neg_lps_trunc * ((t(matrix(rep(neg_lps,length(limit_i)),ncol=length(limit_i))) <= 
                                  t(matrix(rep(limit_i,each=length(neg_lps)),ncol=length(limit_i)))))
  # these are the values that are not in the intersection with H_alpha
  neg_lps_trunc[!(t(matrix(rep(neg_lps,length(limit_i)),ncol=length(limit_i))) <= 
                    t(matrix(rep(limit_i,each=length(neg_lps)),ncol=length(limit_i))))] = NA
  
  #neg_lps_trunc = neg_lps[neg_lps<=limit_i]
  neg_lps_trunc_mu_hat = rowSums(sum_dummy)/rowSums(!is.na(neg_lps_trunc))
  neg_lps_trunc_sigma_hat = sqrt(rowSums((!is.na(neg_lps_trunc))*(sum_dummy - neg_lps_trunc_mu_hat)^2)/
                                   (rowSums(!is.na(neg_lps_trunc)) - 1))
  neg_lps_trunc_sigma_hat[neg_lps_trunc_sigma_hat==0] = 1e-6 # sd can't be 0
  neg_lps_trunc_mu_tilde = (neg_lps_trunc_mu_hat - mean_chisq)*(!is.na(neg_lps_trunc))
  neg_lps_trunc_sigma_tilde = (neg_lps_trunc_sigma_hat / sd_chisq) 
  neg_lps_norm_chisq = (((neg_lps_trunc-neg_lps_trunc_mu_hat)/neg_lps_trunc_sigma_tilde)+mean_chisq)#*(neg_lps_trunc > 0)
  
  sort_non_na = function(x){
    x[!is.na(x)] = rank(x[!is.na(x)])
    return(x)
  }
  
  # jitter the sequence so that there are no repetitions
  #neg_lps_norm_chisq = neg_lps_norm_chisq + matrix(runif(length(neg_lps_norm_chisq),min=-1e-10,max=1e-10),ncol=ncol(neg_lps_norm_chisq))
  
  #cdf_chisq_mat = matrix(pchisq(c(neg_lps_norm_chisq), df=mean_chisq),ncol=ncol(neg_lps_norm_chisq))
  cdf_chisq_mat = matrix(pchisq(c(neg_lps_norm_chisq), df=mean_chisq)/pchisq(max(neg_lps_norm_chisq[!is.na(neg_lps_norm_chisq)]), df=mean_chisq),ncol=ncol(neg_lps_norm_chisq))
  ecdf_mat = t(apply(neg_lps_norm_chisq,1,function(x) sort_non_na(x)))/rowSums(!is.na(neg_lps_trunc))
  kolm_dist = apply(abs(cdf_chisq_mat-ecdf_mat),1,function(x) max(x[!is.na(x)]))
  
  length(thinned_sequence)
  opt_index = thinned_sequence[which.min(kolm_dist)]
  #trunc_i(opt_index, plot_res=TRUE)
  #browser()
  # this is the truncation parameter that minimizes the Kolmogorov distance
  
  #limit = -quantile(neg_lps,0.99)
  limit = -sort(neg_lps,decreasing=TRUE)[opt_index]
  return(limit)
}

extend_param = function(param_test, G){
  if((ncol(param_test)/G) > 1){
    if(G==2){
      param_test_extended = cbind(param_test,1-param_test[,ncol(param_test)])
    } else{
      param_test_extended = cbind(param_test,1-rowSums(param_test[,(ncol(param_test)-(G-2)):ncol(param_test)]))
    }
  } else{
    param_test_extended = param_test
  }
  return(param_test_extended)
}

compute_thames = function(ellipse,params,lps,G,iters,type,logpost,num_R,seed=2024,prior_sampler=NULL,y=NULL,limit=NULL,num_var_g=3,sims=NULL){
  # Computes the THAMES
  # type: "simple" is the efficient version of the mixture THAMES
  #       "standard" is the vanilla THAMES
  #       "permutations" is the mixtures of G! THAMES, each concentrated on 
  #       a different posterior mode
  #browser()
  
  if(type=="mc"){
    #browser()
    naive_mc_likelihood=prior_sampler(y=y, G=G, iters=2*iters)
    naive_mc_point_est = log(mean(exp(naive_mc_likelihood - max(naive_mc_likelihood)))) + max(naive_mc_likelihood)
    naive_mc_log_se = log(sd(exp(naive_mc_likelihood - max(naive_mc_likelihood)))/sqrt(2*iters)) + max(naive_mc_likelihood)
    naive_mc_cv = exp(naive_mc_log_se - naive_mc_point_est)
    
    trunc_quantile = function(p,ratio){
      alpha = - 1/ratio
      (qnorm(p=pnorm(alpha)+p*(1-pnorm(alpha))) * (ratio))+1
    } #Calculates quantiles of the truncated normal
    
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
  }
  
  if(type=="bridge"){
    #browser()
    if(ncol(params)/G>1){
      logpost_bridge = function(param.row,data) logpost(rbind(c(param.row,1-sum(param.row[(ncol(params)-(G-1)+1):ncol(params)])),
                                                              c(param.row,1-sum(param.row[(ncol(params)-(G-1)+1):ncol(params)]))),G)[1]
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
  }
  #browser()
  theta_hat = ellipse$theta_hat
  sigma_hat = ellipse$sigma_hat
  c_opt = ellipse$c_opt
  d_par = length(theta_hat)
  graph=NULL

  n_simuls = 2*iters# 10000
  
  #browser()
  
  #simsmat <- apply(sims, 3, c)
  simsmat = matrix(c(sims),nrow=dim(sims)[1])
  #sims_test = array(c(simsmat),dim=dim(sims))
  # TODO should I compute the ellipse directly here?
  #browser()
  if(dim(sims)[3]==1){
    mu_post = colMeans(simsmat[1:iters,])
    sigma_post = cov(simsmat[1:iters,])
  } else{
    mu_post = colMeans(simsmat[1:iters,-ncol(simsmat)])
    sigma_post = cov(simsmat[1:iters,-ncol(simsmat)])
  }
  
  # sigma_post = sigma_hat
  # mu_post = theta_hat
  
  inv_post_var = sparsediscrim::solve_chol(sigma_post)
  #browser()
  set.seed(seed)
  counter = 0
  log_cor = -Inf
  
  in_ellipse = TRUE
  # sometimes, the ellipse contains no point within itself
  # TODO: fix this
  # browser()
  c_opt_old = c_opt
  while(is.infinite(log_cor) & in_ellipse){
    #browser()
    c_opt = c_opt_old / 2^counter
    print(c_opt)
    param_test = runif_in_ellipsoid(n_simuls, inv_post_var, c_opt) + 
      t(matrix(rep(mu_post,n_simuls),ncol=n_simuls))
    param_test_extended = extend_param(param_test, G)
    sims_test = array(c(param_test_extended),dim=dim(sims))
    param_test = melt_sims(sims_test,G,num_R)[,-(ncol(param_test)+1)]
    param_test_extended = extend_param(param_test, G)
    
    lps_test = logpost(param_test_extended,G)
    
    log_cor = log(mean(lps_test>limit))
    counter = counter + 1
    
    #browser()
    params_centered = params[(iters+1):(2*iters),] - t(matrix(rep(theta_hat,iters),ncol=iters))
    radius = c_opt
    
    # check if ellipse empty; set mean to mode of second half if not
    in_E = rowSums((params_centered %*% sparsediscrim::solve_chol(ellipse$sigma_hat)) * params_centered) <=radius^2
    empty_ellipse =  (sum(in_E) == 0)
    if(empty_ellipse){
      theta_hat = params[(iters+1):(2*iters),][which.max(lps[(iters+1):(2*iters)]),]
      mu_post = theta_hat
      ellipse$theta_hat = theta_hat
      c_opt = sqrt(ncol(params)+1)
      log_cor = -Inf
      counter = 0
    }
    
    if(type=="simple"){
      #browser()
      if(c_opt_old == c_opt){
        scaling = get_lda_scaling(G,sims[(iters+1):(2*iters),,])
      }
      #browser()
      graph_and_non_I_set = calc_non_I_set(scaling, G, sims_test, num_R, c_opt)
      scaling$non_I_set = graph_and_non_I_set$non_I_set
      if(c_opt_old == c_opt){
        graph = graph_and_non_I_set$graph
      }
      #print(graph_and_non_I_set$complexity_limit_estim)
      #browser()
      
      # the set has to include at least two components
      if((length(scaling$non_I_set)==(G-1))){
        graph_and_non_I_set$complexity_limit_estim = factorial(G)
        scaling$non_I_set = scaling$non_I_set[-1]
      }
      
      if((graph_and_non_I_set$complexity_limit_estim>1000000)){
        #browser()
        log_cor = -Inf
      }
    }
    
  }
  #browser()
  #plot(5000:n_simuls,log(cumsum(lps_test>limit)[(5000:n_simuls)]/(5000:n_simuls)))
  # log_cor = log(mean(apply(param_test, 1,function(x) bound_gmm(x,G=G))))
  # print(log_cor)
  # inverse_joint_test = numeric(dim(param_test)[1])
  # if(!(limit== -Inf)){
  #   param_test = param_test[apply(param_test, 1,function(x) bound_gmm(x,G=G) ),]
  #   inverse_joint_test = apply(param_test,1, function(x) -lp_gmm_marginal(theta=c(1-sum(x[1:(G-1)]),x),y=y,G=G,m=m,R=R))
  #   log_cor2 = log(mean(inverse_joint_test < -limit))
  #   print(log_cor2)
  #   log_cor = log_cor + log_cor2
  #   #browser()
  # } else{
  #   inverse_joint_test=0
  # }
  #browser()

  if(type=="simple"){

    theta_hat_extended = extend_param(rbind(theta_hat,theta_hat),G)[1,]
    
    mu_post_extended = extend_param(rbind(mu_post,mu_post),G)[1,]
    sims_theta_hat_extended = array(c(mu_post_extended),dim=c(1,dim(sims)[2:3]))
    # theta_hat_extended = c(theta_hat,1-sum(theta_hat[(2*G+1):(3*G-1)]))
    
    param_test_f_transform = matrix(reorder_by_lda(scaling, G, sims_test)$W,ncol=G)
    
    ### BEGIN TODO REMOVE ? ###
    #browser()
    graphmat = graph_and_non_I_set$graphmat
    delta_mat = matrix(0,nrow=nrow(graphmat),ncol=ncol(graphmat))
    for(g1 in 1:(nrow(delta_mat)-1)){
      for(g2 in (g1+1):(nrow(delta_mat))){
        delta_mat[g1,g2] = (mean(param_test_f_transform[,g1] < param_test_f_transform[,g2]) == 1)
      }
    }
    delta_mat = delta_mat * (1-graph_and_non_I_set$graphmat)

    adj_matrix = delta_mat
    adj_list <- lapply((1:nrow(adj_matrix))[rowSums(adj_matrix)>0], function(i) c(i,which(adj_matrix[i, ] == 1)))
    #browser()
    perms = alltopsorts_recursion(G, adj_list)
    gc()
    rm()
    ### END TODO REMOVE ? ###
    
    # theta_hat_f_transform = reorder_by_lda(scaling,G,sims_theta_hat_extended)$W
    # sort_indices = sort(theta_hat_f_transform, index.return=TRUE)$ix
    # #browser()
    # param_test_f_transform_sorted = param_test_f_transform[,sort_indices]
    # ranges = apply(param_test_f_transform_sorted,2,range)
    # 
    # #browser()
    # log_len_perm_estim = sum(log(sapply(sapply(1:G, function(s) ranges[2,s]<=ranges[1,s:G]),function(t) sum(!t))))
    # 
    # params_f_transform = param_test_f_transform#reorder_by_lda(params, scaling, G)
    # params_sorted = param_test_f_transform_sorted#params_f_transform[,sort_indices]
    # 
    # len_perm_estim_func = Vectorize(function(complexity_limit, return_ranges=FALSE){
    #   
    #   # TODO PUT BACK ?
    #   #param_ranges = apply(params_sorted[1:iters,],2,function(param) mean(param)+sd(param)*c(-complexity_limit,complexity_limit))
    #   
    #   param_ranges = apply(params_sorted,2,function(param) mean(param)+sd(param)*c(-complexity_limit,complexity_limit))
    #   
    #   param_ranges[1,] = sapply(1:G, function(s) max(param_ranges[1,s], ranges[1,s]))
    #   param_ranges[2,] = sapply(1:G, function(s) min(param_ranges[2,s], ranges[2,s]))
    #   # any ranges more spread out than the range of the Monte Carlo sample make no sense
    #   
    #   log_len_perm_estim = sum(log(sapply(sapply(1:G, function(s) param_ranges[2,s]<=param_ranges[1,s:G]),function(t) sum(!t))))
    #   if(return_ranges){
    #     return(param_ranges)
    #   } else{
    #     return(log_len_perm_estim)
    #   }
    # })
    # complexity_limits = (1:100)/10
    # #browser()
    # lens = len_perm_estim_func(complexity_limits)
    # if(max(lens)>10){
    #   select_max = which((lens*(lens<10)) == max((lens*(lens<10)), na.rm = TRUE))
    #   ranges = matrix(len_perm_estim_func(complexity_limits[select_max[length(select_max)]],return_ranges=TRUE),nrow=2)
    # }
    # 
    # counter=0
    # perms = list()
    # for(i in 1:G){
    #   counter = counter + 1
    #   perms[[counter]] = c(i)
    # }
    # #browser()
    # G_size = 1
    # while(length(perms[[1]])<G){
    #   perm_length_old = length(perms)
    #   counter = perm_length_old
    #   for(perm in perms){
    #     successors = (1:G)[-perm]
    #     successors = successors[!sapply(successors, function(s) sum(ranges[2,s]<ranges[1,perm])>=1)]
    #     if(length(successors)!=0){
    #       for(k in successors){
    #         counter = counter + 1
    #         perms[[counter]] = c(perm,k)
    #       }
    #     }
    #   }
    #   #browser()
    #   G_size = G_size + 1
    #   delete = sapply(perms,function(perm) length(perm)<G_size)
    #   perms[delete] = NULL # remove the old sample
    #   if(length(perms)>factorial(10)){
    #     perms = perms[1:factorial(10)] # number of perms could be too high
    #     print("WARNING: Number of permutations larger than 10 factorial!")
    #   }
    # }
    # if(length(perms)==0){
    #   browser()
    # }
    # if(length(perms[[1]])==1){
    #   perms = list(perms) # this is if the list of permutations contains only 1 element
    # }
  }
  #browser()
  
  if(type=="permutations"){
    #browser()
    # thames (all permutations)
    d_par <- length(theta_hat)
    
    lml_thames_marginal <- try(thames_mixture(
      lps = lps[(iters+1):(2*iters)],
      params = (params[(iters+1):(2*iters),1:(dim(params)[2])]),
      theta_hat = theta_hat,
      sigma_hat = sigma_hat,
      d_par = d_par,
      c_opt = c_opt,
      G = G,
      perms = permn(G),
      limit=limit,
      num_R=num_R,
      num_var_g=num_var_g
    ))
    #browser()
    
    list(log_zhat_inv_L=lml_thames_marginal$log_zhat_inv_L - log_cor,
         log_zhat_inv=lml_thames_marginal$log_zhat_inv - log_cor,
         log_zhat_inv_U=lml_thames_marginal$log_zhat_inv_U - log_cor,
         log_cor=log_cor,
         len_perms=length(permn(G)),
         alpha=mean(-lps< -limit),
         c_opt=c_opt,
         etas=param_test)
  } else if((type=="standard")){
    # thames without permutations (single mode)
    #browser()
    lml_thames_marginal <- thames_mixture(
      lps = lps[(iters+1):(2*iters)],
      params = (params[(iters+1):(2*iters),1:(dim(params)[2])]),
      theta_hat = theta_hat,
      sigma_hat = sigma_hat,
      d_par = d_par,
      c_opt = c_opt,
      G = G,
      perms = list(1:G),
      limit=limit,
      num_R=num_R,
      num_var_g=num_var_g
    )
    list(log_zhat_inv_L=lml_thames_marginal$log_zhat_inv_L - log_cor + lfactorial(G),
         log_zhat_inv=lml_thames_marginal$log_zhat_inv - log_cor + lfactorial(G),
         log_zhat_inv_U=lml_thames_marginal$log_zhat_inv_U - log_cor + lfactorial(G),
         log_cor=log_cor,
         alpha=mean(-lps< -limit),
         c_opt=c_opt,
         etas=param_test,
         len_perms=1
    )
    # lml_thames_single_marginal <- thames(
    #   lps = lps[(iters+1):(2*iters)],
    #   params = params[(iters+1):(2*iters),1:(dim(params)[2])],
    #   theta_hat = theta_hat,
    #   sigma_hat = sigma_hat,
    #   d_par = d_par,
    #   c_opt = c_opt,
    #   limit = -limit
    # )
    #browser()
    # list(log_zhat_inv_L=lml_thames_single_marginal$log_zhat_inv_L + log_cor - log(factorial(G)),
    #      log_zhat_inv=lml_thames_single_marginal$log_zhat_inv + log_cor - log(factorial(G)),
    #      log_zhat_inv_U=lml_thames_single_marginal$log_zhat_inv_U + log_cor - log(factorial(G)),
    #      log_cor=log_cor,
    #      alpha=mean(-lps< -limit),
    #      etas=param_test,
    #      len_perms=1
    #      )
  } else if(type=="simple"){
    #browser()
    d_par <- length(theta_hat)
    # lml_thames_marginal <- try(thames_mixture(
    #   lps = lps[(iters+1):(2*iters)],
    #   params = (params[(iters+1):(2*iters),1:(dim(params)[2])]),
    #   theta_hat = theta_hat,
    #   sigma_hat = sigma_hat,
    #   d_par = d_par,
    #   c_opt = c_opt,
    #   G = G,
    #   perms = permn(G),
    #   limit=limit,
    #   num_R=num_var_g
    # ))
    lml_thames_marginal <- thames_mixture_simple(
      lps = lps[(iters+1):(2*iters)],
      params = (params[(iters+1):(2*iters),1:(dim(params)[2])]),
      sims=sims[(iters+1):(2*iters),,],
      theta_hat = theta_hat,
      sigma_hat = sigma_hat,
      mu_post=mu_post,
      sigma_post=sigma_post,
      d_par = d_par,
      c_opt = c_opt,
      G = G,
      perms = perms,
      limit=limit,
      scaling = scaling,
      num_R=num_R,
      num_var_g=num_var_g
    )
    #browser()
    list(log_zhat_inv_L=lml_thames_marginal$log_zhat_inv_L - log_cor,
         log_zhat_inv=lml_thames_marginal$log_zhat_inv - log_cor,
         log_zhat_inv_U=lml_thames_marginal$log_zhat_inv_U - log_cor,
         log_cor=log_cor,
         len_perms=length(perms),
         alpha=mean(-lps< -limit),
         c_opt=c_opt,
         etas=param_test,
         graph=graph)
  }
}
