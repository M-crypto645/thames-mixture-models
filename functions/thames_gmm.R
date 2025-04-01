library(corpcor)
library(uniformly)
library(gor)
library(igraph)
library(parallel)
#source('thames_function.R')

compute_nobile_identity = function(logZhatGminus1,p0hat_value,G,dirichlet_vec,n){
  return(lgamma(sum(dirichlet_vec))+
           lgamma(n+sum(dirichlet_vec[length(dirichlet_vec)]))-
           lgamma(sum(dirichlet_vec[-length(dirichlet_vec)]))-
           lgamma(n+sum(dirichlet_vec))+
           logZhatGminus1-
           log(p0hat_value))
}

# # evaluates W for a set of values (data) 
# # given meanhat, sigmahat, non_I_set (learned using QDA)
param_i_qda_linearized = function(g, data, sims, meanhat, sigmahat, non_I_set){
  G = length(sigmahat)
  #browser() # don't forget the non_I_set
  if((length(dim(sims))==3) & (dim(sims)[3]>1)){
    testmat = sims[,g,-dim(sims)[3]]
    postprob_mat = sapply(1:G,function(g_2) mvtnorm::dmvnorm(testmat,
                                                             mean=meanhat[,g_2],
                                                             sigma=sigmahat[[g_2]],log=TRUE))
  } else{
    testvec = array(sims,dim=c(dim(sims)[1:2],1))[,g,]
    postprob_mat = sapply(1:G,function(g_2) dnorm(testvec,
                                                  mean=meanhat[g_2],
                                                  sd=sqrt(sigmahat[[g_2]]),log=TRUE))
    browser()
  }

  if(is.matrix(postprob_mat)){
    postprob_mat[,non_I_set] = -Inf
    # normalize to deal with potential numeric issues
    maxlogrows = do.call(pmax, c(as.data.frame(postprob_mat)))
    postprob_mat = postprob_mat - maxlogrows
    
    postprob_mat_normalized = exp(postprob_mat) / rowSums(exp(postprob_mat))
    
    # taking the argmax (linearized)
    helper_mat = t(matrix(rep(1:G,nrow(postprob_mat)),nrow=G)) * (postprob_mat==0)
    g_hat = do.call(pmax, c(as.data.frame(helper_mat)))
    
    maxlogrowsnormalized = do.call(pmax, c(as.data.frame(postprob_mat_normalized)))
    Wmat_row_g = (g_hat + (1-maxlogrowsnormalized))
  } else{
    postprob_mat[non_I_set] = -Inf
    maxlogrows = max(postprob_mat)
    postprob_mat = postprob_mat - maxlogrows
    postprob_mat_normalized = exp(postprob_mat) / sum(exp(postprob_mat))
    
    g_hat = which.max(postprob_mat)
    
    maxlogrowsnormalized = max(postprob_mat_normalized)
    Wmat_row_g = (g_hat + (1-maxlogrowsnormalized))
  }
  
  return(Wmat_row_g)
  
}

calc_non_I_set = function(scaling, G, sims, num_R, c_opt){
  #browser()
  meanhat = scaling$meanhat
  sigmahat = scaling$sigmahat
  if(dim(sims)[3]>1){
    data = sapply(1:(dim(sims)[3]-1), function(r) c(sims[,,r]))
    param_test = melt_sims(sims,G,num_R)#[,-(ncol(param_test)+1)]
    simsmat = matrix(c(sims),nrow=dim(sims)[1])
    mu_post = colMeans(simsmat[1:(dim(sims)[1]/2),-ncol(simsmat)])
    sigma_post_inverse = sparsediscrim::solve_chol(cov(simsmat[1:(dim(sims)[1]/2),-ncol(simsmat)]))
  } else{
    data = matrix(c(sims[,,1]),ncol=1)
    meanhat = sapply(0:(G-1), function(g) mean(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]))
    sigmahat = sapply(0:(G-1), function(g) var(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]),simplify = FALSE)
    #browser()
    param_test = melt_sims(sims,G,num_R)#[,-(ncol(param_test)+1)]
    simsmat = matrix(c(sims),nrow=dim(sims)[1])
    mu_post = colMeans(simsmat[1:(dim(sims)[1]/2),])
    sigma_post_inverse = sparsediscrim::solve_chol(cov(simsmat[1:(dim(sims)[1]/2),]))
  }
  
  #browser()
  library(quadprog)
  
  graphmat = diag(G)
  
  # num_cores <- detectCores()
  # num_use_cores = min(c(num_cores-2,9))
  for(g1 in 1:(G-1)){
    for(g2 in (g1+1):G){
      Amat = matrix(0,nrow=nrow(sigma_post_inverse),ncol=ncol(sigma_post_inverse))
      Amat_func = function(u){
        sims_identification = array(0,dim=c(1,dim(sims)[-1]))
        sims_identification[1,g1,u] = 1
        sims_identification[1,g2,u] = -1
        component_pos = matrix(c(sims_identification),nrow=1)[1,]
        if(dim(sims)[3]==1){
          return(component_pos)    
        } else{
          return(component_pos[-length(component_pos)])    
        }
      }
      #browser()

      if(dim(sims)[3]==1){
        # Amat_eqs = t(simplify2array(mclapply(1:dim(sims)[3],Amat_func,
        #                                      mc.cores = num_use_cores)))
        Amat_eqs = t(simplify2array(lapply(1:dim(sims)[3],Amat_func)))
      } else{
        # Amat_eqs = t(simplify2array(mclapply(1:(dim(sims)[3]-1),Amat_func,
        #                                      mc.cores = num_use_cores)))
        Amat_eqs = t(simplify2array(lapply(1:(dim(sims)[3]-1),Amat_func)))
      }

      Amat[1:nrow(Amat_eqs),] = Amat_eqs
      maxnorm_ellipse = 2*try(quadprog::solve.QP(Dmat=sigma_post_inverse,
                                           dvec=sigma_post_inverse%*%mu_post,
                                           Amat = t(Amat),
                                           bvec = rep(0,ncol(Amat)),
                                           meq = ncol(Amat)))$value + 
        t(mu_post)%*%sigma_post_inverse%*%t(t(mu_post))
      graphmat[g1,g2] = 0 + (maxnorm_ellipse <= (c_opt^2))
      graphmat[g2,g1] = 0 + (maxnorm_ellipse <= (c_opt^2))
    }
  }
  #browser()
  
  diag(graphmat)=0
  graph = graph_from_adjacency_matrix(graphmat,mode="undirected")
  non_I_set = gor::build_cover_greedy(graph)$set
  V(graph)$color = rep("blue",G)
  V(graph)$color[non_I_set] = "red"
  plot(graph)
  
  # num_cores <- detectCores()
  # num_use_cores = min(c(num_cores-2,9))
  # if(num_use_cores==0){
  #   browser()
  # }
  # Wmat = simplify2array(mclapply(1:G,
  #                                 function(g) param_i_qda_linearized(g, data, sims, meanhat, sigmahat, non_I_set=non_I_set),
  #                                 mc.cores = num_use_cores))
  Wmat = simplify2array(lapply(1:G,
                               function(g) param_i_qda_linearized(g, data, sims, meanhat, sigmahat, non_I_set=non_I_set)))
  W=c(Wmat)
  
  ranges = sapply(0:(G-1),function(i) range(W[(i*dim(sims)[1]+1):((i+1)*dim(sims)[1])]))
  #ranges = unique(ranges[,which(apply(ranges,2,diff)<1)][1,])
  #browser()
  fixed_params = unique(ranges[1,][which(apply(ranges,2,diff)<1)])
  complexity_limit_estim = exp(lfactorial(G) - lfactorial(length(fixed_params)))
  
  return(list(graph=graph,
              non_I_set=non_I_set,
              complexity_limit_estim=complexity_limit_estim,
              graphmat=graphmat))
}

# get a new Gxiters matrix such that its features can be easily discriminated
reorder_by_lda = function(scaling, G, sims){
  #browser()
  meanhat = scaling$meanhat
  sigmahat = scaling$sigmahat
  non_I_set = scaling$non_I_set
  if((length(dim(sims))==3) & (dim(sims)[3]>1)){
    data = sapply(1:(dim(sims)[3]-1), function(r) c(sims[,,r]))
  } else{
    #browser()
    data = matrix(c(array(sims,dim=c(dim(sims)[1:2],1))[,,1]),ncol=1)
  }
  
  # num_cores <- detectCores()
  # num_use_cores = min(c(num_cores-2,9))
  # if(num_use_cores == 0){
  #   browser()
  # }
  # Wmat = simplify2array(mclapply(1:G,
  #                                function(g) param_i_qda_linearized(g, data, sims, meanhat, sigmahat, non_I_set=non_I_set),
  #                                mc.cores = num_use_cores))
  Wmat = simplify2array(lapply(1:G,
                               function(g) param_i_qda_linearized(g, data, sims, meanhat, sigmahat, non_I_set=non_I_set)))
  W=c(Wmat)
  plot(W)
  return(list(W=W))
}

# Return lda scaling such that we can order by the first principal component
get_lda_scaling = function(G, sims){
  #browser()
  if(length(dim(sims))==3){
    data = sapply(1:(dim(sims)[3]-1), function(r) c(sims[,,r]))
    meanhat = sapply(0:(G-1), function(g) colMeans(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]))
    sigmahat = sapply(0:(G-1), function(g) cov(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]),simplify = FALSE)
  } else{
    sims = array(c(sims),dim=c(dim(sims),1))
    data = matrix(c(sims[,,1]),ncol=1)
    meanhat = sapply(0:(G-1), function(g) mean(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]))
    sigmahat = sapply(0:(G-1), function(g) var(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]),simplify = FALSE)
  }


  ### BEGIN QDA ###
  #test_qda = qda(df_lda,grouping=rep(1:G,rep(dim(parms)[1],G)))
  # data_qda = matrix(nrow=nrow(data),ncol=ncol(data))
  
  #plot(sapply(1:(G*dim(sims)[1]),param_i_qda))
  ### END TODO REMOVE THIS ? ###
  
  scaling = list(meanhat=meanhat,sigmahat=sigmahat)
  return(scaling)
  #return(list(qda_scaling=qda_scaling,lda_scaling=lda_scaling))
}

#Calculates quantiles of the truncated normal
trunc_quantile = function(p,ratio){
  alpha = - 1/ratio
  (qnorm(p=pnorm(alpha)+p*(1-pnorm(alpha))) * (ratio))+1
}

# check if a sample lies in the ellipsoid A
inA <- function(theta,theta_hat,sigma_svd,radius){
  theta_tilde =  sigma_svd$d^(-1/2) * (t(sigma_svd$v) %*% (theta-theta_hat))
  return((sum(theta_tilde^2) < radius^2))
}


# function to calculate THAMES estimate of log marginal likelihood given:
# lps: vector of unnormalized log posterior values of length n_samples
# params: matrix of parameter posterior samples of dimension n_samples * d_par
# n_samples: integer, number of posterior samples
# d_par: integer, number of parameters
# c_opt: real number, radius to use for defining the ellipsoid A
# f_correct: boolean, correct the radius for sampling error?
# theta_hat: vector of length d_par, center of the ellipsoid A
# use_mean: boolean, use posterior mean or posterior median as theta_hat?
# sigma_hat: parameter covariance matrix of size d_par
# cov_shrink: boolean, use shrinkage estimation of posterior covariance matrix?
thames_mixture <- function(
    lps,params,perms,theta_hat,sigma_hat,d_par,G,
    n_samples = NULL, c_opt = NULL, f_correct = FALSE,limit=-Inf, num_R=1,num_var_g=3
){
  #browser()
  if(is.null(n_samples)){
    n_samples <- dim(params)[1]
  }
  if(is.null(c_opt)){
    c_opt <- sqrt(d_par+1)
  }
  
  # calculate SVD of sigma_hat
  sigma_svd = svd(sigma_hat)
  
  # calculate log(det(sigma_hat))
  log_det_sigma_hat = sum(log(sigma_svd$d))
  
  # calculate radius of A
  if(f_correct){
    chisq_cutoff <- pchisq(c_opt^2, df = d_par)
    f_cutoff <- (n_samples*d_par/(n_samples-d_par+1))*
      qf(chisq_cutoff, df1 = d_par, df2 = n_samples-d_par+1)
    radius <- sqrt(f_cutoff)
  }else{
    radius <- c_opt
  }
  
  # calculate volume of A
  logvolA = d_par*log(radius)+(d_par/2)*log(pi)+log_det_sigma_hat/2-lgamma(d_par/2+1)
  
  # which permuted samples are in A?
  # params_full <- cbind(pis,mus,sigmas,params)
  num_inA <- rep(0,n_samples)
  #browser()
  ### TODO REMOVE THIS ###
  if(d_par==2 | d_par==3){
    theta_hat_extended = theta_hat
    inv_sigma_hat_extended = solve(sigma_hat)
  } else{
    theta_hat_extended = c(theta_hat,0)
    inv_sigma_hat_extended = rbind(cbind(solve(sigma_hat),rep(0,d_par)),rep(0,d_par+1))
    #to make up for the lower dimensionality
    
    if(G==2){
      params = cbind(params,1-params[,ncol(params)])
    } else{
      params = cbind(params,1-rowSums(params[,(ncol(params)-(G-2)):ncol(params)]))
    }
    
  }
  ### END TODO REMOVE THIS ###
  #browser()
  for(l in perms){
    shift = calc_shift(l,num_var_g,num_R,G)
    #shift = c(rep(l,num_R))+rep(seq(0,(num_R-1)*G,G),each=G)
    #print(shift)
    #browser()
    theta_hat_total = theta_hat_extended[shift]
    inv_sigma_hat_total = inv_sigma_hat_extended[shift,shift]
    
    #I do not use inA because I want to use the (faster) matrix multiplication
    params_centered = params - t(matrix(rep(theta_hat_total,n_samples),ncol=n_samples))
    num_inA = num_inA + 
      (rowSums((params_centered %*% inv_sigma_hat_total) * params_centered)
       <=radius^2)
    
    ### Things I used for my slides ###
    
    # # print(sum(rowSums((params_centered %*% inv_sigma_hat_total) * params_centered)
    # #           <=radius^2))
    # params_shifted = params[,shift]
    # params_centered = params_shifted - t(matrix(rep(theta_hat_extended,n_samples),ncol=n_samples))
    # inA_shifted = (rowSums((params_centered %*% inv_sigma_hat_extended) * params_centered)
    #                <=radius^2)
    # if(sum(inA_shifted)>0){
    #   shift_extra = shift
    # }
    # 
    # if(sum(rowSums((params_centered %*% inv_sigma_hat_total) * params_centered)
    #        <=radius^2)>=12){
    #   l_extra=l
    # }
  }
  
  #browser()
  
  ### Stuff I plotted for my slides ###
  
  # i=5
  # j=6
  # plot(params[,i],params[,j],type="p",ylim=c(0,40),xlim=c(0,40),xlab="theta_5",ylab="theta_6")
  # 
  # #lines(params[,j],params[,i],type="p",col="black")
  # 
  # params_centered = params - t(matrix(rep(theta_hat_extended,n_samples),ncol=n_samples))
  # inA = (rowSums((params_centered %*% inv_sigma_hat_extended) * params_centered)
  #                <=radius^2)
  # lines(params[inA,i],params[inA,j],type="p",col="red")
  # plot(params[,i],params[,j],type="p",ylim=c(0,40),xlim=c(0,40),xlab="theta_5",ylab="theta_6")
  # params_shifted = params[,shift_extra]
  # params_centered = params_shifted - t(matrix(rep(theta_hat_extended,n_samples),ncol=n_samples))
  # inA_shifted = (rowSums((params_centered %*% inv_sigma_hat_extended) * params_centered)
  #                <=radius^2)
  # #plot(params[,i],params[,j],type="p",ylim=c(0,1),xlim=c(0,1))
  # lines(params_shifted[,i],params_shifted[,j],type="p",col="black")
  # lines(params[inA,i],params[inA,j],type="p",col="red")
  # lines(params_shifted[inA_shifted,i],params_shifted[inA_shifted,j],type="p",col="blue")
  #browser()
  # calculate zhat
  log_zhat_inv  = log(mean(exp(-(lps-max(lps)))*num_inA*(lps>limit)))-logvolA-max(lps)-lfactorial(G)  
  
  # estimate ar(1) model for lps
  lp_ar <- ar(exp(-(lps-max(lps)))*num_inA*(lps>limit), order.max=1)
  phi <- lp_ar$partialacf[1]
  
  # correct standard error for autocorrelation
  standard_error <- sd(exp(-lps+max(lps))*num_inA**(lps>limit))/
    ((1-phi)*sqrt(n_samples)*mean(exp(-lps+max(lps))*num_inA**(lps>limit)))
  
  
  #browser()
  # calculate 95% lower bound
  # log_zhat_inv_L <- log_zhat_inv + log(max(0,1 - 1.96*standard_error))
  log_zhat_inv_L <- log_zhat_inv + log(trunc_quantile(0.025,standard_error))
  
  # calculate 95% upper bound
  # log_zhat_inv_U <- log_zhat_inv + log(1 + 1.96*standard_error)
  log_zhat_inv_U <- log_zhat_inv + log(trunc_quantile(0.975,standard_error))
  
  return(list(theta_hat = theta_hat,
              sigma_hat = sigma_hat, 
              sigma_svd = sigma_svd, 
              log_det_sigma_hat = log_det_sigma_hat, 
              logvolA = logvolA, num_inA = num_inA, log_zhat_inv = log_zhat_inv,
              log_zhat_inv_L = log_zhat_inv_L,
              log_zhat_inv_U = log_zhat_inv_U,
              se = standard_error,
              phi = phi, radius = radius, c_opt = c_opt,
              d_par = d_par, G=G
  ))
}

calc_shift = function(sort_indices,num_var_g, num_R,G){
  if(num_R>1){
    #browser()
    if(num_var_g == 2*num_R+1){
      blocks = c(rep(1:G,each=num_R),
                 rep(1:G,each=num_R),1:G)
      block_mat = sapply(1:G, function(s) which(blocks == s))
      block_mat_shifted = block_mat[,sort_indices]
      shift = c(c(block_mat_shifted[1:num_R,]),
                c(block_mat_shifted[-(1:num_R),][1:num_R,]),
                c(block_mat_shifted[nrow(block_mat),]))
    } else{
      blocks = c(rep(1:G,each=num_R),
                 rep(1:G,each=num_R*(num_R-1)/2+num_R),1:G)
      block_mat = sapply(1:G, function(s) which(blocks == s))
      block_mat_shifted = block_mat[,sort_indices]
      shift = c(c(block_mat_shifted[1:num_R,]),
                c(block_mat_shifted[-(1:num_R),][1:(num_R*(num_R-1)/2+num_R),]),
                c(block_mat_shifted[nrow(block_mat),]))
    }

  } else{
    shift = c(rep(sort_indices,num_var_g))+rep(seq(0,(num_var_g-1)*G,G),each=G)
  }
  return(shift)
}

# function to calculate THAMES estimate of log marginal likelihood given:
# lps: vector of unnormalized log posterior values of length n_samples
# params: matrix of parameter posterior samples of dimension n_samples * d_par
# n_samples: integer, number of posterior samples
# d_par: integer, number of parameters
# c_opt: real number, radius to use for defining the ellipsoid A
# f_correct: boolean, correct the radius for sampling error?
# theta_hat: vector of length d_par, center of the ellipsoid A
# use_mean: boolean, use posterior mean or posterior median as theta_hat?
# sigma_hat: parameter covariance matrix of size d_par
# cov_shrink: boolean, use shrinkage estimation of posterior covariance matrix?
thames_mixture_simple <- function(
    lps,params,sims,perms,theta_hat,sigma_hat,mu_post,sigma_post,d_par,G,
    n_samples = NULL, c_opt = NULL, f_correct = FALSE,limit=-Inf,scaling=c(0,1,0),
    num_R=1, num_var_g=3
){
  #browser()
  if(is.null(n_samples)){
    n_samples <- dim(params)[1]
  }
  if(is.null(c_opt)){
    c_opt <- sqrt(d_par+1)
  }
  # calculate SVD of sigma_hat
  sigma_svd = svd(sigma_hat)
  
  # calculate log(det(sigma_hat))
  log_det_sigma_hat = sum(log(sigma_svd$d))
  
  # calculate radius of A
  if(f_correct){
    chisq_cutoff <- pchisq(c_opt^2, df = d_par)
    f_cutoff <- (n_samples*d_par/(n_samples-d_par+1))*
      qf(chisq_cutoff, df1 = d_par, df2 = n_samples-d_par+1)
    radius <- sqrt(f_cutoff)
  }else{
    radius <- c_opt
  }
  
  # calculate volume of A
  logvolA = d_par*log(radius)+(d_par/2)*log(pi)+log_det_sigma_hat/2-lgamma(d_par/2+1)
  
  # which permuted samples are in A?
  # params_full <- cbind(pis,mus,sigmas,params)
  num_inA <- rep(0,n_samples)
  
  #theta_hat_extended = c(0,theta_hat)
  #inv_sigma_hat_extended = rbind(rep(0,d_par+1),cbind(rep(0,d_par),solve(sigma_hat)))
  
  if(num_var_g>1){
    theta_hat_extended = c(theta_hat,0)
    inv_sigma_hat_extended = rbind(cbind(solve(sigma_hat),rep(0,d_par)),rep(0,d_par+1))
  } else{
    theta_hat_extended = theta_hat
    inv_sigma_hat_extended = solve(sigma_hat)
  }
  
  #to make up for the lower dimensionality
  
  #browser()
  #params = matrix(c(sims),nrow=dim(sims)[1])
  params = extend_param(params,G)
  
  # if(G==2){
  #   params = cbind(params,1-params[,ncol(params)])
  # } else{
  #   params = cbind(params,1-rowSums(params[,(ncol(params)-(G-2)):ncol(params)]))
  # }
  
  #browser()
  
  theta_hat_f = extend_param(rbind(mu_post,mu_post),G)[1,]
  sims_theta_hat_f = array(c(theta_hat_f),dim=c(1,G,num_var_g))
  theta_hat_f_transform = reorder_by_lda(scaling, G, sims_theta_hat_f)$W
  #theta_hat_f_transform = reorder_by_lda(c(1-sum(theta_hat[1:(G-1)]),theta_hat),scaling, G)
  sort_indices = sort(theta_hat_f_transform,index.return=TRUE)$ix
  #sort_indices = sort(theta_hat[G:(2*G-1)],index.return=TRUE)$ix
  #sort_indices = 1:G

  shift = calc_shift(sort_indices,num_var_g, num_R,G)
  #shift = c(rep(sort_indices,num_R))+rep(seq(0,(num_R-1)*G,G),each=G)
  theta_hat_total = theta_hat_extended[shift]
  inv_sigma_hat_total = inv_sigma_hat_extended[shift,shift]
  #to make up for the lower dimensionality

  #I do not use inA because I want to use the (faster) matrix multiplication
  sorted_params = params[,shift]
  #params_sims = array(c(params),dim=c(nrow(params),G,num_var_g))
  params_f_transform = matrix(reorder_by_lda(scaling,G,sims)$W,ncol=G)
  #browser()
  # if(length((1:G)[-id_trunc])==1){
  #   perms = c((1:G)[-id_trunc])
  # } else{
  #   perms = permn((1:G)[-id_trunc])
  # }
  # if(length(id_trunc)==G){
  #   perms=list(sort_indices)
  # }
  # 
  for(l in perms){
    #Only shift the part that is not truncated
    #browser()
    l = sort_indices[l]
    shift = calc_shift(l,num_var_g, num_R,G) #c(rep(l,num_R))+rep(seq(0,(num_R-1)*G,G),each=G)
    #print(shift)
    
    # check that order fits
    sorted_params_f_transform = params_f_transform[,sort_indices]
    
    #TODO PUT BACK?
    #cor3=rep(1,n_samples)
    # TODO
    # for(g in (1:(G-1))){
    #   cor3 = cor3 * (sorted_params_f_transform[,g]<sorted_params_f_transform[,g+1])
    # }
    #browser()
    # check if included in permuted A
    theta_hat_total = theta_hat_extended[shift]
    inv_sigma_hat_total = inv_sigma_hat_extended[shift,shift]
    
    #I do not use inA because I want to use the (faster) matrix multiplication
    params_centered = sorted_params - t(matrix(rep(theta_hat_total,n_samples),ncol=n_samples))
    num_inA = num_inA + 
      (rowSums((params_centered %*% inv_sigma_hat_total) * params_centered)
       <=radius^2)#*cor3 TODO PUT BACK?
  }
  #browser()
  # params_centered = sorted_params - t(matrix(rep(theta_hat_total,n_samples),ncol=n_samples))
  # num_inA = (rowSums((params_centered %*% inv_sigma_hat_total) * params_centered)
  #    <=radius^2)
  print(sum(num_inA))

  #browser()
  # calculate zhat
  (log_zhat_inv  = log(mean(exp(-(lps-max(lps)))*num_inA*(lps>limit)))-logvolA-max(lps)-lfactorial(G))  
  
  # estimate ar(1) model for lps
  lp_ar <- ar(exp(-(lps-max(lps)))*num_inA*(lps>limit), order.max=1)
  phi <- lp_ar$partialacf[1]
  
  # correct standard error for autocorrelation
  standard_error <- sd(exp(-lps+max(lps))*num_inA**(lps>limit))/
    ((1-phi)*sqrt(n_samples)*mean(exp(-lps+max(lps))*num_inA**(lps>limit)))
  
  
  
  # calculate 95% lower bound
  # log_zhat_inv_L <- log_zhat_inv + log(max(0,1 - 1.96*standard_error))
  log_zhat_inv_L <- log_zhat_inv + log(trunc_quantile(0.025,standard_error))
  
  # calculate 95% upper bound
  # log_zhat_inv_U <- log_zhat_inv + log(1 + 1.96*standard_error)
  log_zhat_inv_U <- log_zhat_inv + log(trunc_quantile(0.975,standard_error))
  
  return(list(theta_hat = theta_hat,
              sigma_hat = sigma_hat, 
              sigma_svd = sigma_svd, 
              log_det_sigma_hat = log_det_sigma_hat, 
              logvolA = logvolA, num_inA = num_inA, log_zhat_inv = log_zhat_inv,
              log_zhat_inv_L = log_zhat_inv_L,
              log_zhat_inv_U = log_zhat_inv_U,
              se = standard_error,
              phi = phi, radius = radius, c_opt = c_opt,
              d_par = d_par, G=G
  ))
}