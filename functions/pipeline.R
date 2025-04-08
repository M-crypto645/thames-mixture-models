source("functions/thames_gmm_funcs.R")
library(tidyverse)
library(reshape2)
library(mclust)

#' @title thames mixtures pipeline
#'
#' @description Can be used for uni- and multivariate Gaussian mixtures
#' @param num_sims      [int>0]             number of iid replicates from THAMES
#' @param logpost       [function]          unnormalized logposterior density
#' @param loglik_partial[function]          partial log likelihood (for stephens)
#' @param p0hat         [function]          conditional proportion of empty components
#' @param G_list        [vector of int>0]   different values of G (groups) 
#' @param iters         [int>0]             2*iters=number of sims from posterior
#' @param relabel_algs  [list of functions] algorithms from label.switching
#' @param ellipse_algs  [list of functions] determines center, scale-matrix, and 
#'                                          radius of the ellipsoid E
#' @param samplers      [list of functions] MCMC samplers, take 'init' as input
#' @param num_R         [int>0]             dimension of the data
#' @param num_var_g     [int>0]             number of parameters per mixture component
#' @param init          [list, optional]    initial value for the MCMC sampler                                       
#' @return  output      [dataframe]all different values of the thames
#'                                 for the different settings
#'          lps         [array]    logposterior values of length 2*iters
#'          y_mat       [array]    the data          
#'          alpha       [array]    optimal truncation level for the HPD region
#'          etas        [array]    Monte Carlo sample on the ellipsoid E
#'          thetas      [array]    posterior sample of the parameter vector
#'                                 (before relabelling)
#'          I_map       [array]    MAP estimate of the allocation vector
#'          thetastars  [array]    the relabelled thetas
thames_pipeline = function(num_sims, logposty, loglik_partial, G_list, iters,
                           relabel_algs, ellipse_algs, thames_algs,
                           samplers, num_R,num_var_g, init=NULL,prior_sampler=NULL,
                           sims=NULL,params=NULL,lps=NULL,truth=NULL,y=NULL,p0hat=NULL){
  #browser()
  res=array(dim=c(num_sims,
                  length(G_list),
                  length(samplers),
                  length(relabel_algs),
                  length(ellipse_algs),
                  length(thames_algs),6))
  sampler_names = numeric(length(samplers))
  # this is only to find out the dimension of the parameters
  # sampler_output = samplers[[1]](G_list[1], iters, init, seed=1)
  # sims = sampler_output$results
  
  # these will only vary with g, since g changes the dimension
  thetas_list = list()
  thetastars_list = list()
  sims_list = list()
  I_map_list = list()
  etas_list = list()  
  graphs = list()
  
  # this is assuming there is less data points than MCMC iterations
  y_mat = array(dim=c(length(G_list), num_sims,
                      length(samplers),2*iters))
  lps_mat = array(dim=c(length(G_list), num_sims,
                        length(samplers),2*iters))
  
  for(g in seq_along(G_list)){
    for(i in (1:num_sims)){
      for(s in seq_along(samplers)){
        # g=11
        # i=4
        # s=1
        print(paste0("g: ",g))
        print(paste0("i: ",i))
        print(paste0("s: ",s))
        ### STEPS: (1) Simulate (2) relabel (3) set params (4) calculate THAMES
        #browser()
        ## (1) Simulate (already done once before the loop, hence the if clause)
        # if(!((i==1)&(s==1)&(g==1))){
        #browser()
        if(is.null(params)){
          sampler_output = samplers[[s]](G_list[g], iters, init, seed=i)
          #browser()
          # 1==mean(sapply(1:1000,function(s) mean(sampler_output$alloc_vec[s,]==sampler_output$alloc_vec[1,])))
          
          ### BEGIN TODO REMOVE ###
          
          # sims = theta_full
          # iters = 80000/2
          
          ### END TODO REMOVE ###
          
          sims = sampler_output$results
          #}
          #browser()
          truth = sampler_output$truth
          #save(truth,file=paste0('data/res_G15R5_gaussmulti_truth','.Rda'))
          
          d_plus_one = ncol(melt_sims(sims, G_list[g], num_R)) # minus one because of the pis
          if(num_var_g==1){
            d_plus_one = d_plus_one + 1 # pis are known in this case
          }
          
          ### PLOT1: first parameter (before relabelling)###
          #browser()
          ggplot(data=as.data.frame(cbind(seq_along(sims[,1,1]),sims[,1,1])),aes(x=V1,y=V2)) + geom_point() + labs(x="T",y=expression(theta[1][1]))
          #plot_list = list()
          #for()
          
          # sims[,1,]
          # sims[,2,]
          ### END PLOT1 ###
          
          if(is.null(sampler_output$alloc_vec)){
            z_dummy = matrix(rep(1,n*2*iters),ncol=n)
            # some samplers do not return a matrix of allocation vectors
            # Note: this means that the ECR algorithm cant be used for the sampler
          } else{
            z_dummy = sampler_output$alloc_vec
          }
        } else{
          #browser()
          z_dummy=NULL
          d_plus_one = ncol(params) + 1
          sampler_output = list(params=params[1:(2*iters),],lps=lps[1:(2*iters)],
                                truth=truth,y=y,name="unknown")
          sims = sims[1:(2*iters),,]
        }
        thetas = array(dim=c(num_sims,
                             length(samplers),2*iters, d_plus_one-1))
        thetastars = array(dim=dim(thetas))
        
        I_map = array(dim=c(num_sims,
                            length(samplers), length(sampler_output$alloc_vec[1,])))
        
        etas = array(dim=c(num_sims,
                           length(samplers),2*iters, d_plus_one-1))
        if(mean(is.na(y_mat))==1){
          y_mat = array(dim=c(dim(y_mat)[-length(dim(y_mat))],length(c(sampler_output$y))))
        }
        #y_mat[,,,] = y_mat[,,,seq_along(c(sampler_output$y))]
        #browser()
        # overloading the logpost function (y not needed after here)
        logpost = function(theta,G) logposty(theta,G,sampler_output$y)
        sampler_names[s] = sampler_output$name

        for(j in (1:length(relabel_algs))){
          ## (2) relabel (see relabel_algs for a list of relabelling algorithms)
          #browser()
          if(is.null(params)){
            if((relabel_algs[j]=="ECR")&(is.null(sampler_output$alloc_vec))){
              next # the ECR algorithm only works if the allocation vector was sampled
            }
            unrelabelled_params = melt_sims(sims, G_list[g],num_R)
            #browser()
            lps = logpost(unrelabelled_params,G_list[g])
            #browser()
            new_labels = relabel(lps, loglik_partial, relabel_algs[j],
                                 sims,G_list[g],2*iters,z_dummy = z_dummy, y=sampler_output$y)
            #browser()
            
            ## (3) set params (choose relabeling algorithm and ellipse)
            # note: the last parameter will be removed due to simplex-constraints
            
            relab = relabel_params(sims,new_labels,G_list[g],relabel_algs[j],num_R)
            params <- relab$params
            sims = relab$sims
            if(num_var_g>1){
              params = params[,-ncol(params)]
              unrelabelled_params = unrelabelled_params[,-ncol(unrelabelled_params)]
            } # else I assume the proportions are known
            sampler_output$lps = lps
            #browser()
          } else{
            #browser()
            params = sampler_output$params
            unrelabelled_params = params
          }
          
          ### PLOT2: first parameter (after relabelling)###
          ggplot(data=as.data.frame(cbind(seq_along(params[,1]),params[,1])),aes(x=V1,y=V2)) + geom_point() + labs(x="T",y=expression(theta[1][1]))
          #ggplot(data=as.data.frame(cbind(seq_along(params[,1]),params[,11])),aes(x=V1,y=V2)) + geom_point() + labs(x="T",y=expression(theta[1][1]))
        
          # browser()
          # sims[,1,]
          # sims[,2,]
          ### END PLOT2 ###
          skipthames = FALSE
          if(!is.null(p0hat)){
            p0hat_value = p0hat(sampler_output$y, extend_param(params,G_list[g]), G_list[g])
            if(p0hat_value > 1/(2*iters)){
              res[i,g,s,j,,,] = p0hat_value
              skipthames = TRUE
              thames_res = list(etas=NA)
            }
          }
          
          if(!skipthames){
            #browser()
            # determine the truncation level alpha (named "limit" in the code)
            limit = chisq_find_limit(sampler_output$lps,d_par=ncol(params))

            for(k in (1:length(ellipse_algs))){
              #browser()
              # k=2
              print(paste0("k: ",k))
              ellipse = try(compute_ellipse(params,ellipse_algs[k],iters,sampler_output$lps,limit))
              #browser()
              
              if(is.character(ellipse)){
                next # sometimes the minVol stuff does not work
              }
              #browser()
              
              ## (4) compute thames
              for(a in (1:length(thames_algs))){
                # a=1
                #browser()
                print(paste0("a: ",a))
                
                #library(thames)
                #thames(lps=lps,params = params)-lfactorial(G)
                #if()
                
                # browser()
                if(length(dim(ellipse$sigma_hat))==2){
                  thames_res = try(compute_thames(ellipse,params,sampler_output$lps,G_list[g],iters,
                                                  thames_algs[a],
                                                  logpost,num_R,num_var_g = num_var_g,prior_sampler=prior_sampler,
                                                  y=sampler_output$y, limit=max(c(limit,median(lps))),sims=sims))
                  #browser()
                  if(!is.null(thames_res$graph)){
                    graphs[[length(graphs)+1]] = thames_res$graph
                  }
                  #browser()
                } else{
                  #browser()
                  thames_is_infitnite = TRUE
                  # if(is.character(thames_res)){
                  #   next
                  # }
                  log_zhat_inv = -Inf
                  scale_c = rep(1,dim(ellipse$sigma_hat)[3])
                  while(thames_is_infitnite){
                    #browser()
                    thames_res = try(compute_thames(list(theta_hat=ellipse$theta_hat[,1],
                                                         sigma_hat=ellipse$sigma_hat[,,1],
                                                         c_opt=ellipse$c_opt/scale_c[1]),params,sampler_output$lps,G_list[g],iters,
                                                    thames_algs[a], logpost,num_R,prior_sampler=prior_sampler,
                                                    y=sampler_output$y, limit=limit))
                    thames_res_log_zhat_inv_components = sapply(1:dim(ellipse$sigma_hat)[3], function (s)
                      try(compute_thames(list(theta_hat=ellipse$theta_hat[,s],
                                              sigma_hat=ellipse$sigma_hat[,,s],
                                              c_opt=ellipse$c_opt/scale_c[s]),params,sampler_output$lps,G_list[g],iters,
                                         thames_algs[a], logpost,num_R,prior_sampler=prior_sampler,
                                         y=sampler_output$y, limit=limit)$log_zhat_inv))
                    thames_res_log_zhat_inv_components_radii = sapply(1:dim(ellipse$sigma_hat)[3], function (s)
                      try(compute_thames(list(theta_hat=ellipse$theta_hat[,s],
                                              sigma_hat=ellipse$sigma_hat[,,s],
                                              c_opt=ellipse$c_opt/scale_c[s]),params,sampler_output$lps,G_list[g],iters,
                                         thames_algs[a], logpost,num_R,prior_sampler=prior_sampler,
                                         y=sampler_output$y, limit=limit)$c_opt))
                    if(is.null(thames_res_log_zhat_inv_components_radii[[1]])){
                      thames_is_infitnite=FALSE
                      next
                    }
                    too_small_radius = sapply(1:dim(ellipse$sigma_hat)[3], function(s) is.na(as.numeric(thames_res_log_zhat_inv_components[s])))
                    if(sum(too_small_radius)>0){
                      #browser()
                      thames_is_infitnite = TRUE
                      scale_c[too_small_radius] = scale_c[too_small_radius] / 2
                    } else{
                      # browser()
                      # ellipse_limit = -Inf
                      
                      # all ellipses have to have the same radius 
                      # and at least 50 percent of the mass is covered
                      mode_radii = max(thames_res_log_zhat_inv_components_radii)
                      # mode_radii = as.numeric(names(which.max(table(thames_res_log_zhat_inv_components_radii))))
                      ellipse_condition_01 = abs(thames_res_log_zhat_inv_components_radii - mode_radii) <= 1e-14
                      normalized_props_01 = ellipse$props * ellipse_condition_01 / sum(ellipse$props[ellipse_condition_01])
                      ellipse_limit = sort(normalized_props_01)[cumsum(sort(normalized_props_01))>0.5][1] 
                      normalized_props = normalized_props_01
                      
                      # REMOVE ?
                      ellipse_condition_02 = normalized_props_01>=ellipse_limit
                      ellipse_condition = ellipse_condition_02
                      normalized_props = normalized_props_01 * ellipse_condition / sum(normalized_props_01[ellipse_condition])
                      
                      # ellipse_limit = sort(ellipse$props)[cumsum(sort(ellipse$props))>0.5][1] 
                      # ellipse_condition_02 = ellipse$props>=ellipse_limit
                      # ellipse_condition = ellipse_condition_02
                      # normalized_props = ellipse$props * ellipse_condition / sum(ellipse$props[ellipse_condition])
                      # choose such that at least 50 percent of the mass is covered
                      
                      
                      log_zhat_inv = log(sum(normalized_props*exp(thames_res_log_zhat_inv_components)))
                      print(thames_res_log_zhat_inv_components)
                      print(log_zhat_inv)
                      #log_zhat_inv = log(sum(ellipse$props*exp(thames_res_log_zhat_inv_components)))
                    }
                    if(!is.infinite(log_zhat_inv)){
                      thames_is_infitnite = FALSE
                    } else{
                      scale_c[which(is.infinite(thames_res_log_zhat_inv_components))] = 
                        scale_c[which(is.infinite(thames_res_log_zhat_inv_components))] * 2 # reduce radius such that each ellipse is computable
                    }
                  }
                  #browser()
                  #print(thames_res$log_zhat_inv)
                  thames_res$log_zhat_inv = log_zhat_inv
                }
                
                
                if(is.character(thames_res)){
                  next
                }
                #browser()
                res[i,g,s,j,k,a,] = c(-thames_res$log_zhat_inv_L-sampler_output$truth,
                                      -thames_res$log_zhat_inv-sampler_output$truth,
                                      -thames_res$log_zhat_inv_U-sampler_output$truth,
                                      thames_res$log_cor,
                                      thames_res$len_perms,
                                      thames_res$alpha)
              }
            }
          }
          if(j<length(relabel_algs)){
            params=NULL
          }
        }
        #browser()
        if((!is.array(sampler_output))&(!is.character(thames_res))){
          #browser()
          thetas[i,s,,] = unrelabelled_params
          thetastars[i,s,,] = params
          I_map[i,s,] = z_dummy[which.max(sampler_output$lps),]
          etas[i,s,,] = thames_res$etas
        }
        lps_mat[g,i,s,] =  sampler_output$lps
        y_mat[g,i,s,] = sampler_output$y
        params=NULL
      }
      #browser()
    }
    #browser()
    thetas_list[[g]] = thetas
    thetastars_list[[g]] = thetastars
    sims_list[[g]] = sims
    I_map_list[[g]] = I_map
    etas_list[[g]] = etas
  }
  df_output = melt(res)
  names(df_output) = c("sim", "G", "sampler",
                                "relabalg",
                                "ellipsealg","thamesalg","estimate","value")
  df_output$G = c("",sapply(as.numeric(df_output$G), function(s) G_list[s]))[-1]
  
  for(i in seq_along(samplers)){
    df_output$sampler[df_output$sampler==as.character(i)] = sampler_names[i]
  }
  
  for(i in seq_along(relabel_algs)){
    df_output$relabalg[df_output$relabalg==as.character(i)] = relabel_algs[i]  
  }
  for(i in seq_along(ellipse_algs)){
    df_output$ellipsealg[df_output$ellipsealg==as.character(i)] = ellipse_algs[i]  
  }
  for(i in seq_along(thames_algs)){
    df_output$thamesalg[df_output$thamesalg==as.character(i)] = thames_algs[i]
  }
  estims = c("upper","estim","lower","logcor","permlen","alpha","fitted_K_mclust")
  for(i in seq_along(estims)){
    df_output$estimate[df_output$estimate==as.character(i)] = estims[i]
  }
  #browser()
  return(list(output=res[,,,,,,-(4:6)],
              df_output=df_output,
              lps=lps_mat,
              alpha=res[,1,1,1,1,1,6],
              etas=etas_list,
              thetas=thetas_list,
              thetastars=thetastars_list,
              sims=sims_list,
              I_map=I_map_list,
              y_mat=y_mat,
              graphs=graphs,
              truth=truth))
}