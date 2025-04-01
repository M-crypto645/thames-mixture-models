# Fit the THAMES to the galaxies, enzyme, and acidity datasets
rm(list=ls())
if(strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]!="thames_mixtures"){
  setwd("thames_mixtures")
}
#source('galaxies_funcs_squares.R')
source('functions/galaxies_funcs_squares.R')
source('functions/pipeline.R')
#source('functions/thames_gmm.R')
#source("functions/thames_gmm_funcs.R")
pacman::p_load(rstan,label.switching,combinat)
#options(mc.cores = parallel::detectCores())
# remove.packages("dplyr")
# install.packages("dplyr",dependencies=TRUE)
library(multimode)

### GALAXIES ###
seed = 2024
# Celeux et al. (2018) galaxy data results
lml_true_galaxies <- data.frame(G = 2:6, lml = c(-235.37,-226.85,-226.11,-225.83,-225.75))

# load in galaxy data
y <- multimode::galaxy/1000
n <- length(y)

# prior specification (just like in Celeux et al. (2018))
R <- diff(range(y))
m <- mean(range(y))

# unnormalized log-posterior density
logposty = function(theta,G,y) lp_gmm_marginal_transform(y,theta,G,m,R)
loglik_partial = function(sims,G) loglik_gmm_partial_transform(y,sims,G)

# pipeline settings
G_list <- 2:6 # 2 more than in Celeux et al. (2018)
iters = 5000 # T/2
RELABEL_ALG_NAMES = c("ECR", "STEPHENS","PRA", 
                      "AIC","SJW") # names from label.switching 
relabel_algs = c("ECR","STEPHENS")#relabel(relabel_algs[j],sims,n,y,R,m,G,2*iters,z_dummy = z_dummy)
# options: c("PRA", "AIC","STEPHENS","SJW")

ellipse_algs = list(standard="standard")

#ellipse_algs = list(max_mode_mclust="max_mode_mclust", standard="standard",
#                    robust="robust",min_vol="min_vol_60%HPD")  

# options: robust="robust")#,min_vol="min_vol_60%HPD")
thames_algs = list(simple="simple",permutations="permutations")#,permutations="permutations")
# options: standard_pluslogGfac="standard_pluslogGfac",simple="simple")
samplers = list(jags=function(g, iters, init, seed) sim_jags(n=n,x=y,R2=R^2,m=m,k1=g,
                                                       iters=iters, init=init, seed=seed),
                stan=function(g, iters, init, seed) sim_stan(n,y,R,m,g,iters,init,seed))
num_sims=20
#num_sims = 20

results_galaxies = thames_pipeline(num_sims, logposty, loglik_partial, G_list, iters, 
                                   relabel_algs, ellipse_algs, thames_algs,
                                   samplers,num_R=1,num_var_g = 3) 
save(results_galaxies,file=paste0('data/res_galaxies','.Rda'))
write.csv(data.frame(lapply(results_galaxies$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_galaxies.csv')

results_galaxies_1sim = thames_pipeline(1, logposty, loglik_partial, 2:15, 50000,
                                   c("ECR"), ellipse_algs, list(simple="simple"),
                                   list(jags=samplers[[1]]),num_R=1,num_var_g = 3) 

save(results_galaxies_1sim,file=paste0('data/res_galaxies_1sim','.Rda'))
load(file=paste0('data/res_galaxies_1sim','.Rda'))

write.csv(data.frame(lapply(results_galaxies_1sim$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_galaxies_1sim.csv')

df_1sim = read.csv('data/res_galaxies_1sim.csv')
df_1sim[df_1sim$estimate=="estim",]$value

graphs = results_galaxies_1sim$graphs
save(graphs,file=paste0('data/res_galaxies_1sim_graphs','.Rda'))
load(paste0('data/res_galaxies_1sim_graphs','.Rda'))
co = sapply(seq_along(graphs),function(s) sum(V(graphs[[s]])$color == "blue") - sum(V(graphs[[s]])$color == "red"))
#load(file=paste0('data/res_galaxies','.Rda'))

### GALAXIES ###

### ENZYME ###
seed = 2024

# Celeux et al. (2018) enzyme data results
lml_true_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-74.2,-74.3,-75,-76.7))

# load in enzyme data
y <- multimode::enzyme
n <- length(y)

# prior specification (just like in Celeux et al. (2018))
R <- diff(range(y))
m <- mean(range(y))

# unnormalized log-posterior density
logposty = function(theta,G,y) lp_gmm_marginal_transform(y,theta,G,m,R)
loglik_partial = function(sims,G) loglik_gmm_partial_transform(y,sims,G)

# old options
#radii = function(d) sqrt(d+1)*exp(0/d) # radius from Metodiev et al (2024)
#limits = list(function(lps) quantile(lps,.01))

results_enzyme = thames_pipeline(num_sims, logposty, loglik_partial, G_list, iters, 
                                 relabel_algs, ellipse_algs, thames_algs,
                                 samplers,num_R=1,num_var_g = 3) 
save(results_enzyme,file=paste0('data/res_enzyme','.Rda'))
write.csv(data.frame(lapply(results_enzyme$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_enzyme.csv')
### ENZYME ###

### ACIDITY ###
seed = 2024

# Celeux et al. (2018) acidity data results
lml_true_acidity <- data.frame(G = 2:6, lml = c(-199.4,-198.2,-198.3,-199.0,-200.1))

# load in acidity data
y <- multimode::acidity
n <- length(y)

# prior specification (just like in Celeux et al. (2018))
R <- diff(range(y))
m <- mean(range(y))

# unnormalized log-posterior density
logposty = function(theta,G,y) lp_gmm_marginal_transform(y,theta,G,m,R)
loglik_partial = function(sims,G) loglik_gmm_partial_transform(y,sims,G)

source('functions/galaxies_funcs_squares.R')
results_acidity = thames_pipeline(num_sims, logposty, loglik_partial, G_list, iters, 
                                   relabel_algs, ellipse_algs, thames_algs,
                                   samplers,num_R=1,num_var_g = 3) # TODO test init parameter
save(results_acidity,file=paste0('data/res_acidity','.Rda'))
write.csv(data.frame(lapply(results_acidity$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_acidity.csv')
### ACIDITY ###


##### WORKBENCH #####

### ENZYME ###
# seed = 2024
# 
# # Celeux et al. (2018) enzyme data results
# lml_true_enzyme <- data.frame(G = 2:6, lml = c(-76.5,-74.2,-74.3,-75,-76.7))
# 
# # load in enzyme data
# library(multimode)
# y <- multimode::enzyme
# n <- length(y)
# 
# # prior specification (just like in Celeux et al. (2018))
# R <- diff(range(y))
# m <- mean(range(y))
# 
# # unnormalized log-posterior density
# logposty = function(theta,G,y) lp_gmm_marginal(y,theta,G,m,R)
# loglik_partial = function(sims,G) loglik_gmm_partial(y,sims,G)
# 
# # pipeline settings
# G_list <- 6 # 2 more than in Celeux et al. (2018)
# iters = 5000 # T/2
# RELABEL_ALG_NAMES = c("ECR", "STEPHENS","PRA", 
#                       "AIC","SJW") # names from label.switching 
# relabel_algs = c("ECR","STEPHENS")#relabel(relabel_algs[j],sims,n,y,R,m,G,2*iters,z_dummy = z_dummy)
# # options: c("PRA", "AIC","STEPHENS","SJW")
# ellipse_algs = list(max_mode_mclust="max_mode_mclust", standard="standard")  
# # options: robust="robust")#,min_vol="min_vol_60%HPD")
# thames_algs = list(permutations="permutations", simple="simple")#,permutations="permutations")
# # options: standard_pluslogGfac="standard_pluslogGfac",simple="simple")
# samplers = list(#stan=function(g, iters, init, seed) sim_stan(n,y,R,m,g,iters,init,seed)),
#   jags=function(g, iters, init, seed) sim_jags(n=n,x=y,R2=R^2,m=m,k1=g,
#                                                iters=iters, init=init, seed=seed))
# num_sims=20
# results_enzyme = thames_pipeline(num_sims, logposty, loglik_partial, G_list, iters, 
#                                  relabel_algs, ellipse_algs, thames_algs,
#                                  samplers, num_R=3) # TODO test init parameter
# 
# boxplot(results_enzyme$df_output[(results_enzyme$df_output$relabalg=="ECR") & 
#                            (results_enzyme$df_output$thamesalg=="permutations") & 
#                            (results_enzyme$df_output$ellipsealg=="max_mode_mclust") &
#                            (results_enzyme$df_output$estimate =="estim"),]$value)
# abline(h=-76.7)
# 
# #save(results_enzyme,file=paste0('data/res_enzyme','.Rda'))
# write.csv(data.frame(lapply(results_enzyme$df_output,
#                             as.character), stringsAsFactors=FALSE),
#           file='data/res_enzyme.csv')
### ENZYME ###

# genuine_mode_1 = colMeans(params[1000:2000,])
# genuine_mode_2 = colMeans(params[-(1:3000),])
# 
# library(ggplot2)
# 
# 
# # Create the plot
# test=y
# mode = genuine_mode_2
# mix_dnorm = function(x) sum(mode[5:6]*dnorm(x,mean=mode[1:2],sd=mode[3:4]))
# ggplot(data = as.data.frame(test), aes(x = test)) +
#   geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) + # Histogram
#   stat_function(fun = Vectorize(mix_dnorm), color = "red", size = 1) # Density curve

# save(res,file=paste0('res_full_galaxies_jags','.Rda'))
# save(sims_list,file=paste0('sims_list_full_galaxies_jags','.Rda'))
# save(params_list_list,file=paste0('params_list_full_galaxies_jags','.Rda'))
# 
# load(file=paste0('res_full_bigG_jags','.Rda'))
# df_simple = as.data.frame(cbind(23:27,-sapply(1:3,function(t) sapply(1:5,function(s) res[s,1,1,1,1,1,1,t]))))
# names(df_simple) = c("G","upper","log_marginal_likelihood","lower")
# df_simple$G = c("23","24","25","26","27")
# #dev.off()
# # ggplot(df_simple, aes(G, y=thames)) + geom_point()+
# #   geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)
# plot_2=ggplot(df_simple[-c(1,2),], aes(G, y=log_marginal_likelihood,color="THAMES")) + geom_point(size=2.5)+
#   geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)+ylab("log marginal likelihood") + theme(text = element_text(size =16))
# 
# load(file=paste0('res_full_smallG_jags','.Rda'))
# df_simple = as.data.frame(cbind(2:6,-sapply(1:3,function(t) sapply(1:5,function(s) res[s,1,1,1,1,1,1,t])),
#                                 -sapply(1:3,function(t) sapply(1:5,function(s) res[s,1,1,1,1,2,1,t]))))
# names(df_simple) = c("G","upper","thames","lower","upper_perm","thames_perm","lower_perm")
# #dev.off()
# ggplot(df_simple, aes(G, y=thames)) + geom_point(aes(color="thames"))+
#   geom_errorbar(aes(ymin=lower, ymax=upper,color="thames"), width=.2,
#                 position=position_dodge(0.05)) + geom_point(aes(G,y=thames_perm,color="thames with permutations"))+
#   geom_errorbar(aes(ymin=lower_perm, ymax=upper_perm,color="thames with permutations"), width=.2,
#                 position=position_dodge(0.05))
# ggplot(df_simple[-1,], aes(G, y=thames)) + geom_point(aes(color="thames"))+
#   geom_errorbar(aes(ymin=lower, ymax=upper,color="thames"), width=.2,
#                 position=position_dodge(0.05)) + geom_point(aes(G,y=thames_perm,color="thames with permutations"))+
#   geom_errorbar(aes(ymin=lower_perm, ymax=upper_perm,color="thames with permutations"), width=.2,
#                 position=position_dodge(0.05))
# 
# load(file=paste0('res_full_galaxies_jags','.Rda'))
# dim(res)
# lml_true <- data.frame(G = 2:6, lml = c(-235.36,-226.93,-226.14,-225.84,-225.9))
# df=as.data.frame(cbind(-res[1:5,,,,,1,1,1:3],res[1:5,,,,,1,1,5],lml_true))
# names(df)=c("upper","log_marginal_likelihood","lower","complexity","G","lml")
# 
# plot_1=ggplot(df,aes(x=G,y=log_marginal_likelihood))+geom_point(size=2.5,aes(color=" THAMES"))+
#   geom_point(size=2.5,aes(x=G,y=lml,color="bridge"))+ylab("log marginal likelihood") + theme(text = element_text(size =16))
#library(gridExtra)
#grid.arrange(plot_1,plot_2)

# dim(sims_list[[1]])
# #For G=4
# mean(apply(z_list[[3]],1,function(z) length(unique(z))<4))
# 
# plot(-c(res[1,1,1,1,1,1,1,2], res[2,1,1,1,1,1,1,2],res[3,1,1,1,1,1,1,2]),ylim=c(8720,8735))
# vec_approx=lgamma(G_list)+lgamma(G_list-1+n)-lgamma(G_list-1)-lgamma(G_list+n)
# points(c(0,cumsum(vec_approx[2:3]))-res[1,1,1,1,1,1,1,2],col="red")
# # So far: Only works when underfitting
# load(paste0('res_full_bigG_jags','.Rda'))
# load(paste0('sims_list_full_bigG_jags','.Rda'))
# load(paste0('params_list_full_bigG_jags','.Rda'))
# plot(c(10,11,12),-c(res[1,1,1,1,1,1,1,1], res[2,1,1,1,1,1,1,1], res[3,5,1,1,1,1,1,1]),
#      xlab="log marginal likelihood", ylab="G",main="truncated THAMES")
# G_list