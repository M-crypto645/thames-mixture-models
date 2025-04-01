rm(list=ls())
if(strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]!="thames_mixtures"){
  setwd("thames_mixtures")
}

library(mvtnorm)
library(sparsediscrim)
library(CholWishart)
library(parallel)
library(mvnfast)
library(LaplacesDemon)
library(tictoc)
library(mclust)

source("functions/pipeline.R")
source("functions/true_marglik_gaussmult_funcs.R")

### BEGIN liver dataset ###

# real values
n = 345
G = 15#3#30
R = 5#6#6#27#6#57

# specification of the parameters from which the data is simulated
dist = 100
mustars_func = function(g) (t(matrix(rep(seq(dist,dist*g,dist),R),ncol=R,nrow=g)))
sigmastars_func = function(g) simplify2array(sapply(1:g,function(s) diag(R),simplify = FALSE))
taustars_func = function(g) rep(1/g,g)

# prior hyperparameters
# (the others depend on y and are thus defined within the function)
alpha_0_constant = 1
nu = 2

iters = 10000
burn_in = 2000

source("liver_data.R")
y = liver[,-(6:7)]
co = rep(0,15)
margliks = rep(0,15)

results_T10000_G15 = thames_pipeline(num_sims=1,
                                  logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                  p0hat=p0hat_gaussmulti_vii,
                                  loglik_partial=NULL,
                                  G_list=15, iters=iters/2,
                                  relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                  thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                  num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                  prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                  samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G15 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=15, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[15] = 3 - (15-3)
save(results_T10000_G15,file=paste0('data/res_liver_G15_T10000','.Rda'))
load(file=paste0('data/res_liver_G15_T10000','.Rda'))
results_T10000_G15$df_output

results_T10000_G14 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=14, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G14 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=14, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[14] = 3 - (14-3)
save(results_T10000_G14,file=paste0('data/res_liver_G14_T10000','.Rda'))
load(file=paste0('data/res_liver_G14_T10000','.Rda'))
results_T10000_G14$df_output

results_T10000_G13 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=13, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G13 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=13, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[13] = 3 - (13-3)
save(results_T10000_G13,file=paste0('data/res_liver_G13_T10000','.Rda'))
load(file=paste0('data/res_liver_G13_T10000','.Rda'))
results_T10000_G13$df_output

results_T10000_G12 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=12, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G12 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=12, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[12] = 3 - (12-3)
save(results_T10000_G12,file=paste0('data/res_liver_G12_T10000','.Rda'))
load(file=paste0('data/res_liver_G12_T10000','.Rda'))
results_T10000_G12$df_output

results_T10000_G11 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=11, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G11 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=11, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[11] = 3 - (11-3)
save(results_T10000_G11,file=paste0('data/res_liver_G11_T10000','.Rda'))
load(file=paste0('data/res_liver_G11_T10000','.Rda'))
results_T10000_G11$df_output

results_T10000_G10 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=10, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G10 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=10, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[10] = 3 - (10-3)
save(results_T10000_G10,file=paste0('data/res_liver_G10_T10000','.Rda'))
load(file=paste0('data/res_liver_G10_T10000','.Rda'))
results_T10000_G10$df_output

results_T10000_G09 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=9, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G09 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=9, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[9] = 3 - (9-3)
save(results_T10000_G09,file=paste0('data/res_liver_G09_T10000','.Rda'))
load(file=paste0('data/res_liver_G09_T10000','.Rda'))
results_T10000_G09$df_output

results_T10000_G08 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=8, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G08 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=8, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[8] = 3 - (8-3)
save(results_T10000_G08,file=paste0('data/res_liver_G08_T10000','.Rda'))
load(file=paste0('data/res_liver_G08_T10000','.Rda'))
results_T10000_G08$df_output

results_T10000_G07 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=7, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G07 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=7, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[7] = 3 - (7-3)
save(results_T10000_G07,file=paste0('data/res_liver_G07_T10000','.Rda'))
load(file=paste0('data/res_liver_G07_T10000','.Rda'))
results_T10000_G07$df_output

results_T10000_G06 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=6, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G06 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=6, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[6] = 3 - (6-3)
save(results_T10000_G06,file=paste0('data/res_liver_G06_T10000','.Rda'))
load(file=paste0('data/res_liver_G06_T10000','.Rda'))
results_T10000_G06$df_output

results_T10000_G05 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=5, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# run to see overlap graph
# results_T10000_G05 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=5, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

co[5] = 4 - (5-4)
save(results_T10000_G05,file=paste0('data/res_liver_G05_T10000','.Rda'))
load(file=paste0('data/res_liver_G05_T10000','.Rda'))
results_T10000_G05$df_output

results_T10000_G4 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=4, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
co[4] = 4 - (4-4)
save(results_T10000_G4,file=paste0('data/res_liver_G4_T10000','.Rda'))
load(file=paste0('data/res_liver_G4_T10000','.Rda'))
results_T10000_G4$df_output

results_T10000_G03 = thames_pipeline(num_sims=1,
                                    logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                    p0hat=p0hat_gaussmulti_vii,
                                    loglik_partial=NULL,
                                    G_list=3, iters=iters/2,
                                    relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                    thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                    num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                    prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                    samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                  mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                  init,2*iters,burn_in=burn_in,seed=seed,y=y)))
co[3] = 3 - (3-3)
save(results_T10000_G03,file=paste0('data/res_liver_G03_T10000','.Rda'))
load(file=paste0('data/res_liver_G03_T10000','.Rda'))
results_T10000_G03$df_output

results_T10000_G02 = thames_pipeline(num_sims=1,
                                     logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                     p0hat=p0hat_gaussmulti_vii,
                                     loglik_partial=NULL,
                                     G_list=2, iters=iters/2,
                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                     num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                   mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                   init,2*iters,burn_in=burn_in,seed=seed,y=y)))
co[2] = 2 - (2-2)
save(results_T10000_G02,file=paste0('data/res_liver_G02_T10000','.Rda'))
load(file=paste0('data/res_liver_G02_T10000','.Rda'))
results_T10000_G02$df_output

marglik_res = rep(0,14)
marglik_res[1] = results_T10000_G02$df_output$value[results_T10000_G02$df_output$estimate=="estim"]
marglik_res[2] = results_T10000_G03$df_output$value[results_T10000_G03$df_output$estimate=="estim"]
marglik_res[3] = results_T10000_G4$df_output$value[results_T10000_G4$df_output$estimate=="estim"]
marglik_res[4] = compute_nobile_identity(logZhatGminus1 = marglik_res[3],
                                         p0hat_value = results_T10000_G05$df_output$value[1],
                                         dirichlet_vec = rep(1,5),
                                         n=nrow(y))
marglik_res[5] = compute_nobile_identity(logZhatGminus1 = marglik_res[4],
                                         p0hat_value = results_T10000_G05$df_output$value[1],
                                         dirichlet_vec = rep(1,6),
                                         n=nrow(y))
marglik_res[6] = compute_nobile_identity(logZhatGminus1 = marglik_res[5],
                                         p0hat_value = results_T10000_G05$df_output$value[1],
                                         dirichlet_vec = rep(1,7),
                                         n=nrow(y))
marglik_res[7] = compute_nobile_identity(logZhatGminus1 = marglik_res[6],
                                         p0hat_value = results_T10000_G05$df_output$value[1],
                                         dirichlet_vec = rep(1,8),
                                         n=nrow(y))
marglik_res[8] = compute_nobile_identity(logZhatGminus1 = marglik_res[7],
                                         p0hat_value = results_T10000_G05$df_output$value[1],
                                         dirichlet_vec = rep(1,9),
                                         n=nrow(y))
marglik_res[9] = compute_nobile_identity(logZhatGminus1 = marglik_res[8],
                                         p0hat_value = results_T10000_G05$df_output$value[1],
                                         dirichlet_vec = rep(1,10),
                                         n=nrow(y))
marglik_res[10] = compute_nobile_identity(logZhatGminus1 = marglik_res[9],
                                         p0hat_value = results_T10000_G05$df_output$value[1],
                                         dirichlet_vec = rep(1,11),
                                         n=nrow(y))
marglik_res[11] = compute_nobile_identity(logZhatGminus1 = marglik_res[10],
                                          p0hat_value = results_T10000_G05$df_output$value[1],
                                          dirichlet_vec = rep(1,12),
                                          n=nrow(y))
marglik_res[12] = compute_nobile_identity(logZhatGminus1 = marglik_res[11],
                                          p0hat_value = results_T10000_G05$df_output$value[1],
                                          dirichlet_vec = rep(1,13),
                                          n=nrow(y))
marglik_res[13] = compute_nobile_identity(logZhatGminus1 = marglik_res[12],
                                          p0hat_value = results_T10000_G05$df_output$value[1],
                                          dirichlet_vec = rep(1,14),
                                          n=nrow(y))
marglik_res[14] = compute_nobile_identity(logZhatGminus1 = marglik_res[13],
                                          p0hat_value = results_T10000_G05$df_output$value[1],
                                          dirichlet_vec = rep(1,15),
                                          n=nrow(y))
marglik_res = matrix(marglik_res,nrow=1)
colnames(marglik_res) = c("",2:15)[-1]
write.csv(marglik_res,"data/marglik_res_liver.csv")
read.csv("data/marglik_res_liver.csv")

# co values can be read of from the graph visualization
co = matrix(co[-1],nrow=1)
write.csv(co,"data/co_liver.csv")
read.csv("data/co_liver.csv")

pints = sort(unique(liver[,6]))
confusion_mat_liver = rbind(sapply(pints, function(s) sum((results_T10000_G4$I_map[[1]][1,1,]==1)&(liver[,6]==s))),
                            sapply(pints, function(s) sum((results_T10000_G4$I_map[[1]][1,1,]==2)&(liver[,6]==s))),
                            sapply(pints, function(s) sum((results_T10000_G4$I_map[[1]][1,1,]==3)&(liver[,6]==s))),
                            sapply(pints, function(s) sum((results_T10000_G4$I_map[[1]][1,1,]==4)&(liver[,6]==s))))
colnames(confusion_mat_liver) = c("",pints)[-1]
rownames(confusion_mat_liver) = c("1","2","3","4")

write.csv(confusion_mat_liver,"data/confusion_mat_liver.csv")
read.csv("data/confusion_mat_liver.csv")

liver_post_mean_vec = colMeans(results_T10000_G4$thetastars[[1]][1,1,,1:20])
liver_post_mean = as.data.frame(rbind(liver_post_mean_vec[1:5],
                                      liver_post_mean_vec[6:10],
                                      liver_post_mean_vec[11:15],
                                      liver_post_mean_vec[16:20]))
colnames(liver_post_mean) = c("mcv","alkphos","sgpt","sgot","gammagt")
# mcv     corpuscular volume
# alkphos alkaline phosphotase
# sgpt    alanine aminotransferase
# sgot    aspartate aminotransferase
# gammagt gamma-glutamyl transpeptidase
write.csv(liver_post_mean,"data/liver_post_mean.csv")
read.csv("data/liver_post_mean.csv")


# run to see overlap graph
# results_T10000_G15 = thames_pipeline(num_sims=1,
#                                      logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                      loglik_partial=NULL,
#                                      G_list=G, iters=iters/2,
#                                      relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                      thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                      num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
#                                      prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                      samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
#                                                                                                                    mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                    init,2*iters,burn_in=burn_in,seed=seed,y=y)))

# co[15] = 6 - 9
# save(results_T10000_G15,file=paste0('data/res_liver_G15_T10000','.Rda'))
# load(file=paste0('data/res_liver_G15_T10000','.Rda'))
# results_T10000_G15$df_output

### END liver dataset ###

### BEGIN SIMULATION STUDY G=15, R=5 ###

# real values
n = 345
G = 15#3#30
R = 5#6#6#27#6#57

# specification of the parameters from which the data is simulated
dist = 100
mustars_func = function(g) (t(matrix(rep(seq(dist,dist*g,dist),R),ncol=R,nrow=g)))
sigmastars_func = function(g) simplify2array(sapply(1:g,function(s) diag(R),simplify = FALSE))
taustars_func = function(g) rep(1/g,g)

# prior hyperparameters
# (the others depend on y and are thus defined within the function)
alpha_0_constant = 1
nu = 2

iters = 200000
burn_in = 2000
results_T200000 = thames_pipeline(num_sims=1,
                                  logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                  p0hat=p0hat_gaussmulti_vii,
                                  loglik_partial=NULL,
                                  G_list=G, iters=iters/2,
                                  relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                  thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                  num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                  prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                  samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                            mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                            init,2*iters,burn_in=burn_in,seed=seed)))
save(results_T200000,file=paste0('data/res_G15R5_gaussmulti_T200000','.Rda'))
load(file=paste0('data/res_G15R5_gaussmulti_T200000','.Rda'))
#load(file=paste0('data/res_G15R5_gaussmulti_T200000_truth','.Rda'))
#load(file=paste0('data/res_G15R5_gaussmulti_T200000_sims','.Rda'))
results_T200000$df_output

y=matrix(results_T200000$y_mat[1,1,1,],ncol=R)
params = results_T200000$thetastars[[1]][1,1,,]
lps = results_T200000$lps[1,1,1,]
truth = results_T200000$truth
sims = results_T200000$sims[[1]]
# truth = log_true_marglik_gauss_multi_vii(y=matrix(results_T200000$y_mat[1,1,1,],ncol=R), nu=nu,
#                                     lambda=calc_lambda_vii(matrix(results_T200000$y_mat[1,1,1,],ncol=R),nu),
#                                     alpha_0=rep(alpha_0_constant,G), kappa_0=calc_kappa_0(matrix(results_T200000$y_mat[1,1,1,],ncol=R)),
#                                     beta=calc_beta(matrix(results_T200000$y_mat[1,1,1,],ncol=R)), mustars = mustars_func(G))
#sims=? TODO

iters = 100000
results_T100000 = thames_pipeline(num_sims=1,
                                  logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                  p0hat=p0hat_gaussmulti_vii,
                                  loglik_partial=NULL,
                                  G_list=G, iters=iters/2,
                                  relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                  thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                  num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                  prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                  samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                                mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                init,2*iters,burn_in=burn_in,seed=seed)),
                                  params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T100000,file=paste0('data/res_G15R5_gaussmulti_T100000','.Rda'))
# load(file=paste0('data/res_G15R5_gaussmulti_T100000','.Rda'))
# results_T100000$df_output

iters = 50000
results_T50000 = thames_pipeline(num_sims=1,
                                 logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                 p0hat=p0hat_gaussmulti_vii,
                                 loglik_partial=NULL,
                                 G_list=G, iters=iters/2,
                                 relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                 thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                 num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                 prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                 samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                               mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                               init,2*iters,burn_in=burn_in,seed=seed)),
                                 params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T50000,file=paste0('data/res_G15R5_gaussmulti_T50000','.Rda'))
# load(file=paste0('data/res_G15R5_gaussmulti_T50000','.Rda'))
# results_T50000$df_output

iters = 20000
results_T20000 = thames_pipeline(num_sims=1,
                                 logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                 p0hat=p0hat_gaussmulti_vii,
                                 loglik_partial=NULL,
                                 G_list=G, iters=iters/2,
                                 relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                 thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                 num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                 prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                 samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                               mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                               init,2*iters,burn_in=burn_in,seed=seed)),
                                 params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T20000,file=paste0('data/res_G15R5_gaussmulti_T20000','.Rda'))
# load(file=paste0('data/res_G15R5_gaussmulti_T20000','.Rda'))
# results_T20000$df_output

iters = 10000
results_T10000 = thames_pipeline(num_sims=1,
                                 logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                 p0hat=p0hat_gaussmulti_vii,
                                 loglik_partial=NULL,
                                 G_list=G, iters=iters/2,
                                 relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                 thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                 num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                 prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                 samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                               mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                               init,2*iters,burn_in=burn_in,seed=seed)),
                                 params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T10000,file=paste0('data/res_G15R5_gaussmulti_T10000','.Rda'))
# load(file=paste0('data/res_G15R5_gaussmulti_T10000','.Rda'))
# results_T10000$df_output

iters = 5000
results_T5000 = thames_pipeline(num_sims=1,
                                logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                p0hat=p0hat_gaussmulti_vii,
                                loglik_partial=NULL,
                                G_list=G, iters=iters/2,
                                relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                              mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                              init,2*iters,burn_in=burn_in,seed=seed)),
                                params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T5000,file=paste0('data/res_G15R5_gaussmulti_T5000','.Rda'))
# load(file=paste0('data/res_G15R5_gaussmulti_T5000','.Rda'))
# results_T5000$df_output

iters = 2000
results_T2000 = thames_pipeline(num_sims=1,
                                logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                p0hat=p0hat_gaussmulti_vii,
                                loglik_partial=NULL,
                                G_list=G, iters=iters/2,
                                relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                              mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                              init,2*iters,burn_in=burn_in,seed=seed)),
                                params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T2000,file=paste0('data/res_G15R5_gaussmulti_T2000','.Rda'))
# load(file=paste0('data/res_G15R5_gaussmulti_T2000','.Rda'))
# results_T2000$df_output

iters = 1000
results_T1000 = thames_pipeline(num_sims=1,
                                logposty=function(thetas, G, y) logposty_gaussmulti_vii(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                p0hat=p0hat_gaussmulti_vii,
                                loglik_partial=NULL,
                                G_list=G, iters=iters/2,
                                relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                num_R=R,num_var_g=1+2*R,init=function(y) 1, # init is set within the function
                                prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_vii(n,nu,rep(alpha_0_constant,g),
                                                                                                              mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                              init,2*iters,burn_in=burn_in,seed=seed)),
                                params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T1000,file=paste0('data/res_G15R5_gaussmulti_T1000','.Rda'))
# load(file=paste0('data/res_G15R5_gaussmulti_T1000','.Rda'))
# results_T1000$df_output

results_T200000$df_output$T = 200000
results_T100000$df_output$T = 100000
results_T50000$df_output$T = 50000
results_T20000$df_output$T = 20000
results_T10000$df_output$T = 10000
results_T5000$df_output$T = 5000
results_T2000$df_output$T = 2000
results_T1000$df_output$T = 1000

df_output = rbind(results_T1000$df_output, 
                  results_T2000$df_output,
                  results_T5000$df_output,
                  results_T10000$df_output,
                  results_T20000$df_output,
                  results_T50000$df_output,
                  results_T100000$df_output,
                  results_T200000$df_output)

write.csv(data.frame(lapply(df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_truemarglikgaussmultiG15R5.csv')

### END SIMULATION STUDY G=15, R=5 ###

### BEGIN SIMULATION STUDY G=5, R=6 ###

n = 200
G = 5
R = 6

# specification of the parameters from which the data is simulated
dist = 100
mustars_func = function(g) (t(matrix(rep(seq(dist,dist*g,dist),R),ncol=R,nrow=g)))
sigmastars_func = function(g) simplify2array(sapply(1:g,function(s) diag(R),simplify = FALSE))
taustars_func = function(g) rep(1/g,g)

# prior hyperparameters
# (the others depend on y and are thus defined within the function)
alpha_0_constant = 1
nu = R

iters = 200000
burn_in = 2000
results_T200000 = thames_pipeline(num_sims=1,
                                  logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                  p0hat=p0hat_gaussmulti,
                                  loglik_partial=NULL,
                                  G_list=G, iters=iters/2,
                                  relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                  thames_algs = c("simple","bridge"),
                                  num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                  prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                  samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                            mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                            init,2*iters,burn_in=burn_in,seed=seed)))

# save(results_T200000,file=paste0('data/res_G5R6_gaussmulti_T200000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T200000','.Rda'))
# results_T200000$df_output

params = results_T200000$thetastars[[1]][1,1,,]
sims = results_T200000$sims[[1]]


lps = results_T200000$lps[1,1,1,]
truth = log_true_marglik_gauss_multi(y=matrix(results_T200000$y_mat[1,1,1,],ncol=6), nu=nu,
                                     lambda=calc_lambda(matrix(results_T200000$y_mat[1,1,1,],ncol=6),nu,G),
                                     alpha_0=rep(1,5), kappa_0=calc_kappa_0(matrix(results_T200000$y_mat[1,1,1,],ncol=6)),
                                     beta=calc_beta(matrix(results_T200000$y_mat[1,1,1,],ncol=6)), mustars = mustars_func(5))

y=matrix(results_T200000$y_mat[1,1,1,],ncol=6)

iters = 100000
results_T100000 = thames_pipeline(num_sims=1,
                          logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                          p0hat=p0hat_gaussmulti,
                          loglik_partial=NULL,
                          G_list=G, iters=iters/2,
                          relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                          thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                          num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                          prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                          samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                          mustars_func(g),sigmastars_func(g),taustars_func(g),
                                          init,2*iters,burn_in=burn_in,seed=seed)),
                          params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T100000,file=paste0('data/res_G5R6_gaussmulti_T100000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T100000','.Rda'))
# results_T100000$df_output

iters = 50000
results_T50000 = thames_pipeline(num_sims=1,
                                 logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                 p0hat=p0hat_gaussmulti,
                                 loglik_partial=NULL,
                                 G_list=G, iters=iters/2,
                                 relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                 thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                 num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                 prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                 samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                                       mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                       init,2*iters,burn_in=burn_in,seed=seed)),
                                 params=params,lps=lps,truth=truth,y=y,sims=sims)
# save(results_T50000,file=paste0('data/res_G5R6_gaussmulti_T50000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T50000','.Rda'))
# results_T50000$df_output

iters = 20000
results_T20000 = thames_pipeline(num_sims=1,
                                 logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                 p0hat=p0hat_gaussmulti,
                                 loglik_partial=NULL,
                                 G_list=G, iters=iters/2,
                                 relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                 thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                 num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                 prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                 samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                           mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                           init,2*iters,burn_in=burn_in,seed=seed)),
                                 params=params,lps=lps,truth=truth,y=y,sims=sims)

# save(results_T20000,file=paste0('data/res_G5R6_gaussmulti_T20000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T20000','.Rda'))
# results_T20000$df_output

iters = 10000
results_T10000 = thames_pipeline(num_sims=1,
                                 logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                 p0hat=p0hat_gaussmulti,
                                 loglik_partial=NULL,
                                 G_list=G, iters=iters/2,
                                 relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                 thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                 num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                 prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                 samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                           mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                           init,2*iters,burn_in=burn_in,seed=seed)),
                                 params=params,lps=lps,truth=truth,y=y,sims=sims)

# save(results_T10000,file=paste0('data/res_G5R6_gaussmulti_T10000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T10000','.Rda'))
# results_T10000$df_output

iters = 5000
results_T5000 = thames_pipeline(num_sims=1,
                                logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                p0hat=p0hat_gaussmulti,
                                loglik_partial=NULL,
                                G_list=G, iters=iters/2,
                                relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                thames_algs = c("simple","bridge"),
                                num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, 
                                prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                                      mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                      init,2*iters,burn_in=burn_in,seed=seed)),
                                params=params,lps=lps,truth=truth,y=y,sims=sims)

# save(results_T5000,file=paste0('data/res_G5R6_gaussmulti_T5000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T5000','.Rda'))
# results_T5000$df_output

iters = 2000
results_T2000 = thames_pipeline(num_sims=1,
                                logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                p0hat=p0hat_gaussmulti,
                                loglik_partial=NULL,
                                G_list=G, iters=iters/2,
                                relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                          mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                          init,2*iters,burn_in=burn_in,seed=seed)),
                                params=params,lps=lps,truth=truth,y=y,sims=sims)

# save(results_T2000,file=paste0('data/res_G5R6_gaussmulti_T2000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T2000','.Rda'))
# results_T2000$df_output

iters = 1000
results_T1000 = thames_pipeline(num_sims=1,
                                logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                p0hat=p0hat_gaussmulti,
                                loglik_partial=NULL,
                                G_list=G, iters=iters/2,
                                relabel_algs = c("ECR"), ellipse_algs = c("standard"),#c("min_vol_60%HPD"),
                                thames_algs = c("simple","bridge"),#c(permutations="permutations",simple="simple"),
                                num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                          mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                          init,2*iters,burn_in=burn_in,seed=seed)),
                                params=params,lps=lps,truth=truth,y=y,sims=sims)

# save(results_T1000,file=paste0('data/res_G5R6_gaussmulti_T1000','.Rda'))
# load(file=paste0('data/res_G5R6_gaussmulti_T1000','.Rda'))
# results_T1000$df_output 

results_T1000$df_output$T = 1000
results_T2000$df_output$T = 2000
results_T5000$df_output$T = 5000
results_T10000$df_output$T = 10000
results_T20000$df_output$T = 20000
results_T50000$df_output$T = 50000
results_T100000$df_output$T = 100000
results_T200000$df_output$T = 200000

df_output = rbind(results_T1000$df_output, 
                  results_T2000$df_output,
                  results_T5000$df_output,
                  results_T10000$df_output,
                  results_T20000$df_output,
                  results_T50000$df_output,
                  results_T100000$df_output,
                  results_T200000$df_output)

write.csv(data.frame(lapply(df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_truemarglikgaussmulti.csv')

### END SIMULATION STUDY G=5, R=6 ###

### BEGIN banknotes dataset ###
data("banknote")
y = banknote[, -1]
iters=10000
burn_in = 2000
sampler = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                     mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                     init,2*iters,burn_in=burn_in,seed=seed))
results_T10000_G2 = thames_pipeline(num_sims=1,
                                    logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                    p0hat=p0hat_gaussmulti,
                                    loglik_partial=NULL,
                                    G_list=2, iters=iters/2,
                                    relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                    thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                    num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                    prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                    samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                                          mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                          init,2*iters,burn_in=burn_in,seed=seed,y=y)))

# results_T10000_G2$df_output
# save(results_T10000_G2,file=paste0('data/res_banknotes_G2_T10000','.Rda'))
# load(file=paste0('data/res_banknotes_G2_T10000','.Rda'))

results_T10000_G3 = thames_pipeline(num_sims=1,
                                    logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                    p0hat=p0hat_gaussmulti,
                                    loglik_partial=NULL,
                                    G_list=3, iters=iters/2,
                                    relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                    thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                    num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                    prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                    samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                                          mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                          init,2*iters,burn_in=burn_in,seed=seed,y=y)))

# results_T10000_G3$df_output
# save(results_T10000_G3,file=paste0('data/res_banknotes_G3_T10000','.Rda'))
# load(file=paste0('data/res_banknotes_G3_T10000','.Rda'))

# run to visualize CO
# results_T10000_G4 = thames_pipeline(num_sims=1,
#                                     logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                     loglik_partial=NULL,
#                                     G_list=4, iters=iters/2,
#                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                     num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
#                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
#                                                                                                                           mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                           init,2*iters,burn_in=burn_in,seed=seed,y=y)))

results_T10000_G4 = thames_pipeline(num_sims=1,
                                    logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                    p0hat=p0hat_gaussmulti,
                                    loglik_partial=NULL,
                                    G_list=4, iters=iters/2,
                                    relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                    thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                    num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                    prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                    samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                                          mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                          init,2*iters,burn_in=burn_in,seed=seed,y=y)))

# save(results_T10000_G4,file=paste0('data/res_banknotes_G4_T10000','.Rda'))
# load(file=paste0('data/res_banknotes_G4_T10000','.Rda'))

# run to visualize CO
# results_T10000_G5 = thames_pipeline(num_sims=1,
#                                     logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
#                                     loglik_partial=NULL,
#                                     G_list=5, iters=iters/2,
#                                     relabel_algs = c("ECR"), ellipse_algs = c("standard"),
#                                     thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
#                                     num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
#                                     prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
#                                     samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
#                                                                                                                           mustars_func(g),sigmastars_func(g),taustars_func(g),
#                                                                                                                           init,2*iters,burn_in=burn_in,seed=seed,y=y)))

results_T10000_G5 = thames_pipeline(num_sims=1,
                                    logposty=function(thetas, G, y) logposty_gaussmulti_transformed(thetas,G,y,nu,rep(alpha_0_constant,G)),
                                    p0hat=p0hat_gaussmulti,
                                    loglik_partial=NULL,
                                    G_list=5, iters=iters/2,
                                    relabel_algs = c("ECR"), ellipse_algs = c("standard"),
                                    thames_algs = c("simple"),#c(permutations="permutations",simple="simple"),
                                    num_R=R,num_var_g=1+R+R*(R-1)/2+R,init=function(y) 1, # init is set within the function
                                    prior_sampler = function(y, G, iters) prior_sampler_marglik_gaussmulti(y, nu, iters, rep(alpha_0_constant,G)),
                                    samplers = list(function(g, iters, init, seed) y_theta_sampler_gaussmulti_transformed(n,nu,rep(alpha_0_constant,g),
                                                                                                                          mustars_func(g),sigmastars_func(g),taustars_func(g),
                                                                                                                          init,2*iters,burn_in=burn_in,seed=seed,y=y)))
# save(results_T10000_G5,file=paste0('data/res_banknotes_G5_T10000','.Rda'))
# load(file=paste0('data/res_banknotes_G5_T10000','.Rda'))

marglik_res = rep(0,4)
marglik_res[1] = results_T10000_G2$df_output$value[results_T10000_G2$df_output$estimate=="estim"]
marglik_res[2] = results_T10000_G3$df_output$value[results_T10000_G3$df_output$estimate=="estim"]
marglik_res[3] = compute_nobile_identity(logZhatGminus1 = marglik_res[2],
                                         p0hat_value = results_T10000_G4$df_output$value[1],
                                         dirichlet_vec = rep(1,4),
                                         n=nrow(y))
marglik_res[4] = compute_nobile_identity(logZhatGminus1 = marglik_res[3],
                                         p0hat_value = results_T10000_G5$df_output$value[1],
                                         dirichlet_vec = rep(1,5),
                                         n=nrow(y))
marglik_res = matrix(marglik_res,nrow=1)
colnames(marglik_res) = c("2","3","4","5")
write.csv(marglik_res,"data/marglik_res_banknote.csv")
read.csv("data/marglik_res_banknote.csv")

# co values can be read of from the graph visualization
co = rep(0,4)
co[1] = 2
co[2] = 3
co[3] = 3 - 1 
co[4] = 3 - 2
co = matrix(oic,nrow=1)
write.csv(oic,"data/oic_banknote.csv")
read.csv("data/oic_banknote.csv")

confusion_mat_banknote = rbind(c(sum((results_T10000_G3$I_map[[1]][1,1,]==1)&(banknote[,1]=="genuine")),
                                 sum((results_T10000_G3$I_map[[1]][1,1,]==1)&(banknote[,1]=="counterfeit"))),
                               c(sum((results_T10000_G3$I_map[[1]][1,1,]==2)&(banknote[,1]=="genuine")),
                                 sum((results_T10000_G3$I_map[[1]][1,1,]==2)&(banknote[,1]=="counterfeit"))),
                               c(sum((results_T10000_G3$I_map[[1]][1,1,]==3)&(banknote[,1]=="genuine")),
                                 sum((results_T10000_G3$I_map[[1]][1,1,]==3)&(banknote[,1]=="counterfeit"))))
colnames(confusion_mat_banknote) = c("genuine","counterfeit")
rownames(confusion_mat_banknote) = c("1","2","3")

write.csv(confusion_mat_banknote,"data/confusion_mat_banknote.csv")
read.csv("data/confusion_mat_banknote.csv")
### END banknotes dataset ###