# Fit the THAMES to the galaxies, enzyme, and acidity datasets
rm(list=ls())
source('functions/galaxies_funcs_squares.R')
source('functions/pipeline.R')
pacman::p_load(rstan,label.switching,combinat)
library(thamesmix)
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
logposty = function(sims,y) lp_gmm_marginal_transform(y,sims,m,R)
loglik_partial = function(sims,G) loglik_gmm_partial_transform(y,sims,G)

# pipeline settings
G_list <- 2:6 # 2 more than in Celeux et al. (2018)
iters = 5000 # T/2
num_var_g = 3 # number of variables per component
RELABEL_ALG_NAMES = c("ECR", "STEPHENS","PRA", 
                      "AIC","SJW") # names from label.switching 
relabel_algs = c("ECR","STEPHENS")#relabel(relabel_algs[j],sims,n,y,R,m,G,2*iters,z_dummy = z_dummy)
# options: c("PRA", "AIC","STEPHENS","SJW")

ellipse_algs = list(standard="standard")

#ellipse_algs = list(max_mode_mclust="max_mode_mclust", standard="standard",
#                    robust="robust",min_vol="min_vol_60%HPD")  

# options: robust="robust")#,min_vol="min_vol_60%HPD")
thames_algs = list(permutations="permutations",simple="simple")#,permutations="permutations")
# options: standard_pluslogGfac="standard_pluslogGfac",simple="simple")
samplers = list(jags=function(g, iters, init, seed) sim_jags(n=n,x=y,R2=R^2,m=m,k1=g,
                                                       iters=iters, init=init, seed=seed),
                stan=function(g, iters, init, seed) sim_stan(n,y,R,m,g,iters,init,seed))
num_sims=20

results_galaxies = thames_pipeline(num_sims, logposty, loglik_partial, G_list, iters, 
                                   relabel_algs, ellipse_algs, thames_algs,
                                   samplers,num_R=1,num_var_g = num_var_g) 
# save(results_galaxies,file=paste0('data/res_galaxies','.Rda'))
write.csv(data.frame(lapply(results_galaxies$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_galaxies.csv')

# average number of times the exact evaluation and the approximation differed 
sapply(2:6,function(g) mean((results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="permutations")&(results_galaxies$df_output$estimate=="estim"),]$value-results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="simple")&(results_galaxies$df_output$estimate=="estim"),]$value)[!is.na(results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="permutations")&(results_galaxies$df_output$estimate=="estim"),]$value-results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="simple")&(results_galaxies$df_output$estimate=="estim"),]$value)]==0))
# average size of the error
sapply(2:6,function(g) mean(abs(results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="permutations")&(results_galaxies$df_output$estimate=="estim"),]$value-results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="simple")&(results_galaxies$df_output$estimate=="estim"),]$value)[!is.na(results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="permutations")&(results_galaxies$df_output$estimate=="estim"),]$value-results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="simple")&(results_galaxies$df_output$estimate=="estim"),]$value)]))
# maximum size of the error
sapply(2:6,function(g) max(abs(results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="permutations")&(results_galaxies$df_output$estimate=="estim"),]$value-results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="simple")&(results_galaxies$df_output$estimate=="estim"),]$value)[!is.na(results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="permutations")&(results_galaxies$df_output$estimate=="estim"),]$value-results_galaxies$df_output[(results_galaxies$df_output$G==g)&(results_galaxies$df_output$thamesalg=="simple")&(results_galaxies$df_output$estimate=="estim"),]$value)]))

# for visualizations

for(g in seq_along(G_list)){
  df_galaxies = as.data.frame(results_galaxies$Ws[[g]])
  names(df_galaxies) = c(paste0(rep("mu",G_list[g]),1:G_list[g]),
                         paste0(rep("sigmasqu",G_list[g]),1:G_list[g]),
                         paste0(rep("tau",G_list[g]-1),1:(G_list[g]-1)),
                         paste0(rep("W",G_list[g]),1:G_list[g]),
                         paste0(rep("muGrid",G_list[g]),1:G_list[g]),
                         paste0(rep("sigmasquGrid",G_list[g]),1:G_list[g]),
                         paste0(rep("WGrid",G_list[g]),1:G_list[g]))
  write.csv(df_galaxies,file=paste0("data/thetaW_galaxies_G",G_list[g],".csv"))
}

# for(g in seq_along(G_list)){
#   df_galaxies = as.data.frame(cbind(results_galaxies$sims[[g]][,,1],
#                                     exp(results_galaxies$sims[[g]][,,2]),
#                                     results_galaxies$sims[[g]][,,3],
#                                     results_galaxies$Ws[[60*(g-1)+1]]))
#   names(df_galaxies) = c(paste0(rep("mu",G_list[g]),1:G_list[g]),
#                          paste0(rep("sigmasqu",G_list[g]),1:G_list[g]),
#                          paste0(rep("tau",G_list[g]),1:G_list[g]),
#                          paste0(rep("W",G_list[g]),1:G_list[g]))
#   write.csv(df_galaxies,file=paste0("data/thetaW_galaxies_G",G_list[g],".csv"))
# }

results_galaxies_1sim = thames_pipeline(1, logposty, loglik_partial, 2:15, 50000,
                                   c("ECR"), ellipse_algs, list(simple="simple"),
                                   list(jags=samplers[[1]]),num_R=1,num_var_g = 3) 

#save(results_galaxies_1sim,file=paste0('data/res_galaxies_1sim','.Rda'))
#load(file=paste0('data/res_galaxies_1sim','.Rda'))

write.csv(data.frame(lapply(results_galaxies_1sim$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_galaxies_1sim.csv')

df_1sim = read.csv('data/res_galaxies_1sim.csv')
df_1sim[df_1sim$estimate=="estim",]$value

graphs = results_galaxies_1sim$graphs
save(graphs,file=paste0('data/res_galaxies_1sim_graphs','.Rda'))
load(paste0('data/res_galaxies_1sim_graphs','.Rda'))
co = results_galaxies_1sim$co
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
logposty = function(sims,y) lp_gmm_marginal_transform(y,sims,m,R)
loglik_partial = function(sims,G) loglik_gmm_partial_transform(y,sims,G)

# old options
#radii = function(d) sqrt(d+1)*exp(0/d) # radius from Metodiev et al (2024)
#limits = list(function(lps) quantile(lps,.01))

results_enzyme = thames_pipeline(num_sims, logposty, loglik_partial, G_list, iters, 
                                 relabel_algs, ellipse_algs, thames_algs,
                                 samplers,num_R=1,num_var_g = 3) 
#save(results_enzyme,file=paste0('data/res_enzyme','.Rda'))
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
logposty = function(sims,y) lp_gmm_marginal_transform(y,sims,m,R)
loglik_partial = function(sims,G) loglik_gmm_partial_transform(y,sims,G)

source('functions/galaxies_funcs_squares.R')
results_acidity = thames_pipeline(num_sims, logposty, loglik_partial, G_list, iters, 
                                   relabel_algs, ellipse_algs, thames_algs,
                                   samplers,num_R=1,num_var_g = 3) # TODO test init parameter
#save(results_acidity,file=paste0('data/res_acidity','.Rda'))
write.csv(data.frame(lapply(results_acidity$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_acidity.csv')
### ACIDITY ###