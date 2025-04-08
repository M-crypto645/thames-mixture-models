# Fit the THAMES to the galaxies, enzyme, and acidity datasets
rm(list=ls())
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
# save(results_galaxies,file=paste0('data/res_galaxies','.Rda'))
write.csv(data.frame(lapply(results_galaxies$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_galaxies.csv')

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
co = sapply(seq_along(graphs),function(s) sum(V(graphs[[s]])$color == "blue") - sum(V(graphs[[s]])$color == "red"))
co
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
logposty = function(theta,G,y) lp_gmm_marginal_transform(y,theta,G,m,R)
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