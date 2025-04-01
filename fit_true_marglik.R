# simulate from the Gaussian mixture model with known proportions and variances
# and compute the marginal likelihood estimators
rm(list=ls())
if(strsplit(getwd(),"/")[[1]][length(strsplit(getwd(),"/")[[1]])]!="thames_mixtures"){
  setwd("thames_mixtures")
}
source("functions/galaxies_funcs_squares.R")
source("functions/true_marglik_funcs.R")          
#source("functions/thames_gmm_funcs.R")
source("functions/pipeline.R")
library(reshape2)
library(combinat)

### PRETTIER SOLUTION ###
n = 10
rho = c(0,0.5,1)
iters = 10000
offset = 4

# setting 1: true model (fitted G=2, true G=2)
results_G2 = thames_pipeline(num_sims=50, 
                                  logposty= function(thetas, G, y) logposty(thetas, G, y, sigma_tilde=1, taus_tilde=c(1/2,1/2), mus_tilde=c(0,0)), 
                                  loglik_partial=NULL, 
                                  G_list=2, iters=iters/2, 
                                   relabel_algs=c("ECR"), ellipse_algs=c("standard"),
                                  thames_algs=c(simple="simple",mc="mc",bridge="bridge",standard="standard"),
                                  num_R=1,num_var_g=1, init=function(y) mean(range(y))*c(1,1),
                                  prior_sampler=function(y, G, iters) prior_sampler_marglik(y=y,mus_tilde=c(0,0), sigma_tilde=1, taus_tilde=c(1/2,1/2), iters=iters),
                                  samplers=list(function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0), 
                                                                                            sigma_tilde=1, 
                                                                                            taus_tilde=c(1/2,1/2),
                                                                                            mustars=offset+(2*(1-rho[1])+1)*c(-1,1),
                                                                                            sigmastar=1,
                                                                                            taustars=c(1/3,2/3),
                                                                                            init,2*iters, name="truemodel_G2_rho0",seed=seed),
                                                function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0), 
                                                                                            sigma_tilde=1, 
                                                                                            taus_tilde=c(1/2,1/2),
                                                                                            mustars=offset+(2*(1-rho[2])+1)*c(-1,1),
                                                                                            sigmastar=1,
                                                                                            taustars=c(1/3,2/3),
                                                                                            init,2*iters, name="truemodel_G2_rho.5",seed=seed),
                                                function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0), 
                                                                                            sigma_tilde=1, 
                                                                                            taus_tilde=c(1/2,1/2),
                                                                                            mustars=offset+(2*(1-rho[3])+1)*c(-1,1),
                                                                                            sigmastar=1,
                                                                                            taustars=c(1/3,2/3),
                                                                                            init,2*iters, name="truemodel_G2_rho1",seed=seed)))


# setting 2: underfitting (fitted G=2, true G=3)
results_underfitting = thames_pipeline(num_sims=50, 
                                  logposty= function(thetas, G, y) logposty(thetas, G, y, sigma_tilde=1, taus_tilde=c(1/2,1/2), mus_tilde=c(0,0)), 
                                  loglik_partial=NULL, 
                                  G_list=2, iters=iters/2, 
                                  relabel_algs=c("ECR"), ellipse_algs=c("standard"),
                                  thames_algs=c(mc="mc",bridge="bridge",simple="simple",standard="standard"),
                                  num_R=1,num_var_g=1,init=function(y) mean(range(y))*c(1,1),
                                  prior_sampler=function(y, G, iters) prior_sampler_marglik(y=y,mus_tilde=c(0,0), sigma_tilde=1, taus_tilde=c(1/2,1/2), iters=iters),
                                  samplers=list(function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0), 
                                                                                            sigma_tilde=1, 
                                                                                            taus_tilde=c(1/2,1/2),
                                                                                            mustars=offset+(2*(1-rho[1])+1)*c(-1,0,1),
                                                                                            sigmastar=1,
                                                                                            taustars=c(2/6,1/6,3/6),
                                                                                            init,2*iters, name="underfitting_rho0",seed=seed),
                                                function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0), 
                                                                                            sigma_tilde=1, 
                                                                                            taus_tilde=c(1/2,1/2),
                                                                                            mustars=offset+(2*(1-rho[2])+1)*c(-1,0,1),
                                                                                            sigmastar=1,
                                                                                            taustars=c(2/6,1/6,3/6),
                                                                                            init,2*iters, name="underfitting_rho.5",seed=seed),
                                                function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0), 
                                                                                            sigma_tilde=1, 
                                                                                            taus_tilde=c(1/2,1/2),
                                                                                            mustars=offset+(2*(1-rho[3])+1)*c(-1,0,1),
                                                                                            sigmastar=1,
                                                                                            taustars=c(2/6,1/6,3/6),
                                                                                            init,2*iters, name="underfitting_rho1",seed=seed)))


# setting 3: overfitting (fitted G=3, true G=2)
results_overfitting = thames_pipeline(num_sims=50, 
                                       logposty= function(thetas, G, y) logposty(thetas, G, y, sigma_tilde=1, taus_tilde=c(1/3,1/3,1/3), mus_tilde=c(0,0,0)), 
                                       loglik_partial=NULL, 
                                       G_list=3, iters=iters/2, 
                                       relabel_algs=c("ECR"), ellipse_algs=c("standard"),
                                       thames_algs=c(mc="mc",bridge="bridge",simple="simple",standard="standard"),
                                       num_R=1,num_var_g=1,init=function(y) mean(range(y))*c(1,1,1),
                                       prior_sampler=function(y, G, iters) prior_sampler_marglik(y=y,mus_tilde=c(0,0,0), sigma_tilde=1, taus_tilde=c(1/3,1/3,1/3), iters=iters),
                                       samplers=list(function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0,0), 
                                                                                                 sigma_tilde=1, 
                                                                                                 taus_tilde=c(1/3,1/3,1/3),
                                                                                                 mustars=offset+(2*(1-rho[1])+1)*c(-1,1),
                                                                                                 sigmastar=1,
                                                                                                 taustars=c(1/3,2/3),
                                                                                                 init,2*iters, name="overfitting_rho0",seed=seed),
                                                     function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0,0), 
                                                                                                 sigma_tilde=1, 
                                                                                                 taus_tilde=c(1/3,1/3,1/3),
                                                                                                 mustars=offset+(2*(1-rho[2])+1)*c(-1,1),
                                                                                                 sigmastar=1,
                                                                                                 taustars=c(1/3,2/3),
                                                                                                 init,2*iters, name="overfitting_rho.5",seed=seed),
                                                     function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0,0), 
                                                                                                 sigma_tilde=1, 
                                                                                                 taus_tilde=c(1/3,1/3,1/3),
                                                                                                 mustars=offset+(2*(1-rho[3])+1)*c(-1,1),
                                                                                                 sigmastar=1,
                                                                                                 taustars=c(1/3,2/3),
                                                                                                 init,2*iters, name="overfitting_rho1",seed=seed)))

# setting 4: true model (fitted G=3, true G=3)
results_G3 = thames_pipeline(num_sims=50, 
                                      logposty= function(thetas, G, y) logposty(thetas, G, y, sigma_tilde=1, taus_tilde=c(1/3,1/3,1/3), mus_tilde=c(0,0,0)), 
                                      loglik_partial=NULL, 
                                      G_list=3, iters=iters/2, 
                                      relabel_algs=c("ECR"), ellipse_algs=c("standard"),
                                      thames_algs=c(mc="mc",bridge="bridge",simple="simple",standard="standard"),
                                      num_R=1,num_var_g=1,init=function(y) mean(range(y))*c(1,1,1),
                                      prior_sampler=function(y, G, iters) prior_sampler_marglik(y=y,mus_tilde=c(0,0,0), sigma_tilde=1, taus_tilde=c(1/3,1/3,1/3), iters=iters),
                                      samplers=list(function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0,0), 
                                                                                                sigma_tilde=1, 
                                                                                                taus_tilde=c(1/3,1/3,1/3),
                                                                                                mustars=offset+(2*(1-rho[1])+1)*c(-1,0,1),
                                                                                                sigmastar=1,
                                                                                                taustars=c(2/6,1/6,3/6),
                                                                                                init,2*iters, name="truemodel_G3_rho0",seed=seed),
                                                    function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0,0), 
                                                                                                sigma_tilde=1, 
                                                                                                taus_tilde=c(1/3,1/3,1/3),
                                                                                                mustars=offset+(2*(1-rho[2])+1)*c(-1,0,1),
                                                                                                sigmastar=1,
                                                                                                taustars=c(2/6,1/6,3/6),
                                                                                                init,2*iters, name="truemodel_G3_rho.5",seed=seed),
                                                    function(g, iters, init, seed) y_mu_sampler(n=n,mus_tilde=c(0,0,0), 
                                                                                                sigma_tilde=1, 
                                                                                                taus_tilde=c(1/3,1/3,1/3),
                                                                                                mustars=offset+(2*(1-rho[3])+1)*c(-1,0,1),
                                                                                                sigmastar=1,
                                                                                                taustars=c(2/6,1/6,3/6),
                                                                                                init,2*iters, name="truemodel_G3_rho1",seed=seed)))

save(results_G2,file=paste0('data/res_sim1_G2','.Rda'))
save(results_underfitting,file=paste0('data/res_sim1_underfitting','.Rda'))
save(results_overfitting,file=paste0('data/res_sim1_overfitting','.Rda'))
save(results_G3,file=paste0('data/res_sim1_G3','.Rda'))

write.csv(data.frame(lapply(results_G2$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_sim1_G2.csv')
write.csv(data.frame(lapply(results_underfitting$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_sim1_underfitting.csv')
write.csv(data.frame(lapply(results_overfitting$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_sim1_overfitting.csv')
write.csv(data.frame(lapply(results_G3$df_output,
                            as.character), stringsAsFactors=FALSE),
          file='data/res_sim1_G3.csv')