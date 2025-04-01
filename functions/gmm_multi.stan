data {
 int p; //dimension
 int G; //number of mixture components
 int n; //sample size
 vector[p] x[n]; //data
 
 real<lower=0> alpha0; //dirichlet prior hyperparameter
 real<lower=0> nu0; //Wishart prior hyperparameter
 real<lower=0> beta0; // scale hyperparameter
 vector[p] m0; //empirical mean
 cov_matrix[p] W0; //empirical covariance
}

parameters {
 simplex[G] theta; //mixing proportions
 vector[p] mu[G]; //mixture component means
 // cholesky_factor_corr[p] L[G]; //cholesky factor of covariance
 cov_matrix[p] Sigma[G];
}

model {
  real ps[G];
  
  // prior
  theta ~ dirichlet(rep_vector(alpha0,G)); 
  for(g in 1:G){
    mu[g] ~ multi_normal(m0,Sigma[g]*beta0);
    Sigma[g] ~ inv_wishart(nu0, W0);
  }

  // log likelihood
  for (j in 1:n){
    for (g in 1:G){
      ps[g] = log(theta[g])+multi_normal_lpdf(x[j] | mu[g], Sigma[g]); 
      // ps[g] = log(theta[k])+multi_normal_cholesky_lpdf(y[n] | mu[k], L[k]); 
      //increment log probability
    }
    target += log_sum_exp(ps);
  }
}
