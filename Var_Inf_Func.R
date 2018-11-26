library(numDeriv)
library(quantreg)
library(bayesQR)
#Data simulation
set.seed(66)
n <- 1000
x <- runif(n=n,min=0,max=1)
y <- 1 + 2*x + 3*x^2 + rnorm(n=n, mean=0, sd=.6) #beta_p=c(1,2,3)
X <- cbind(1, x, x^2)
colnames(X) <- NULL

#This is the case using prior z ~ Exp(1) and beta_p ~ N(mu_p0, Sigma_p0)
#Approximation of expectation of GIG can be used but this might lead to negative determinant of var-cov matrix
variational_inf_ALD <- function(y, X, p, sigma=1, tol=10^(-5)){
  #Prior dist.
  #The dist. z_i is exponential 1 by default.
  mu_p0 <- rep(mean(y), ncol(X)) #Prior for beta_p
  Sigma_p0 <- diag(rep(var(y), ncol(X)))
  
  #Initial value for q*(beta_p) 
  
  init.mu <- rep(0, ncol(X)) 
  init.Sigma <- diag(ncol(X))
  
  #Paramter in ALD:
  #y_i = x_i*beta + theta*Exp(z_i) +tau*sqrt(sigma*z_i)*N(0,1)
  
  theta <- (1-2*p)/(p*(1-p))
  tau2 <- 2/(p*(1-p))
  
  mu_q <- init.mu
  Sigma_q <- init.Sigma
  #Absolute change of lower bound between two iterations
  delta <- 1 
  # Number of iteration 
  iter <- 0 
  #Previous lower bound
  prev_l <- 0 
  
  while(delta >= tol){
    
    iter <- iter + 1
    #Update parameter in dist. of z_i, q*(z_i) ~ GIG(1/2, a_i, b_i)
    a <- theta^2/(tau2*sigma) + 2*sigma
    #Expectation of t(beta_p)%*%beta_p Square matrix, ncol = ncol(X)
    beta2 <- mu_q%*%t(mu_q) + Sigma_q 
    b <- (y^2 - 2*y*(X%*%mu_q) + diag(X%*%beta2%*%t(X)))/(tau2*sigma)
    
    #Expectation of 1/z_i and log(z_i)
    bessel <- besselK(sqrt(a*b),nu=3/2,expon.scaled=FALSE)/
      besselK(sqrt(a*b),nu=1/2,expon.scaled=FALSE)
    
    e_z <- sqrt(b/a)*bessel
    
    recp_z <- sqrt(a/b)*bessel - 1/b
    
    partial <- grad(func=besselK,sqrt(b*a),nu=1/2,
                    expon.scaled=FALSE,method.args = list(eps = 1e-8, show.details = FALSE))
    
    ln_z <- log(sqrt(b/a)) + 
      partial/besselK(sqrt(b*a),nu=1/2,expon.scaled=FALSE)
    
    e_z <- as.vector(e_z)
    recp_z <- as.vector(recp_z)
    ln_z <- as.vector(ln_z)
    
    #Update parameter in dist. of beta_p, q*(beta_p) ~ N(mu_q, Sigma_q)
    Sigma_q <- solve(t(X)%*%(X*recp_z)/(tau2*sigma) + solve(Sigma_p0))
    
    M <- t(X*y*recp_z)/(tau2*sigma)
    
    mu_q <- Sigma_q%*%(apply(M, 1, sum)- 
                         apply(t(X), 1, sum)*theta/(tau2*sigma)+ solve(Sigma_p0)%*%mu_p0)
    
    #Calculate the lower bound for each iteration
    E_log_y <- sum(ln_z)
    
    E_logq_beta <- -1/2*log(det(Sigma_q))
    
    E_logp_beta <- t(mu_q-mu_p0)%*%solve(Sigma_p0)%*%(mu_q-mu_p0)+ 
      sum(diag(Sigma_q%*%solve(Sigma_p0)))
    
    E_logq_z <- -1/2*sum(ln_z) - 1/2*sum(a*recp_z+b*e_z)
    
    E_logp_z <- sum(-sigma*e_z)
    
    l <- E_log_y + E_logq_beta + E_logp_beta +E_logq_z +E_logp_z
    
    delta <- abs(l - prev_l)
    prev_l <- l
  }
  
  quan_var_inf <- list(mu_q,Sigma_q,iter,p)
  names(quan_var_inf) <- c("mu_q","Sigma_q","niter","quantile")
  return(quan_var_inf)
}

Var_inf <- function(y,X,p,sigma=1, LASSO=FALSE,tol=10^(-5)){
  beta_p <- matrix(NA, ncol=length(p),nrow=ncol(X))
  for(i in 1:length(p)){
    beta_p[,i] <- apply(variational_inf_ALD(y,X,p[i], sigma=sigma, tol=tol)$mu_q, 1, mean)
  }
  colnames(beta_p) <- as.character(p[1:length(p)])
  return(beta_p)
}

#Computing Mean Squre Error of the variational inference
MSE <- function(y,X, p, beta, beta_p, mu, sd){
  mse <- as.numeric()
  for(i in 1:length(p)){
    quan <- p[i]
    mse[i] <- sum((X%*%beta + qnorm(quan,mean=mu,sd=sd) - X%*%beta_p[,i])^2)}
  return(mse)
}

#Comparasion
#The results are pretty close especially when p is around 0.5. When p is moves closer to 0 or 1, the approximation
#might be slightly different

model <- bayesQR(y~X[,2:3],quantile=c(0.05,0.25,0.5,0.75,0.95),alasso = FALSE, ndraw = 500)
summary(model)

model1_woLASSO <- Var_inf(y,X,c(0.05,0.25,0.50,0.75,0.95))
model1_woLASSO
MSE(y,X,c(0.05,0.25,0.50,0.75,0.95),c(1,2,3),model1_woLASSO,mu=0,sd=0.6)
system.time( Var_inf(y,X,c(0.05,0.25,0.50,0.75,0.95)))

#The result given by package bayesQR where gibbs sampling is used instead
model2 <- bayesQR(y~X[,2:3],quantile=c(0.05,0.25,0.5,0.75,0.95),alasso = FALSE, ndraw = 500)
summary(model2)
#The result given by package quantreg
model3 <- rq(y~X[,2:3],tau=c(0.05,0.25,0.5,0.75,0.95))
MSE(y,X,c(0.05,0.25,0.50,0.75,0.95),c(1,2,3),model3$coefficients,mu=0,sd=0.6)