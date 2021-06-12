#setwd("~/Research/logit_occ/OfficeCluster")
#-------------------------------------------------------------------------------

#Windows version
#-------------------------------------------------------------------------------
require(parallel)
detectCores()

require(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)
#-------------------------------------------------------------------------------

#Packages and functions that should be loaded
#-------------------------------------------------------------------------------
require(Rcpp)
require(RcppArmadillo)
require(jagsUI)
require(RcppOccupancy2)
require(coda)

require(rstan)
rstan_options(auto_write = TRUE)
#-------------------------------------------------------------------------------

#Bayesian jags code
#-------------------------------------------------------------------------------
writeLines("
           model{
           #likelihood part
           for (i in 1:n) #loop over each site
           {
           #loop over the number of visits per site
           #make general later
           #ie the number of surveys should be generalized
           #beta vector and alpha vector should be generalized

           #state occupancy
           zb[i] ~ dbern(pz[i])
           logit(pz[i]) <- beta0 + beta1*X[i]

           for (j in 1:J) #loop over the number of visits per site #make general later
           {
           #observation process
           Y[i,j] ~ dbern(py[i,j])
           py[i,j] <- zb[i]*pd[i,j]
           logit(pd[i,j])<- alpha0 + alpha1*W[i,j,2]  #The W should be generalized later!!!
           }#end loop over surveys
           }#end loop over sites

           alpha0 ~ dnorm(0, 0.001)
           alpha1 ~ dnorm(0, 0.001)

           beta0 ~ dnorm(0, 0.001)
           beta1 ~ dnorm(0, 0.001)

           #occupied <-sum(zb[])
}##model
", con = "occ.txt")
#-------------------------------------------------------------------------------

scode <- "
data {
  int<lower=1> N; //indicates that the variable is an integer with smallest value = 0
  int<lower=1> V; //number of visits
  int<lower=0,upper=1> y[N,V]; //presence absence data

  vector[N] Xmat; // X variable
  real Wmat[N, V]; // W matrix - detection covariates
}

parameters {
  vector[2] beta;  //occupancy slope params
  vector[2] alpha;  //detection slope params
}

transformed parameters {
  vector[N] psi; // prob of of occurrence
  real pij[N, V]; // prob of detection

  //calculate psi_i and pij
  for (isite in 1:N){
    psi[isite] = inv_logit( beta[1] + Xmat[isite]*beta[2] );

    for (ivisit in 1:V){
      pij[isite, ivisit] = inv_logit( alpha[1] + Wmat[isite, ivisit]*alpha[2] );
    }
  }
}

model {
  // local variables to avoid recomputing log(psi) and log(1-psi)
  vector[N] log_psi;
  vector[N] log1m_psi;

  for (isite in 1:N) {
    log_psi[isite] = log(psi[isite]);
    log1m_psi[isite] = log1m(psi[isite]);
  }

  // priors
  alpha ~ normal(0, 1000);
  beta ~ normal(0, 1000);

  // likelihood
  for (isite in 1:N) {

    if (sum(y[isite]) > 0){
      target += log_psi[isite] + bernoulli_lpmf(y[isite]|pij[isite]) ;
    }else {
      target += log_sum_exp(log_psi[isite] + bernoulli_lpmf(y[isite]|pij[isite]),log1m_psi[isite]);
    }
  }//end likelihood contribution
}

"

Model_code_covocc <- stan_model(model_code = scode)

#-------------------------------------------------------------------------------

#MLE function for the single season occupancy model

negLL = function(param) 
{
  beta = param[1:dim(X)[2]]
  #psi = pnorm(X %*% beta)
  psi = as.vector(1/(1+exp(-X %*% beta))) #logistic link function used
  alpha = param[(dim(X)[2]+1):(dim(X)[2]+dim(W)[3])]
  p = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    #p[, j] = pnorm(W[,j,] %*% alpha)
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha)))  
  }
  logL = rep(NA,n)
  for (i in 1:n) {
    yvec = y[i, jind[i,]]
    pvec = p[i, jind[i,]]
    terms = pvec^yvec * (1-pvec)^(1-yvec)
    logL[i] = log( psi[i] * prod(terms)  + ifelse(ysum[i]==0,1,0)*(1-psi[i])  )
  }
  (-1)*sum(logL)
}
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.5", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(0, 1.75)

n=50
J=3
iterations <- 500
up=100000

output1<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
{
     #Simulate data

     set.seed(iter+up)

     #initial generation of the sample data
     x1 = runif(n, -2,2)
     x = (x1 - mean(x1)) / sd(x1)
     X  = cbind(rep(1,n), x)
     psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
     z = rbinom(n, size=1, prob=psi)

     w1 = runif(n*J, -5,5)
     w = (w1 - mean(w1)) / sd(w1)
     w = matrix(w, nrow=n, ncol=J)
     W = array(dim=c(n,J,2))
     W[,,1] = 1
     W[,,2] = w

     p = matrix(nrow=n, ncol=J)
     y = matrix(nrow=n, ncol=J)
     for (j in 1:J) {
          p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
          y[, j] = rbinom(n, size=1, prob=z*p[, j])
     }

     jind = !is.na(y)   # index for non-missing observations of y
     ysum = apply(y,1,sum, na.rm=TRUE)

     Y.eg<-y
     X.eg = as.data.frame(x) #siteCovs
     colnames(X.eg)<-c("X1")
     W1=matrix(NA,nrow=n, ncol=J)
     W1[1:n,]<-W[1:n,,2]
     W.eg.l1<-list(W1=W1)

     #-------------------------------------------------------------------------------
     
     #mle fit
     
     jind = !is.na(y)  # index for non-missing observations of y
     ysum = apply(y,1,sum, na.rm=TRUE)
     betaGuess = rep(0, dim(X)[2])
     alphaGuess = rep(0, dim(W)[3])
     paramGuess = c(betaGuess, alphaGuess)
     fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
     
     mle = fit$par
     
     #-------------------------------------------------------------------------------
     
     if (fit$convergence==0){
       
       #covariance matrix based on the mle fit
       
       vcv = chol2inv(chol(fit$hessian))

       #do the bayesian fit
       #**************************
       # Input data for WinBUGS
       #**************************
       occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
  
       #*********************
       # Initial values
       #*********************
  
       zst<-apply(y,1,max)
       
       
       inits <- function()
       {
         if ( mean((mle[1:2]-beta.param))>10 |  mean((mle[3:4]-alpha.param))>10 ){
           
           list(zb=zst,
                beta0=0,
                beta1=0,
                alpha0=0,
                alpha1=0)
           
         }else{
           
           list(zb=zst,
                beta0=mle[1],
                beta1=mle[2],
                alpha0=mle[3],
                alpha1=mle[4] )
         }
       }
  
       #****************************
       # Paramaters to be monitored
       #****************************
  
       parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
  
       #print("start mcmc")
       nsims<-20000
       niter <- nsims
       nburn <- floor(niter/2)
       nthin <- 1
       nchains <- 1
       nkeep<-niter-nburn
  
  
       t1<-proc.time()
       occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                        n.chains=1, n.iter=niter, n.burnin=nburn,
                        parallel = FALSE, verbose=FALSE)
       t2<-proc.time()
       jagsTimer<-t2[3]-t1[3]
  
       jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
       effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
       effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
       effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
       
       #alpha0 = occ.jags$sims.list$alpha0
       #alpha1 = occ.jags$sims.list$alpha1
       #beta0 = occ.jags$sims.list$beta0
       #beta1 = occ.jags$sims.list$beta1
  
       #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
       rm(occ.jags)
       #-------------------------------------------------------------------------------
  
       design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
  
       #set the formula to use
       form1<- V1 ~ X1 ~ W1
  
       #the prior distributions specified here
       alpha_m<-matrix( mle[3:4], ncol=1)
       beta_m<-matrix( mle[1:2], ncol=1)
       sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
       sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
  
  
       t1_drum5<-proc.time()
            occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
       t2_drum5<-proc.time()
       drumTimer<-t2_drum5[3]-t1_drum5[3]
       
       dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                  effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                  effectiveSize(as.mcmc(occdRUM$beta[1,])),
                  effectiveSize(as.mcmc(occdRUM$beta[2,])))
       
       rm(occdRUM)
       #-------------------------------------------------------------------------------
       
       #Fit using Stan
       N<-n
       V<-J
       Wmat<-W[,,2]
       Xmat<-X[,2]
  
       Stan1<-proc.time()
       occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                             iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
       Stan2<-proc.time()
       Stantimer<-Stan2[3]-Stan1[3]
  
       post_occregfit<-extract(occregfit, permuted = F)
       StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
  
       rm(post_occregfit)
       #-------------------------------------------------------------------------------
  
       c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
     }else{
       c(rep(NA, 15))
     }
}
cat("\n output1 DONE! \n")
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
  
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}

output<-output1
output1<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J3_case1_EURING_revision.RData")
rm(output1) #remove output 1
cat("\n output1 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.5", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(0, 1.75)

n=50
J=5
iterations <- 500
up=500

output2<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                  
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}

output<-output2
output2<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J5_case1_EURING_revision.RData")
rm(output2) 
cat("\n output2 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.5", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(0, 1.75)

n=100
J=3
iterations <- 500
up=5000

output3<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                    
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output3

output3<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J3_case1_EURING_revision.RData")
rm(output3) #remove output 3
cat("\n output3 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.5", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(0, 1.75)

n=100
J=5
iterations <- 500
up=50000

output4<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output4

output4<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J5_case1_EURING_revision.RData")
rm(output4) #remove output 4
cat("\n output4 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.5", "; mean(psi)=0.5", sep="")
alpha.param = c(0, 1.75)
beta.param = c(-0.1, 2.5)

n=50
J=3
iterations <- 500
up=1000

output5<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output5

output5<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J3_case2_EURING_revision.RData")
rm(output5) #remove output 5

cat("\n output5 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.5", "; mean(psi)=0.5", sep="")
alpha.param = c(0, 1.75)
beta.param = c(-0.1, 2.5)

n=50
J=5
iterations <- 500
up=2000

output6<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output6

output6<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J5_case2_EURING_revision.RData")
rm(output6) #remove output 6

cat("\n output6 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.5", "; mean(psi)=0.5", sep="")
alpha.param = c(0, 1.75)
beta.param = c(-0.1, 2.5)

n=100
J=3
iterations <- 500
up=3000

output7<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output7

output7<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J3_case2_EURING_revision.RData")
rm(output7) #remove output 7

cat("\n output7 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.5", "; mean(psi)=0.5", sep="")
alpha.param = c(0, 1.75)
beta.param = c(-0.1, 2.5)

n=100
J=5
iterations <- 500
up=4000

output8<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output8

output8<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J5_case2_EURING_revision.RData")
rm(output8) #remove output 8

cat("\n output8 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.7", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(1.35, 1.75)

n=50
J=3
iterations <- 500
up=0

output13<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      zst<-apply(y,1,max)
                      
                      
                      inits <- function()
                      {
                        if ( mean((mle[1:2]-beta.param))>10 |  mean((mle[3:4]-alpha.param))>10 ){
                          
                          list(zb=zst,
                               beta0=0,
                               beta1=0,
                               alpha0=0,
                               alpha1=0)
                          
                        }else{
                          
                          list(zb=zst,
                               beta0=mle[1],
                               beta1=mle[2],
                               alpha0=mle[3],
                               alpha1=mle[4] )
                        }
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output13

output13<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J3_case5_EURING_revision.RData")
rm(output13) #remove output 13

cat("\n output13 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.7", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(1.35, 1.75)

n=50
J=5
iterations <- 500
up=500

output14<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output14

output14<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J5_case5_EURING_revision.RData")
rm(output14) #remove output 14

cat("\n output14 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.7", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(1.35, 1.75)

n=100
J=3
iterations <- 500
up=5000

output15<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output15

output15<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J3_case5_EURING_revision.RData")
rm(output15) #remove output 15

cat("\n output15 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.7", "; mean(psi)=0.3", sep="")
alpha.param = c(1.35, 1.75)
beta.param = c(-1.85, 2.5)

n=100
J=5
iterations <- 500
up=50000

output16<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output16

output16<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J5_case5_EURING_revision.RData")
rm(output16) #remove output 16

cat("\n output16 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.7", "; mean(psi)=0.5", sep="")
alpha.param = c(1.35, 1.75)
beta.param = c(-0.1, 2.5)

n=50
J=3
iterations <- 500
up=1000

output17<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output17

output17<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J3_case6_EURING_revision.RData")
rm(output17) #remove output 17

cat("\n output17 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=50
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.7", "; mean(psi)=0.5", sep="")
alpha.param = c(1.35, 1.75)
beta.param = c(-0.1, 2.5)

n=50
J=5
iterations <- 500
up=2000

output18<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output18

output18<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n50_J5_case6_EURING_revision.RData")
rm(output18) #remove output 18

cat("\n output18 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=3, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.7", "; mean(psi)=0.5", sep="")
alpha.param = c(1.35, 1.75)
beta.param = c(-0.1, 2.5)

n=100
J=3
iterations <- 500
up=3000

output19<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output19

output19<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J3_case6_EURING_revision.RData")
rm(output19) 

cat("\n output19 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=100
#-------------------------------------------------------------------------------

case<-paste("alpha=c(1.35, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.7", "; mean(psi)=0.5", sep="")
alpha.param = c(1.35, 1.75)
beta.param = c(-0.1, 2.5)

n=100
J=5
iterations <- 500
up=4000

output20<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      inits <- function()
                      {
                        list("zb"=apply(y,1,max),
                             "beta0"=mle[1],
                             "beta1"=mle[2],
                             "alpha0"=mle[3],
                             "alpha1"=mle[4] )
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output20

output20<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n100_J5_case6_EURING_revision.RData")
rm(output20) #remove output 20

cat("\n output20 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#larger data sets
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#Simulations run below
#for j=5, n=500
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.5", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(0, 1.75)

n=500
J=5
iterations <- 500
up=11000

output9<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      zst<-apply(y,1,max)
                      
                      
                      inits <- function()
                      {
                        if ( mean((mle[1:2]-beta.param))>10 |  mean((mle[3:4]-alpha.param))>10 ){
                          
                          list(zb=zst,
                               beta0=0,
                               beta1=0,
                               alpha0=0,
                               alpha1=0)
                          
                        }else{
                          
                          list(zb=zst,
                               beta0=mle[1],
                               beta1=mle[2],
                               alpha0=mle[3],
                               alpha1=mle[4] )
                        }
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output9

output9<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n500_J5_case3_EURING_revision.RData")
rm(output9) #remove output 9

cat("\n output9 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------


#larger data sets
#-------------------------------------------------------------------------------
#Simulations run below
#for j=10, n=500
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-1.85, 2.5)","; mean(p)=0.5", "; mean(psi)=0.3", sep="")
beta.param = c(-1.85, 2.5)
alpha.param = c(0, 1.75)

n=500
J=10
iterations <- 500
up=12000

output10<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      zst<-apply(y,1,max)
                      
                      
                      inits <- function()
                      {
                        if ( mean((mle[1:2]-beta.param))>10 |  mean((mle[3:4]-alpha.param))>10 ){
                          
                          list(zb=zst,
                               beta0=0,
                               beta1=0,
                               alpha0=0,
                               alpha1=0)
                          
                        }else{
                          
                          list(zb=zst,
                               beta0=mle[1],
                               beta1=mle[2],
                               alpha0=mle[3],
                               alpha1=mle[4] )
                        }
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output10

output10<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n500_J10_case3_EURING_revision.RData")
rm(output10) #remove output 10

cat("\n output10 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#Simulations run below
#for j=5, n=500
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.5", "; mean(psi)=0.5", sep="")
alpha.param = c(0, 1.75)
beta.param = c(-0.1, 2.5)

n=500
J=5
iterations <- 500
up=13000

output11<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      zst<-apply(y,1,max)
                      
                      
                      inits <- function()
                      {
                        if ( mean((mle[1:2]-beta.param))>10 |  mean((mle[3:4]-alpha.param))>10 ){
                          
                          list(zb=zst,
                               beta0=0,
                               beta1=0,
                               alpha0=0,
                               alpha1=0)
                          
                        }else{
                          
                          list(zb=zst,
                               beta0=mle[1],
                               beta1=mle[2],
                               alpha0=mle[3],
                               alpha1=mle[4] )
                        }
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output11

output11<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n500_J5_case4_EURING_revision.RData")
rm(output11) #remove output 11

cat("\n output11 DONE! \n")

#open cluster again
cl <- makeCluster(16)
registerDoParallel(cl)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#Simulations run below
#for j=10, n=500
#-------------------------------------------------------------------------------

case<-paste("alpha=c(0, 1.75)", "; beta=c(-0.1, 2.5)","; mean(p)=0.5", "; mean(psi)=0.5", sep="")
alpha.param = c(0, 1.75)
beta.param = c(-0.1, 2.5)

n=500
J=10
iterations <- 500
up=14000

output12<- foreach(iter = 1:iterations,.combine='rbind',
                  .packages = c("jagsUI", "Rcpp", "RcppArmadillo", "RcppOccupancy2", "rstan", "coda"),
                  .inorder=T) %dopar%
                  {
                    #Simulate data
                    
                    set.seed(iter+up)
                    
                    #initial generation of the sample data
                    x1 = runif(n, -2,2)
                    x = (x1 - mean(x1)) / sd(x1)
                    X  = cbind(rep(1,n), x)
                    psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
                    z = rbinom(n, size=1, prob=psi)
                    
                    w1 = runif(n*J, -5,5)
                    w = (w1 - mean(w1)) / sd(w1)
                    w = matrix(w, nrow=n, ncol=J)
                    W = array(dim=c(n,J,2))
                    W[,,1] = 1
                    W[,,2] = w
                    
                    p = matrix(nrow=n, ncol=J)
                    y = matrix(nrow=n, ncol=J)
                    for (j in 1:J) {
                      p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
                      y[, j] = rbinom(n, size=1, prob=z*p[, j])
                    }
                    
                    jind = !is.na(y)   # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    
                    Y.eg<-y
                    X.eg = as.data.frame(x) #siteCovs
                    colnames(X.eg)<-c("X1")
                    W1=matrix(NA,nrow=n, ncol=J)
                    W1[1:n,]<-W[1:n,,2]
                    W.eg.l1<-list(W1=W1)
                    
                    #-------------------------------------------------------------------------------
                    
                    #mle fit
                    
                    jind = !is.na(y)  # index for non-missing observations of y
                    ysum = apply(y,1,sum, na.rm=TRUE)
                    betaGuess = rep(0, dim(X)[2])
                    alphaGuess = rep(0, dim(W)[3])
                    paramGuess = c(betaGuess, alphaGuess)
                    fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
                    
                    mle = fit$par
                    
                    #-------------------------------------------------------------------------------
                    
                    if (fit$convergence==0){
                      
                      #covariance matrix based on the mle fit
                      
                      vcv = chol2inv(chol(fit$hessian))
                      
                      #do the bayesian fit
                      #**************************
                      # Input data for WinBUGS
                      #**************************
                      occ.data <- list(Y=y, X=c(X[,2]), W=W, n=n, J=J)
                      
                      #*********************
                      # Initial values
                      #*********************
                      
                      zst<-apply(y,1,max)
                      
                      
                      inits <- function()
                      {
                        if ( mean((mle[1:2]-beta.param))>10 |  mean((mle[3:4]-alpha.param))>10 ){
                          
                          list(zb=zst,
                               beta0=0,
                               beta1=0,
                               alpha0=0,
                               alpha1=0)
                          
                        }else{
                          
                          list(zb=zst,
                               beta0=mle[1],
                               beta1=mle[2],
                               alpha0=mle[3],
                               alpha1=mle[4] )
                        }
                      }
                      
                      #****************************
                      # Paramaters to be monitored
                      #****************************
                      
                      parameters <- c("beta0", "beta1", "alpha0", "alpha1")#, "occupied")
                      
                      #print("start mcmc")
                      nsims<-20000
                      niter <- nsims
                      nburn <- floor(niter/2)
                      nthin <- 1
                      nchains <- 1
                      nkeep<-niter-nburn
                      
                      
                      t1<-proc.time()
                      occ.jags <- jags(occ.data, inits, parameters, "occ.txt",
                                       n.chains=1, n.iter=niter, n.burnin=nburn,
                                       parallel = FALSE, verbose=FALSE)
                      t2<-proc.time()
                      jagsTimer<-t2[3]-t1[3]
                      
                      jagsESS<-c(effectiveSize(as.mcmc(occ.jags$sims.list$alpha0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$alpha1)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta0)),
                                 effectiveSize(as.mcmc(occ.jags$sims.list$beta1)))
                      
                      #alpha0 = occ.jags$sims.list$alpha0
                      #alpha1 = occ.jags$sims.list$alpha1
                      #beta0 = occ.jags$sims.list$beta0
                      #beta1 = occ.jags$sims.list$beta1
                      
                      #jags_samples<-cbind(alpha0, alpha1, beta0, beta1)
                      rm(occ.jags)
                      #-------------------------------------------------------------------------------
                      
                      design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
                      
                      #set the formula to use
                      form1<- V1 ~ X1 ~ W1
                      
                      #the prior distributions specified here
                      alpha_m<-matrix( mle[3:4], ncol=1)
                      beta_m<-matrix( mle[1:2], ncol=1)
                      sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
                      sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
                      
                      
                      t1_drum5<-proc.time()
                      occdRUM<-dRUMocc(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p)
                      t2_drum5<-proc.time()
                      drumTimer<-t2_drum5[3]-t1_drum5[3]
                      
                      dRUMESS<-c(effectiveSize(as.mcmc(occdRUM$alpha[1,])),
                                 effectiveSize(as.mcmc(occdRUM$alpha[2,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[1,])),
                                 effectiveSize(as.mcmc(occdRUM$beta[2,])))
                      
                      rm(occdRUM)
                      #-------------------------------------------------------------------------------
                      
                      #Fit using Stan
                      N<-n
                      V<-J
                      Wmat<-W[,,2]
                      Xmat<-X[,2]
                      
                      Stan1<-proc.time()
                      occregfit <- sampling(Model_code_covocc, data=c("N","V","Xmat","Wmat","y"),
                                            iter = nsims, chains = 1, warmup = floor(nsims/2), thin = 1, verbose=F)
                      Stan2<-proc.time()
                      Stantimer<-Stan2[3]-Stan1[3]
                      
                      post_occregfit<-extract(occregfit, permuted = F)
                      StanESS<-effectiveSize(post_occregfit[,1,c(3,4,1,2)])
                      
                      rm(post_occregfit)
                      #-------------------------------------------------------------------------------
                      
                      c( jagsESS, dRUMESS, StanESS, jagsTimer, drumTimer, Stantimer)
                    }else{
                      c(rep(NA, 15))
                    }
                  }
stopCluster(cl)

#Now fit PG fit in series
outres<-matrix(0, nrow=iterations, 5)

for (iter in 1:iterations){
  
  #Simulate data
  
  set.seed(iter+up)
  
  #initial generation of the sample data
  x1 = runif(n, -2,2)
  x = (x1 - mean(x1)) / sd(x1)
  X  = cbind(rep(1,n), x)
  psi = as.vector(1/(1+exp(-X %*% beta.param))) ##logistic link function used
  z = rbinom(n, size=1, prob=psi)
  
  w1 = runif(n*J, -5,5)
  w = (w1 - mean(w1)) / sd(w1)
  w = matrix(w, nrow=n, ncol=J)
  W = array(dim=c(n,J,2))
  W[,,1] = 1
  W[,,2] = w
  
  p = matrix(nrow=n, ncol=J)
  y = matrix(nrow=n, ncol=J)
  for (j in 1:J) {
    p[, j] = c(1/(1+exp(-W[,j,] %*% alpha.param)))
    y[, j] = rbinom(n, size=1, prob=z*p[, j])
  }
  
  jind = !is.na(y)   # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  
  Y.eg<-y
  X.eg = as.data.frame(x) #siteCovs
  colnames(X.eg)<-c("X1")
  W1=matrix(NA,nrow=n, ncol=J)
  W1[1:n,]<-W[1:n,,2]
  W.eg.l1<-list(W1=W1)
  
  #-------------------------------------------------------------------------------
  
  #mle fit
  
  jind = !is.na(y)  # index for non-missing observations of y
  ysum = apply(y,1,sum, na.rm=TRUE)
  betaGuess = rep(0, dim(X)[2])
  alphaGuess = rep(0, dim(W)[3])
  paramGuess = c(betaGuess, alphaGuess)
  fit = optim(par=paramGuess, fn=negLL, method='BFGS', hessian=TRUE)
  
  mle = fit$par
  
  
  if (fit$convergence==0){
    
    #covariance matrix based on the mle fit
    
    vcv = chol2inv(chol(fit$hessian))
    
    #print("start mcmc")
    nsims<-20000
    
    design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=Y.eg)
    
    #set the formula to use
    form1<- V1 ~ X1 ~ W1
    
    #the prior distributions specified here
    alpha_m<-matrix( mle[3:4], ncol=1)
    beta_m<-matrix( mle[1:2], ncol=1)
    sigma_inv_beta_p<- solve(vcv[1:2, 1:2]) #prior inverse covariance for beta
    sigma_inv_alpha_p<- solve(vcv[3:4, 3:4]) #prior inverse covariance for alpha
    
    t1_PG<-proc.time()
    occPG<-PGocc4(formula=form1, design_mats=design_mats, ndraws=nsims, alpha_m, beta_m, sigma_inv_alpha_p, sigma_inv_beta_p, 0.5)
    t2_PG<-proc.time()
    PGTimer<-t2_PG[3]-t1_PG[3]
    
    PGESS<-c(effectiveSize(t(occPG[[1]])),
             effectiveSize(t(occPG[[2]])))
    
    rm(occPG)
    
    outres[iter,]<-c( PGESS, PGTimer)
  }else{
    outres[iter,]<-c(rep(NA, 5))
  }
}
output<-output12

output12<- cbind(output[,c(1:8)], outres[,c(1:4)], output[,c(9:12)], output[,c(13:14)], outres[,5], output[,15])
rm(output)

save.image("Logit_simulations_n500_J10_case4_EURING_revision.RData")
rm(output12) #remove output 12

cat("\n output12 DONE! \n")

#-------------------------------------------------------------------------------

