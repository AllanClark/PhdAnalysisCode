#----------FUNCTIONS TO BE USED -----------------------------------------------
# Define negative of log-likelihood function used to calculate MLE of parameter
negLL = function(param) 
{
     #The negative loglikelihood function assuming a probit link fucntion
     beta = param[1:dim(X)[2]]
     psi = pnorm(X %*% beta)
     alpha = param[(dim(X)[2]+1):(dim(X)[2]+dim(W)[3])]
     p = matrix(nrow=n, ncol=J)
     for (j in 1:J) {
          p[, j] = pnorm(W[,j,] %*% alpha)
     }
     logL = vector(length = n)
     for (i in 1:n) {
          yvec = y[i, jind[i,]]
          pvec = p[i, jind[i,]]
          terms = pvec^yvec * (1-pvec)^(1-yvec)
          logL[i] = log( psi[i] * prod(terms)  + ifelse(ysum[i]==0,1,0)*(1-psi[i])  )
     }
     (-1)*sum(logL)
}

#define the DA algorithm function used by Rodriguez and Dorazio
#rewrite the code to account for the fact that different locations might have been 
#sampled a different amount of times!
DA<-function(X, W, Y, J, ndraws=25000, alpha.starts, beta.starts,nburn=0.5*ndraws,thin=1)
{
     #the DA algorithm function used by Rodriguez and Dorazio - code slightly adjusted! 
     #X = occupancy design matrix. The constant should be included. Site specific covariates
     #W = detection design array. Detection specific covariates.
     #Y = the detection/no detection data
     #J = the number of site visits. Here I assume that it is the same for all sites.
     #X and W should not contain any missing values
     # alpha.starts, beta.starts - starting values of the MCMC iterations
     #ndraws = number of samples obtained (Note = burnin + post burn in)
     #half of the samples are discarded
     #nburn<- (ndraws/2)
     #thin is the interval at which the post burnin samples are thinned. default interval is 1
     # Use Gibbs algorithm to estimate summaries of the posterior distribution
     
     # .... initialize Gibbs sampler
     n<-NROW(X) #the number of sites
     XprimeX = t(X) %*% X
     # .... assign parameters of priors
     mu.beta = matrix(rep(0, dim(X)[2]), ncol=1)
     Sigma.beta = matrix(0, nrow=dim(X)[2], ncol=dim(X)[2])
     diag(Sigma.beta) = sqrt(1000)     #  use high value to specify non-informative prior
     SigmaInv.beta = chol2inv(chol(Sigma.beta))
     V.beta = chol2inv(chol(SigmaInv.beta + XprimeX))
     ScaledMu.beta = SigmaInv.beta %*% mu.beta
     
     mu.alpha = matrix(rep(0, dim(W)[3]), ncol=1)
     Sigma.alpha = matrix(0, nrow=dim(W)[3], ncol=dim(W)[3])
     diag(Sigma.alpha) = sqrt(1000)     #  use high value to specify non-informative prior
     SigmaInv.alpha = chol2inv(chol(Sigma.alpha))
     ScaledMu.alpha = SigmaInv.alpha %*% mu.alpha
     
     alpha = matrix(alpha.starts, ncol=1) #set random starts
     beta = matrix(beta.starts, ncol=1) #set random starts
     
     jind = !is.na(Y)   # index for non-missing observations of y
     ysum = apply(Y,1,sum, na.rm=TRUE)
     z = as.integer(ysum>0)
     
     # .... compute Gibbs draws
     post = matrix(nrow=ndraws, ncol=dim(X)[2] + dim(W)[3] + 1 + n)
     
     CR = '\n'
     cat('\n Begin Gibbs sampling:', CR, CR)
     for (draw in 1:ndraws) {
          if (draw == 0.5*ndraws)  cat('\n halfway done..... drawing sample #', draw, CR)
          
          #  draw z
          psi = as.vector(pnorm(X %*% beta))
          p = matrix(nrow=n, ncol=J)
          for (j in 1:J) {
               p[, j] = pnorm(W[,j,] %*% alpha)
          }
          
          z.prob<-matrix(1, nrow=n)
          for (i in 1:n) {
               if(ysum[i]==0) {
                    z.prob[i] = 1/ ( 1 +  (1 - psi[i] )/ (psi[i] * prod(1-p[i,jind[i,]])))
                    z[i] = rbinom(1, size=1, prob=z.prob[i])
               }
          }
          
          # draw v and beta 
          ind = z==1
          vmean = as.vector(X %*% beta)
          #v = rep(NA, dim(X)[1])
          v = vector(length = dim(X)[1])
          v[ind] = rtnorm(sum(ind), mean=vmean[ind], sd=1, lower=0)
          v[!ind] = rtnorm(sum(!ind), mean=vmean[!ind], sd=1, upper=0)
          
          betaMean =  V.beta %*% (ScaledMu.beta + (t(X) %*% v) )
          beta = matrix(mvrnorm(1, betaMean, V.beta), ncol=1)
          
          # draw u and alpha
          ind = z==1
          #  ... extract y values and w vectors of occupied sites, and arrange them in column-major order
          ymat = Y[ind, ]
          Warr = W[ind, ,]
          yind = jind[ind, ]
          yvec = as.vector(ymat[yind[,1], 1])
          Wmat = Warr[yind[,1] ,1, ]
          for (j in 2:J) {
               yvec = c(yvec, as.vector(ymat[yind[,j], j]))
               Wmat = rbind(Wmat, Warr[yind[,j] ,j, ])
          }
          umean = as.vector(Wmat %*% alpha)
          u = vector(length = dim(Wmat)[1])
          ind.y = yvec==1
          u[ind.y] = rtnorm(sum(ind.y), mean=umean[ind.y], sd=1, lower=0)
          u[!ind.y] = rtnorm(sum(!ind.y), mean=umean[!ind.y], sd=1, upper=0)
          
          WprimeW = t(Wmat) %*% Wmat
          V.alpha = chol2inv(chol(SigmaInv.alpha + WprimeW))
          
          alphaMean = V.alpha %*% (ScaledMu.alpha + (t(Wmat) %*% u) )
          alpha = matrix(mvrnorm(1, alphaMean, V.alpha), ncol=1)
          
          post[draw, ] = c(as.vector(beta), as.vector(alpha), mean(z), as.vector(z.prob))
     }
     cat('\n Gibbs sampling is completed!', CR, CR)
     cat("\n ....................................\n")
     
     #discard burnin samples4
     post.burn = post[-c(1:nburn),]
     #thin the postburnin samples
     thin.post = post.burn[seq(1,nrow(post.burn),by=thin),]
     #give names to the column variables
     colnames(thin.post) = c('beta_0', 'beta_1',
                             'alpha_0', 'alpha_1', "PAO", paste("occ",1:n, sep=""))
     return(thin.post)
}

vb_Designs_check<-function(formula, Names)
{
     #perform some checks on the design matrices
     #check that the names in the formula call are in the dataframes provided
     detVars<-all.vars(formula) #return names contained in formula
     if ( (sum(detVars[-1]%in% Names)==length(detVars[-1]))!= 1)
     {
          #if the condition '(sum(detVars[-1]%in% Names)==length(detVars[-1]))'
          #is false print the error below and terminate the function
          stop(print("\n CHECK YOUR FORMULA CALL. \n MISMATCH BETWEEN CALL AND 
                     DESIGN MATRICES. \n ."))
          #stop() 
     }
}

vb_ReqDesigns<-function(formula, design_mats)
{
     
     vb_Designs_check(formula, design_mats$Names)
     #create the W matrix
     W<-model.matrix(as.formula(paste("~",formula[3],sep="")), data=design_mats$W)
     #print(W)
     #create the X matrix
     f_xmat<-paste(formula[[2]])
     X<-model.matrix(as.formula(paste(f_xmat[1],f_xmat[3],sep="")), 
                     data=design_mats$X)
     #print(dim(X))
     list(W=W, X=X)
}

vb_Designs<-function(W, X, y)
{
     #create the required 'response' and 'regressor matrices'
     #using all of the X and W data
     #the output is stored as a named list
     #create the Y matrix that will be used
     
     #Y is a n*J times 1 matrix
     Y<-matrix(na.omit(matrix(t(y), ncol=1)))
     #col.names(Y)<-names(y)
     pres_abs <- apply(y,1,max,na.rm=T)
     #create the W matrix
     W.temp<-NULL
     nv<-length(W)
     for (i in 1:nv){W.temp<-cbind(W.temp, W[[i]])}
     #print(W.temp)
     #need to ask about this line below, doesnt seem to work
     nvisits<-apply(W[[1]],1,function(x){length(na.omit(x))})
     n<-length(nvisits)
     W.out<-NULL
     for (i in 1:n){W.out<-rbind(W.out, matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]) )}
     colnames(W.out)<-names(W)
     list(Y=as.data.frame(Y), X=as.data.frame(X), W=as.data.frame(W.out),
          Names=c( colnames(X), colnames(W.out)), nvisits=nvisits,
          pres_abs=pres_abs)
}

#VB code - Gaussian approximation
vb_model2_prob_la<-function(formula, design_mats, alpha_0, beta_0, 
                            Sigma_alpha_0, Sigma_beta_0, 
                            epsilon=1e-5, Print=FALSE)
{     
     #Date: 2 august 2016
     #Allan E Clark
     #Multi-visit site occupancy model - probit link functions are used
     #Estimation is undertaken using variational bayes and Laplace approximation
     #The assumption made is that the VB distributions of all regression 
     #coefficients are multivariate normal
     #--------------------------------------------------------------------------
     
     #--------------------------------------------------------------------------
     #Arguments
     #--------------------------------------------------------------------------
     #formula<- y~ occupancy covariates  ~ site detection covariates
     #y <- n by J matrix of presence absence data
     #n <- the number of locations
     #J <- the number of visits to the sites
     #Assumed that each site is visited J times
     
     #X <- a dataframe that contains the covariates used to calculate site 
     #occupancy probabilities
     #W <- a named list that contains the site covariates used to calculate 
     #the site detection probabilities
     #W has the same form as an unmarkedFrameOccu object in the 
     #unmarked package
     
     #alpha_0 <- starting value of detection covariate coefficients
     #beta_0 <- starting value of occurence covariate coefficients
     #Sigma_alpha_0, Sigma_beta_0 - the variance covariance matrix of 
     #alpha_0 and beta_0
     #LargeSample<-TRUE - indicates that the number of sites is 'large' 
     #and that an
     #approximation to B(mu, sigma^2) is used instead of integrations
     #---------------------------------------------------------------------
     
     #load the required functions and packages
     require(MASS) 
     require(spatstat)
     
     #Used to calculate <a(X%*%beta)>, <b(X%*%beta)> and <b(W%*%alpha)>
     a_x<-function(x){ pnorm(x, log.p=T)}
     b_x<-function(x){ pnorm(x, log.p=T, lower.tail = FALSE)}
     
     #load the required functions
     bx<-function(x){log(1+exp(x))}
     
     XWX.fun2 <- function(X,w)
     {
          #return( t(X)%*%diag(as.vector(w)))%*%X ) # Much slower
          #return( t(X*as.vector(w))%*%X  )         # Much faster
          return( crossprod(X*as.vector(w),X)  )    # faster than second version
     }
     
     diagSet <- function(d)
     {
          #John Ormerod code
          return(  d*((1:d)-1) + (1:d) )
     }
     
     diagElements <- function(A)
     {
          #John Ormerod code
          return( A[diagSet(nrow(A))] )
     }
     
     trace <- function(A)
     {
          #John Ormerod code
          #return(sum(diag(A)))          # Much slower
          return(sum(diagElements(A)))  # Much faster
     }
     
     X_row.Y_col=function(X,Y)
     {
          #does row by column matrix multiplication
          #ie row_i of matrix X %*% column_j of matrix Y
          #returns a column vector with the required elements
          #appears faster than a simple loop
          index=matrix(1:NROW(X))
          matrix(apply(index,1,function(x, m1=X,m2=Y){m1[x,]%*%m2[,x]}))
     }
     
     vb_Designs_check<-function(formula, Names)
     {
          #perform some checks on the design matrices
          #check that the names in the formula call are in the dataframes provided
          
          detVars<-all.vars(formula)
          
          if ( (sum(detVars[-1]%in% Names)==length(detVars[-1]))!= 1)
          {
               stop(print("\n \n CHECK YOUR FORMULA CALL. \n MISMATCH BETWEEN CALL AND DESIGN MATRICES. \n 
                          i.e. You included objects in the call: '~occupancy variables ~ detection variables' 
                          that does not appear in the design matrices."))
               #stop()
          }
     }
     
     vb_ReqDesigns<-function(formula, design_mats)
     {
          vb_Designs_check(formula, design_mats$Names)
          
          #create the W matrix
          W<-model.matrix(as.formula(paste("~",formula[3],sep="")), data=design_mats$W)
          #print(W)
          
          #create the X matrix
          f_xmat<-paste(formula[[2]])
          X<-model.matrix(as.formula(paste(f_xmat[1],f_xmat[3],sep="")), data=design_mats$X)
          #print(dim(X))
          
          list(W=W, X=X)
     }
     
     sum_Index_vec=function(Vec, starts_2_ends)
     {
          #The following function calulates of the sum of the elements of a vector over
          #various indices where the indices are stored in 'start_2_ends'
          #'start_2_ends' is a 'x by 2' matrix containing starting and
          #ending indices respectively.
          #eg.
          #xx=1:10 #Vec
          #indices<-cbind(c(1,3,7,8), c(2,6,7,10)) #start_2_ends
          #sum_Index_vec(xx, indices) #is the same as
          #matrix( c( sum(xx[1:2]), sum(xx[3:6]), sum(xx[7:7]), sum(xx[8:10]) ) )
          
          matrix(apply(starts_2_ends,1,function(x,y=Vec){sum(Vec[x[1]:x[2]])}))
     }
     
     FUN_Index_vec=function(Vec, starts_2_ends,FUN=sum)
     {
          #The following function applies a function to the elements of a vector over
          #various indices where the indices are stored in 'start_2_ends'
          #'start_2_ends' is a 'x by 2' matrix containing starting and
          #ending indices respectively.
          #here FUN only has one argument!
          #eg.
          #xx=1:10 #Vec
          #indices<-cbind(c(1,3,7,8), c(2,6,7,10)) #start_2_ends
          #FUN_Index_vec(xx, indices,sum) #is the same as
          #matrix( c( sum(xx[1:2]), sum(xx[3:6]), sum(xx[7:7]), sum(xx[8:10]) ) )
          
          matrix(apply(starts_2_ends,1,function(x,y=Vec){eval(FUN(Vec[x[1]:x[2]]))}))
     }
     
     # #At a later stage you could possibly investigate whether numerical integrations work better!
     # exp_a<-function(row_i, nsims=1000, mu, Sigma)
     # {
     #      #row_i is the row vector of a design matrix
     #      #mu is the mean equation of the vb distribution used
     #      #Sigma is the covariance matrix of the vb distribution used
     #  
     #      #the 1000 could be made general! number of bootstrap samples used
     #      #the random variable sampled
     #      rv<-mvrnorm(nsims, mu=mu, Sigma=Sigma)
     #  
     #  #a check that the function is correct!
     #  #jj<-NULL
     #  #for (i in 1:nsims)
     #  #{
     #  #    print(rv[i,])
     #  #    jj[i]<-pnorm(row_i%*%rv[i,], lower.tail=TRUE, log.p=TRUE)
     #  #}
     #  #print(jj)
     #  #print(mean(jj)) - this should be the same as the answer below
     #  
     #      return( mean(apply(rv,1, function(x, y=row_i) pnorm(y%*%x, lower.tail=TRUE, log.p=TRUE))) )
     # }
     # 
     # exp_b<-function(row_i, nsims=1000, mu, Sigma)
     # {
     #      #row_i is the row vector of a design matrix
     #      #mu is the mean equation of the vb distribution used
     #      #Sigma is the covariance matrix of the vb distribution used
     #  
     #      #the 1000 could be made general! number of bootstrap samples used
     #      #the random variable sampled
     #      rv<-mvrnorm(nsims, mu=mu, Sigma=Sigma)
     #  
     #  #a check that the function is correct!
     #  #jj<-NULL
     #  #for (i in 1:nsims)
     #  #{
     #  #    print(rv[i,])
     #  #    jj[i]<-pnorm(row_i%*%rv[i,], lower.tail=TRUE, log.p=TRUE)
     #  #}
     #  #print(jj)
     #  #print(mean(jj)) - this should be the same as the answer below
     #  
     #      return( mean(apply(rv,1, function(x, y=row_i) pnorm(y%*%x, lower.tail=FALSE, log.p=TRUE))) )
     # }
     
     #----------------------------------------------------------------------------------
     
     log_q_alpha<-function(par)
     {
          alpha<-matrix(par, ncol=1)
          alpha_diff<-alpha - alpha_0
          
          mu_alpha<-W_vb%*%alpha
          #sigma_alpha<-sqrt(diagElements(W_vb%*%Sigma_alpha%*%t(W_vb))) #slightly faster (consider loops here)
          
          #mu_alpha.sigma_alpha<-cbind(mu_alpha,sigma_alpha)
          a_Wa<-a_x(mu_alpha)
          b_Wa<-b_x(mu_alpha)
          
          term1<- -crossprod(alpha_diff , SigmaInv_alpha_0)%*%alpha_diff/2
          term2<- t(Y)%*%P_tilde%*%( a_Wa - b_Wa ) + sum( P_tilde%*%b_Wa )
          
          return( -1*( term1 + term2 ) ) #multiply by -1 since we want to maximise the function
     }
     
     log_q_beta<-function(par)
     {
          beta<-matrix(par, ncol=1)
          beta_diff<-beta - beta_0
          
          mu_beta<-X%*%beta
          #sigma_beta<-sqrt(diagElements(X%*%Sigma_beta%*%t(X))) #slightly faster (consider loops here)
          
          #mu_beta.sigma_beta<-cbind(mu_beta,sigma_beta)
          a_Xb<-a_x(mu_beta)
          b_Xb<- b_x(mu_beta)
          
          term1<- -crossprod(beta_diff , SigmaInv_beta_0)%*%beta_diff/2
          term2<- sum( p*( a_Xb - b_Xb ) ) + sum( b_Xb )
          
          return( -1*( term1 + term2 ) ) #multiply by -1 since we want to maximise the function
     }
     #==========================================================================================
     
     
     #Declare certain matrices and constants
     #---------------------------------------
     #req_design_mats<-vb_ReqDesigns(formula=form1, design_mats)
     req_design_mats<-vb_ReqDesigns(formula, design_mats)
     W_vb<-req_design_mats$W
     X<-req_design_mats$X
     Y<-matrix(design_mats$Y$V1)
     
     #a matrix containing how many times each of the sites are visited
     #J<-dim(W)[2]
     V<-design_mats$nvisits #nvisits
     
     n<-NROW(X)
     N<-length(Y)
     pn<-NCOL(X) #number of columns of X
     qn<-NCOL(W_vb) #number of columns of W
     
     const1 <- (pn+qn)*(log(2*pi)+1) #+1 added since we require it for the calculation of the entropy of a mult norm
     
     const2<-determinant(Sigma_alpha_0, logarithm=T)$modulus[1] + determinant(Sigma_beta_0, logarithm=T)$modulus[1]
     const3<- 0.5*(const1+const2)
     
     SigmaInv_alpha_0 <- chol2inv(chol(Sigma_alpha_0))
     SigmaInv_beta_0 <- chol2inv(chol(Sigma_beta_0))
     SigmaInv_times_alpha_0 <- SigmaInv_alpha_0%*%alpha_0
     SigmaInv_times_beta_0 <- SigmaInv_beta_0%*%beta_0
     #---------------------------------------
     
     #the prior mean and covariance matrix of alpha and beta
     #------------------------------------------------------
     #the starting value for the alpha and beta vector (and covariance matrices)
     #--------------------------------------------------------------------------
     
     alpha <- alpha_0
     beta <- beta_0 
     alpha_0<-alpha
     beta_0<-beta
     Sigma_alpha <- Sigma_alpha_0
     Sigma_beta <- Sigma_beta_0
     #Sigma_alpha_inv <-solve(Sigma_alpha_0)
     #Sigma_beta_inv <- solve(Sigma_beta_0)
     
     oldparams<- c(alpha, Sigma_alpha, beta, Sigma_beta)
     nparams <- length(oldparams)
     
     #indentify which at which of the sites the species were observed!
     #----------------------------------------------------------------
     pres_abs<-design_mats$pres_abs
     obs_vec <- which(pres_abs==1) #where they were observed
     not_obs_vec <- which(pres_abs==0) #where they were not observed
     n_obs <- length(obs_vec)
     n_not_obs <- length(not_obs_vec)
     V_not_obs<- V[not_obs_vec] #the nvistis where the species has not been observed
     
     #initial starting values of E(Z_tilde) and E(Z)
     #----------------------------------------------
     
     #calculate <a(X%*%beta)>, <b(X%*%beta)>, <a(W%*%alpha)> and <b(W%*%alpha)>  for all sites
     mu_beta<-X%*%beta
     sigma_beta<-sqrt(diagElements(X%*%Sigma_beta%*%t(X))) #find quick way of doing this!
     mu_beta.sigma_beta<-cbind(mu_beta,sigma_beta)
     Exp.a_Xb<-matrix(apply(mu_beta.sigma_beta,1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     Exp.b_Xb<-matrix(apply(mu_beta.sigma_beta,1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     
     mu_alpha<-W_vb%*%alpha
     #sigma_alpha<-sqrt(diag(W_vb%*%Sigma_alpha%*%t(W_vb))) #find quick way of doing this!
     sigma_alpha<-sqrt(diagElements(W_vb%*%Sigma_alpha%*%t(W_vb))) #slightly faster (consider loops here)
     
     mu_alpha.sigma_alpha<-cbind(mu_alpha,sigma_alpha)
     Exp.b_Wa<-matrix(apply(mu_alpha.sigma_alpha,1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     Exp.a_Wa<-matrix(apply(mu_alpha.sigma_alpha,1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     
     #if the species is observed at a site the element should be 1
     p_tilde <- matrix(1, ncol=1, nrow=N)
     p<-matrix(1, ncol=1, nrow=n)
     
     #the index of the location of the starting values of Z_tilde
     #Z_tilde <- [(z1,z1,...z1),....,(zn,.......,zn)]^T
     #only keep the ones where the species has not been observed
     starts_z <- matrix(c(1+cumsum(V)-V)[not_obs_vec] , ncol=1)
     ends_z<- starts_z+V_not_obs-1
     starts_2_ends_z<-cbind(starts_z, ends_z)
     
     #Indices<-apply(starts_z,1,function(x,J.=3) x:(x+J.-1))
     #c_Indices<-c(Indices)
     #vectorize this!
     c_Indices<-NULL
     for (c_i in 1:n_not_obs){c_Indices=c(c_Indices, starts_z[c_i]:ends_z[c_i])}
     starts_notz <- matrix(c(1+cumsum(V_not_obs)-V_not_obs) , ncol=1)
     ends_notz<- starts_notz+V_not_obs-1
     starts_2_ends_notz<-cbind(starts_notz, ends_notz)
     
     #it2<-apply(X, 1, FUN=exp_a, nsims=50000, mu=beta, Sigma=Sigma_beta) #does using sampling
     #cbind(it1, it2)
     #Let dii be the i-th diagonal element, then
     #dii = (sum over j from 1 to b)(sum over k from 1 to b) [(xij)(xik)(bjk)]
     #where i = 1,, a.
     
     #term1<-apply(X[not_obs_vec,], 1, FUN=exp_a, nsims=1000, mu=beta, Sigma=Sigma_beta)
     #term2<-apply(X[not_obs_vec,], 1, FUN=exp_b, nsims=1000, mu=beta, Sigma=Sigma_beta)
     #term3<-sum_Index_vec( apply(W_vb[c_Indices,], 1, FUN=exp_b, nsims=100000, 
     #                            mu=alpha, Sigma=Sigma_alpha) , starts_2_ends_notz)
     
     #Some checks regarding the gauss hermite quadrature procedure and the bootstrap sampling method
     #---------------------------------------------------------------------------
     
     #the bootstrap and the gauss hermite solutions are not exactly the same
     #but they are close enough!
     #apply(X[1:5,], 1, FUN=exp_a, nsims=100000, mu=beta, Sigma=matrix(c(1,0.5,0.5,1),2,2))
     #matrix(apply(cbind(mu_beta, sqrt(diag(X%*%matrix(c(1,.5,0.5,1),2,2)%*%t(X))) ),1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=5)}), ncol=1)[1:5]
     
     #the differences are larger when the means and the variances are large - but they are still quite close
     #apply(X[1:5,], 1, FUN=exp_a, nsims=100000, mu=beta, Sigma=matrix(c(20,0.5,0.5,50),2,2))
     #matrix(apply(cbind(mu_beta, sqrt(diag(X%*%matrix(c(20,.5,0.5,50),2,2)%*%t(X))) ),1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=5)}), ncol=1)[1:5]
     
     #apply(X[1:5,], 1, FUN=exp_a, nsims=100000, mu=beta*10, Sigma=matrix(c(20,0.5,0.5,50),2,2))
     #matrix(apply(cbind(mu_beta*10, sqrt(diag(X%*%matrix(c(20,.5,0.5,50),2,2)%*%t(X))) ),1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=5)}), ncol=1)[1:5]
     
     #apply(X[1:5,], 1, FUN=exp_b, nsims=100000, mu=beta, Sigma=matrix(c(1,0.5,0.5,1),2,2))
     #matrix(apply(cbind(mu_beta, sqrt(diag(X%*%matrix(c(1,0.5,0.5,1),2,2)%*%t(X))) ),1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=5)}), ncol=1)[1:5]
     
     #c3<-sum_Index_vec( apply(W_vb[c_Indices,], 1, FUN=exp_b, nsims=50000, mu=alpha, Sigma=matrix(c(1,0.5,0.5,1),2,2)) , starts_2_ends_notz)
     #c(c3)
     
     #c(sum_Index_vec(
     #               matrix(apply(cbind(mu_alpha, sqrt(diagElements(W_vb%*%matrix(c(1,0.5,0.5,1),2,2)%*%t(W_vb)))
     #                   ),1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=50)})
     #                   , ncol=1)[c_Indices], starts_2_ends_notz))
     
     #Have a look at the following to compute multivariate Gauss hermite integrals
     #https://www.r-bloggers.com/notes-on-multivariate-gaussian-quadrature-with-r-code/
     #---------------------------------------------------------------------------
     
     term1<-Exp.a_Xb[not_obs_vec]
     term2<-Exp.b_Xb[not_obs_vec]
     term3<-sum_Index_vec( Exp.b_Wa[c_Indices], starts_2_ends_notz)
     ti<-term1 - term2 + term3
     
     p[not_obs_vec] <- 1/(1+ exp(-ti) )
     #the group sizes are all different!
     p_tilde[c_Indices]<-rep( p[not_obs_vec], times=V_not_obs)
     P_tilde<-diag(c(p_tilde))
     
     #change this below!!!
     
     #Marginal likelihood approximation <- E(Log_p(y,z,alpha,beta))- E(q(alpha, beta, z))
     #-----------------------------------------------------------------------------------
     
     Logp<-sum( p*(Exp.a_Xb - Exp.b_Xb) ) + sum( Exp.b_Xb ) + t(Y)%*%P_tilde%*%(Exp.a_Wa - Exp.b_Wa)
     
     #<ln( p(alpha, beta) )> when one is integrating wrt the optimal distribution!
     alpha_diff<-alpha - alpha_0
     beta_diff<-beta - beta_0
     f1<- trace(SigmaInv_alpha_0%*%Sigma_alpha) + crossprod(alpha_diff , SigmaInv_alpha_0)%*%alpha_diff
     f2<- trace(SigmaInv_beta_0%*%Sigma_beta) + crossprod(beta_diff , SigmaInv_beta_0)%*%beta_diff
     Logp<-Logp -0.5*( (pn+qn)*log(2*pi) - determinant(Sigma_alpha_0, logarithm=T)$modulus[1] - determinant(Sigma_beta_0, logarithm=T)$modulus[1] + f1 + f2)
     
     #the E(q(alpha, beta, z)) part
     #------------------------------
     
     #E_log_pz <- sum(p[not_obs_vec]*log(p[not_obs_vec]/(1-p[not_obs_vec])) + log(1-p[not_obs_vec])) #gave numerical errors sometimes
     E_log_pz <-sum(log(p[not_obs_vec]^p[not_obs_vec]) + log( (1-p[not_obs_vec])^(1-p[not_obs_vec]) ))
     E_log_p_alpha_beta <- -0.5*( const1 + determinant(Sigma_alpha, logarithm=T)$modulus[1] + determinant(Sigma_beta, logarithm=T)$modulus[1] )
     Log_ml <- Logp - E_log_p_alpha_beta[1] - E_log_pz
     
     #cat("------------------------------------------------------------------------------------","\n")
     #cat("STARTING VALUES","\n")
     #cat("Log_ml <- ", Log_ml,", alpha <- ", round(alpha,digits<-3) , ", beta <- ", round(beta,digits<-3), "\n")
     #cat("------------------------------------------------------------------------------------","\n")
     
     old_Log_ml<-Log_ml
     
     diff_Log_ml<-1
     its<-0
     #epsilon <- 1e-10
     Breakcounter<-0
     
     while (diff_Log_ml > epsilon)
     {
          its<-1+its
          cat(" VB - Iteration #: ", its,";")
          #cat("------------------------ \n")
          
          if (its == 100)
          {
               #print("Iterations large")
               Breakcounter<-1
               stop("Iterations large")
          }
          
          #Perform the Newton Rhaphson algortihm to calculate alpha and beta at iteration t
          #---------------------------------------------------------------------------------
          #crit <- .5 #used to assess convergence of the parameters - Rule 1
          crit<-0 #used for Rule 3
          
          alpha.optim<-optim(alpha, log_q_alpha, method="BFGS", hessian=TRUE)
          alpha<-alpha.optim$par
          Sigma_alpha<- solve(alpha.optim$hessian)
          
          beta.optim<-optim(beta, log_q_beta, method="BFGS", hessian=TRUE)
          beta<-beta.optim$par
          Sigma_beta<- solve(beta.optim$hessian)
          
          #initial starting values of E(Z_tilde) and E(Z)
          #----------------------------------------------
          
          #calculate <a(X%*%beta)>, <b(X%*%beta)>, <a(W%*%alpha)> and <b(W%*%alpha)>  for all sites
          mu_beta<-X%*%beta
          sigma_beta<-sqrt(diagElements(X%*%Sigma_beta%*%t(X))) #find quick way of doing this!
          mu_beta.sigma_beta<-cbind(mu_beta,sigma_beta)
          Exp.a_Xb<-matrix(apply(mu_beta.sigma_beta,1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
          Exp.b_Xb<-matrix(apply(mu_beta.sigma_beta,1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
          
          mu_alpha<-W_vb%*%alpha
          #sigma_alpha<-sqrt(diag(W_vb%*%Sigma_alpha%*%t(W_vb))) #find quick way of doing this!
          sigma_alpha<-sqrt(diagElements(W_vb%*%Sigma_alpha%*%t(W_vb))) #slightly faster (consider loops here)
          
          mu_alpha.sigma_alpha<-cbind(mu_alpha,sigma_alpha)
          Exp.b_Wa<-matrix(apply(mu_alpha.sigma_alpha,1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
          Exp.a_Wa<-matrix(apply(mu_alpha.sigma_alpha,1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
          
          #if the species is observed at a site the element should be 1
          p_tilde <- matrix(1, ncol=1, nrow=N)
          p<-matrix(1, ncol=1, nrow=n)
          
          term1<-Exp.a_Xb[not_obs_vec]
          term2<-Exp.b_Xb[not_obs_vec]
          term3<-sum_Index_vec( Exp.b_Wa[c_Indices], starts_2_ends_notz)
          ti<-term1 - term2 + term3
          
          p[not_obs_vec] <- 1/(1+ exp(-ti) )
          #the group sizes are all different!
          p_tilde[c_Indices]<-rep( p[not_obs_vec], times=V_not_obs)
          P_tilde<-diag(c(p_tilde))
          
          #Marginal likelihood approximation <- E(Log_p(y,z,alpha,beta))- E(q(alpha, beta, z))
          #-----------------------------------------------------------------------------------
          
          Logp<-sum( p*(Exp.a_Xb - Exp.b_Xb) ) + sum( Exp.b_Xb ) + t(Y)%*%P_tilde%*%(Exp.a_Wa - Exp.b_Wa)
          
          #<ln( p(alpha, beta) )> when one is integrating wrt the optimal distribution!
          alpha_diff<-alpha - alpha_0
          beta_diff<-beta - beta_0
          f1<- trace(SigmaInv_alpha_0%*%Sigma_alpha) + crossprod(alpha_diff , SigmaInv_alpha_0)%*%alpha_diff
          f2<- trace(SigmaInv_beta_0%*%Sigma_beta) + crossprod(beta_diff , SigmaInv_beta_0)%*%beta_diff
          Logp<-Logp -0.5*( (pn+qn)*log(2*pi) - determinant(Sigma_alpha_0, logarithm=T)$modulus[1] - determinant(Sigma_beta_0, logarithm=T)$modulus[1] + f1 + f2)
          
          #the E(q(alpha, beta, z)) part
          #------------------------------
          
          #E_log_pz <- sum(p[not_obs_vec]*log(p[not_obs_vec]/(1-p[not_obs_vec])) + log(1-p[not_obs_vec])) #gave numerical errors sometimes
          E_log_pz <-sum(log(p[not_obs_vec]^p[not_obs_vec]) + log( (1-p[not_obs_vec])^(1-p[not_obs_vec]) ))
          E_log_p_alpha_beta <- -0.5*( const1 + determinant(Sigma_alpha, logarithm=T)$modulus[1] + determinant(Sigma_beta, logarithm=T)$modulus[1] )
          Log_ml <- Logp - E_log_p_alpha_beta[1] - E_log_pz
          
          diff_Log_ml <- abs(Log_ml-old_Log_ml)
          old_Log_ml <- Log_ml
          
          #if (Print==TRUE)
          #{cat("Iteration: ", its,", Log_ml <- ", round(Log_ml,digits<-6), ", alpha <- ", round(alpha,digits<-3) , ", beta <- ", round(beta,digits<-3), "\n")}
     }
     
     list(alpha=alpha, beta=beta, Sigma_alpha=Sigma_alpha, Sigma_beta=Sigma_beta, occup_p=c(p), Log_mla=Log_ml, Breakcounter=Breakcounter)
     } 

#VB code - mixture of Gaussian approximation
vb_model2_gm<-function(formula, design_mats, alpha_0, beta_0, 
                       Sigma_alpha_0, Sigma_beta_0, epsilon=1e-3)
{     
     #Date: 25 October 2016
     #Allan E Clark
     #the following function goes about fitting a Gaussian mixture distribution
     #to the VB optimal distributions of alpha and beta
     
     #Multi-visit site occupancy model - probit link functions are used
     #Estimation is undertaken using variational bayes and Laplace approximation
     #The assumption made is that the VB distributions of all regression 
     #coefficients are multivariate normal
     #--------------------------------------------------------------------------
     
     #--------------------------------------------------------------------------
     #Arguments
     #--------------------------------------------------------------------------
     #formula<- y~ occupancy covariates  ~ site detection covariates
     #y <- n by J matrix of presence absence data
     #n <- the number of locations
     #J <- the number of visits to the sites
     #Assumed that each site is visited J times
     
     #X <- a dataframe that contains the covariates used to calculate site 
     #occupancy probabilities
     #W <- a named list that contains the site covariates used to calculate 
     #the site detection probabilities
     #W has the same form as an unmarkedFrameOccu object in the 
     #unmarked package
     
     #alpha_0 <- starting value of detection covariate coefficients
     #beta_0 <- starting value of occurence covariate coefficients
     #Sigma_alpha_0, Sigma_beta_0 - the variance covariance matrix of 
     #alpha_0 and beta_0
     #LargeSample<-TRUE - indicates that the number of sites is 'large' 
     #and that an
     #approximation to B(mu, sigma^2) is used instead of integrations
     #---------------------------------------------------------------------
     
     #load the required functions and packages
     require(MASS) 
     require(spatstat)
     require(iterLap)
     require(mvtnorm)
     
     #create a vectorized Gauss hermite function
     #the arguments mu and sd are vectorized
     #not that these arguments are vectors and not matrices
     v.gauss.hermite<-Vectorize(gauss.hermite, c("mu", "sd"))
     
     #Used to calculate <a(X%*%beta)>, <b(X%*%beta)> and <b(W%*%alpha)>
     a_x<-function(x){ pnorm(x, log.p=T)}
     b_x<-function(x){ pnorm(x, log.p=T, lower.tail = FALSE)}
     
     #load the required functions
     bx<-function(x){log(1+exp(x))}
     
     XWX.fun2 <- function(X,w)
     {
          #return( t(X)%*%diag(as.vector(w)))%*%X ) # Much slower
          #return( t(X*as.vector(w))%*%X  )         # Much faster
          return( crossprod(X*as.vector(w),X)  )    # faster than second version
     }
     
     diagSet <- function(d)
     {
          #John Ormerod code
          return(  d*((1:d)-1) + (1:d) )
     }
     
     diagElements <- function(A)
     {
          #John Ormerod code
          return( A[diagSet(nrow(A))] )
     }
     
     trace <- function(A)
     {
          #John Ormerod code
          #return(sum(diag(A)))          # Much slower
          return(sum(diagElements(A)))  # Much faster
     }
     
     X_row.Y_col=function(X,Y)
     {
          #does row by column matrix multiplication
          #ie row_i of matrix X %*% column_j of matrix Y
          #returns a column vector with the required elements
          #appears faster than a simple loop
          index=matrix(1:NROW(X))
          matrix(apply(index,1,function(x, m1=X,m2=Y){m1[x,]%*%m2[,x]}))
     }
     
     vb_Designs_check<-function(formula, Names)
     {
          #perform some checks on the design matrices
          #check that the names in the formula call are in the dataframes provided
          
          detVars<-all.vars(formula)
          
          if ( (sum(detVars[-1]%in% Names)==length(detVars[-1]))!= 1)
          {
               stop(print("\n \n CHECK YOUR FORMULA CALL. \n MISMATCH BETWEEN CALL AND DESIGN MATRICES. \n 
                          i.e. You included objects in the call: '~occupancy variables ~ detection variables' 
                          that does not appear in the design matrices."))
               #stop()
          }
     }
     
     vb_ReqDesigns<-function(formula, design_mats)
     {
          vb_Designs_check(formula, design_mats$Names)
          
          #create the W matrix
          W<-model.matrix(as.formula(paste("~",formula[3],sep="")), data=design_mats$W)
          #print(W)
          
          #create the X matrix
          f_xmat<-paste(formula[[2]])
          X<-model.matrix(as.formula(paste(f_xmat[1],f_xmat[3],sep="")), data=design_mats$X)
          #print(dim(X))
          
          list(W=W, X=X)
     }
     
     sum_Index_vec=function(Vec, starts_2_ends)
     {
          #The following function calulates of the sum of the elements of a vector over
          #various indices where the indices are stored in 'start_2_ends'
          #'start_2_ends' is a 'x by 2' matrix containing starting and
          #ending indices respectively.
          #eg.
          #xx=1:10 #Vec
          #indices<-cbind(c(1,3,7,8), c(2,6,7,10)) #start_2_ends
          #sum_Index_vec(xx, indices) #is the same as
          #matrix( c( sum(xx[1:2]), sum(xx[3:6]), sum(xx[7:7]), sum(xx[8:10]) ) )
          
          matrix(apply(starts_2_ends,1,function(x,y=Vec){sum(Vec[x[1]:x[2]])}))
     }
     
     FUN_Index_vec=function(Vec, starts_2_ends,FUN=sum)
     {
          #The following function applies a function to the elements of a vector over
          #various indices where the indices are stored in 'start_2_ends'
          #'start_2_ends' is a 'x by 2' matrix containing starting and
          #ending indices respectively.
          #here FUN only has one argument!
          #eg.
          #xx=1:10 #Vec
          #indices<-cbind(c(1,3,7,8), c(2,6,7,10)) #start_2_ends
          #FUN_Index_vec(xx, indices,sum) #is the same as
          #matrix( c( sum(xx[1:2]), sum(xx[3:6]), sum(xx[7:7]), sum(xx[8:10]) ) )
          
          matrix(apply(starts_2_ends,1,function(x,y=Vec){eval(FUN(Vec[x[1]:x[2]]))}))
     }
     
     log_q_alpha<-function(par)
     {
          
          alpha<-matrix(par, ncol=1)
          alpha_diff<-alpha - alpha_0
          
          mu_alpha<-W_vb%*%alpha
          #sigma_alpha<-sqrt(diagElements(W_vb%*%Sigma_alpha%*%t(W_vb))) #slightly faster (consider loops here)
          
          #mu_alpha.sigma_alpha<-cbind(mu_alpha,sigma_alpha)
          a_Wa<-a_x(mu_alpha)
          b_Wa<-b_x(mu_alpha)
          
          term1<- -crossprod(alpha_diff , SigmaInv_alpha_0)%*%alpha_diff/2
          term2<- t(Y)%*%P_tilde%*%( a_Wa - b_Wa ) + sum( P_tilde%*%b_Wa )
          
          
          return( ( term1 + term2 ) ) 
     }
     
     log_q_beta<-function(par)
     {
          beta<-matrix(par, ncol=1)
          beta_diff<-beta - beta_0
          
          mu_beta<-X%*%beta
          #sigma_beta<-sqrt(diagElements(X%*%Sigma_beta%*%t(X))) #slightly faster (consider loops here)
          
          #mu_beta.sigma_beta<-cbind(mu_beta,sigma_beta)
          a_Xb<-a_x(mu_beta)
          b_Xb<- b_x(mu_beta)
          
          term1<- -crossprod(beta_diff , SigmaInv_beta_0)%*%beta_diff/2
          term2<- sum( p*( a_Xb - b_Xb ) ) + sum( b_Xb )
          
          return( ( term1 + term2 ) ) 
     }
     
     # dmvnorm<-function(x, mu, Sigma)
     # {
     #      #the pdf of a multivariate normal distribution
     #      M<-length(mu)
     #      t1<-det(Sigma)^(-.5)
     #      t2<-(2*pi)^(M/2)
     #      t3<-exp(-.5*t(x-mu)%*%solve(Sigma)%*%(x-mu))
     #      return(t1*t2*t3)
     # }
     
     #==========================================================================================
     
     
     #Declare certain matrices and constants
     #---------------------------------------
     #req_design_mats<-vb_ReqDesigns(formula=form1, design_mats)
     req_design_mats<-vb_ReqDesigns(formula, design_mats)
     W_vb<-req_design_mats$W
     X<-req_design_mats$X
     Y<-matrix(design_mats$Y$V1)
     
     #a matrix containing how many times each of the sites are visited
     #J<-dim(W)[2]
     V<-design_mats$nvisits #nvisits
     
     n<-NROW(X)
     N<-length(Y)
     pn<-NCOL(X) #number of columns of X
     qn<-NCOL(W_vb) #number of columns of W
     
     const1 <- (pn+qn)*(log(2*pi)+1) #+1 added since we require it for the calculation of the entropy of a mult norm
     
     const2<-determinant(Sigma_alpha_0, logarithm=T)$modulus[1] + determinant(Sigma_beta_0, logarithm=T)$modulus[1]
     const3<- 0.5*(const1+const2)
     
     SigmaInv_alpha_0 <- chol2inv(chol(Sigma_alpha_0))
     SigmaInv_beta_0 <- chol2inv(chol(Sigma_beta_0))
     SigmaInv_times_alpha_0 <- SigmaInv_alpha_0%*%alpha_0
     SigmaInv_times_beta_0 <- SigmaInv_beta_0%*%beta_0
     
     #print(alpha_0)
     const4<- t(alpha_0)%*%SigmaInv_alpha_0
     const5<- t(beta_0)%*%SigmaInv_beta_0
     #---------------------------------------
     
     #the prior mean and covariance matrix of alpha and beta
     #------------------------------------------------------
     #the starting value for the alpha and beta vector (and covariance matrices)
     #--------------------------------------------------------------------------
     
     alpha <- alpha_0
     beta <- beta_0 
     alpha_0<-alpha
     beta_0<-beta
     Sigma_alpha <- Sigma_alpha_0
     Sigma_beta <- Sigma_beta_0
     #Sigma_alpha_inv <-solve(Sigma_alpha_0)
     #Sigma_beta_inv <- solve(Sigma_beta_0)
     
     oldparams<- c(alpha, Sigma_alpha, beta, Sigma_beta)
     nparams <- length(oldparams)
     
     #indentify which at which of the sites the species were observed!
     #----------------------------------------------------------------
     pres_abs<-design_mats$pres_abs
     obs_vec <- which(pres_abs==1) #where they were observed
     not_obs_vec <- which(pres_abs==0) #where they were not observed
     n_obs <- length(obs_vec)
     n_not_obs <- length(not_obs_vec)
     V_not_obs<- V[not_obs_vec] #the nvistis where the species has not been observed
     
     #initial starting values of E(Z_tilde) and E(Z)
     #----------------------------------------------
     
     #calculate <a(X%*%beta)>, <b(X%*%beta)>, <a(W%*%alpha)> and <b(W%*%alpha)>  for all sites
     mu_beta<-X%*%beta
     sigma_beta<-sqrt(diagElements(X%*%Sigma_beta%*%t(X))) #find quick way of doing this!
     mu_beta.sigma_beta<-cbind(mu_beta,sigma_beta)
     Exp.a_Xb<-matrix(apply(mu_beta.sigma_beta,1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     Exp.b_Xb<-matrix(apply(mu_beta.sigma_beta,1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     
     mu_alpha<-W_vb%*%alpha
     #sigma_alpha<-sqrt(diag(W_vb%*%Sigma_alpha%*%t(W_vb))) #find quick way of doing this!
     sigma_alpha<-sqrt(diagElements(W_vb%*%Sigma_alpha%*%t(W_vb))) #slightly faster (consider loops here)
     
     mu_alpha.sigma_alpha<-cbind(mu_alpha,sigma_alpha)
     Exp.b_Wa<-matrix(apply(mu_alpha.sigma_alpha,1, function(x){gauss.hermite(b_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     Exp.a_Wa<-matrix(apply(mu_alpha.sigma_alpha,1, function(x){gauss.hermite(a_x, mu=x[1], sd=x[2],order=25)}), ncol=1)
     
     #if the species is observed at a site the element should be 1
     p_tilde <- matrix(1, ncol=1, nrow=N)
     p<-matrix(1, ncol=1, nrow=n)
     
     #the index of the location of the starting values of Z_tilde
     #Z_tilde <- [(z1,z1,...z1),....,(zn,.......,zn)]^T
     #only keep the ones where the species has not been observed
     starts_z <- matrix(c(1+cumsum(V)-V)[not_obs_vec] , ncol=1)
     ends_z<- starts_z+V_not_obs-1
     starts_2_ends_z<-cbind(starts_z, ends_z)
     
     #Indices<-apply(starts_z,1,function(x,J.=3) x:(x+J.-1))
     #c_Indices<-c(Indices)
     #vectorize this!
     c_Indices<-NULL
     for (c_i in 1:n_not_obs){c_Indices=c(c_Indices, starts_z[c_i]:ends_z[c_i])}
     starts_notz <- matrix(c(1+cumsum(V_not_obs)-V_not_obs) , ncol=1)
     ends_notz<- starts_notz+V_not_obs-1
     starts_2_ends_notz<-cbind(starts_notz, ends_notz)
     
     term1<-Exp.a_Xb[not_obs_vec]
     term2<-Exp.b_Xb[not_obs_vec]
     term3<-sum_Index_vec( Exp.b_Wa[c_Indices], starts_2_ends_notz)
     ti<-term1 - term2 + term3
     
     p[not_obs_vec] <- 1/(1+ exp(-ti) )
     #the group sizes are all different!
     p_tilde[c_Indices]<-rep( p[not_obs_vec], times=V_not_obs)
     P_tilde<-diag(c(p_tilde))
     
     #cat("------------------------------------------------------------------------------------\n")
     #cat("initialization: Gaussian - 1 mode \n")
     #cat("-------------------------------------------------------------------------------------\n")
     
     #put this inside a while loop
     lower.bound.diff<-1 #greater than epsilon which will be small in general
     lower.bound.old<-0 #initialize (random value)
     its<-0
     Breakcounter<-0
     
     while ( lower.bound.diff > epsilon)
     {
          its<-1+its
          cat(" GM - Iteration #: ", its,";")
          
          if (its == 100)
          {
               #print("Iterations large")
               Breakcounter<-1
               stop("Iterations large")
          }
          
          #obtain VB posterior for alpha and beta assuming that the optimal distributions are Gaussian mixtures
          
          #if (its==0)
          #{
          #     iter_gm.a<-iterLap(post=log_q_alpha, startVals=rbind(c(alpha_0)))
          #     comps.a<-length(iter_gm.a$weights) #the number of components of the Gaussian mixture - alpha
          
          #     iter_gm.b<-iterLap(post=log_q_beta, startVals=rbind(c(beta_0)))
          #     comps.b<-length(iter_gm.b$weights)#} #the number of components of the Gaussian mixture - beta
          #}else
          #{
          iter_gm.a<-iterLap(post=log_q_alpha, startVals=rbind(c(0,0)))
          comps.a<-length(iter_gm.a$weights) #the number of components of the Gaussian mixture - alpha
          #mean.beta
          iter_gm.b<-iterLap(post=log_q_beta, startVals=rbind(c(0,0)))
          comps.b<-length(iter_gm.b$weights)#} #the number of components of the Gaussian mixture - beta
          #}
          #print(time.iterlap)
          
          #cat("\n iter_gm.a$weights =")
          #print(iter_gm.a$weights)
          #cat("\n iter_gm.a$means =")
          #print(iter_gm.a$means)
          #print(iter_gm.a$sigmas)
          
          #time.moments<-system.time({
          #calculate <a(X%*%beta)>, <b(X%*%beta)>, <a(W%*%alpha)> and <b(W%*%alpha)>
          #for all sites
          #-------------------------------------------------------------------------
          #generates a n by comps.a matrix where each column relates to entries 
          #X%*%iter_gm.a$means[i,] for i=1,...,comps.a. These are the means of the 
          #linear combination of a Gausssian mixture
          mu_beta<-apply(iter_gm.b$means, 1, function(x) X%*%x) 
          
          #generates a n by comps.a matrix where each column relates to entries 
          #sqrt( diag( X%*%iter_gm.a$sigmas[[i]] ) ) for i=1,...,comps.a. 
          #These are the sds of the linear combination of a Gausssian mixture
          sigma_beta<-sapply(iter_gm.b$sigmas, function(x) sqrt(diagElements(X%*%x%*%t(X))) )
          Exp.a_Xb<-matrix(v.gauss.hermite(f=a_x, mu=c(mu_beta), sd=c(sigma_beta), order=50), ncol=comps.b)%*%matrix(iter_gm.b$weights, ncol=1)
          Exp.b_Xb<-matrix(v.gauss.hermite(f=b_x, mu=c(mu_beta), sd=c(sigma_beta), order=50), ncol=comps.b)%*%matrix(iter_gm.b$weights, ncol=1)
          
          mu_alpha<-apply(iter_gm.a$means, 1, function(x) W_vb%*%x) 
          sigma_alpha<-sapply(iter_gm.a$sigmas, function(x) sqrt(diagElements(W_vb%*%x%*%t(W_vb))) )
          Exp.a_Wa<-matrix(v.gauss.hermite(f=a_x, mu=c(mu_alpha), sd=c(sigma_alpha), order=50), ncol=comps.a)%*%matrix(iter_gm.a$weights, ncol=1)
          Exp.b_Wa<-matrix(v.gauss.hermite(f=b_x, mu=c(mu_alpha), sd=c(sigma_alpha), order=50), ncol=comps.a)%*%matrix(iter_gm.a$weights, ncol=1)
          #})
          #print(time.moments)
          
          #Calculate mean and covariance matrices for alpha and beta
          #-------------------------------------------------------------------------
          mean.alpha<- matrix(apply(iter_gm.a$means*iter_gm.a$weights, 2, sum) , ncol=1)
          #cat("\n dim(mean.alpha)=", dim(mean.alpha))
          mean.beta<- matrix(apply(iter_gm.b$means*iter_gm.b$weights, 2, sum), ncol=1)
          #cat("\n mean.beta=",mean.beta)
          
          if (comps.a==1)
          {
               E.alpha.alpha<-iter_gm.a$weights[1]*( iter_gm.a$sigmas[[1]] + crossprod(t(iter_gm.a$means[1,])) )#initialise
          }else
          {
               E.alpha.alpha<-iter_gm.a$weights[1]*( iter_gm.a$sigmas[[1]] + crossprod(t(iter_gm.a$means[1,])) )
               for (k in 2:comps.a)
               {
                    E.alpha.alpha<-E.alpha.alpha + iter_gm.a$weights[k]*( iter_gm.a$sigmas[[k]] + crossprod(t(iter_gm.a$means[k,])) )
               }
          }
          #cat("\n E.alpha.alpha=",E.alpha.alpha)
          
          if (comps.b==1)
          {
               E.beta.beta<-iter_gm.b$weights[1]*( iter_gm.b$sigmas[[1]] + crossprod(t(iter_gm.b$means[1,])) )#initialise
          }else
          {
               E.beta.beta<-iter_gm.b$weights[1]*( iter_gm.b$sigmas[[1]] + crossprod(t(iter_gm.b$means[1,])) )
               for (k in 2:comps.b)
               {
                    E.beta.beta<-E.beta.beta + iter_gm.b$weights[k]*( iter_gm.b$sigmas[[k]] + crossprod(t(iter_gm.b$means[k,])) )
               }
          }
          #cat("\n E.beta.beta=",E.beta.beta)
          
          #Update the occupancy probabilities for each of the sites
          #--------------------------------------------------------------------------
          #if the species is observed at a site the element should be 1
          p_tilde <- matrix(1, ncol=1, nrow=N)
          p<-matrix(1, ncol=1, nrow=n)
          
          #the index of the location of the starting values of Z_tilde
          #Z_tilde <- [(z1,z1,...z1),....,(zn,.......,zn)]^T
          #only keep the ones where the species has not been observed
          starts_z <- matrix(c(1+cumsum(V)-V)[not_obs_vec] , ncol=1)
          ends_z<- starts_z+V_not_obs-1
          starts_2_ends_z<-cbind(starts_z, ends_z)
          
          #Indices<-apply(starts_z,1,function(x,J.=3) x:(x+J.-1))
          #c_Indices<-c(Indices)
          #vectorize this!
          c_Indices<-NULL
          for (c_i in 1:n_not_obs){c_Indices=c(c_Indices, starts_z[c_i]:ends_z[c_i])}
          starts_notz <- matrix(c(1+cumsum(V_not_obs)-V_not_obs) , ncol=1)
          ends_notz<- starts_notz+V_not_obs-1
          starts_2_ends_notz<-cbind(starts_notz, ends_notz)
          
          term1<-Exp.a_Xb[not_obs_vec]
          term2<-Exp.b_Xb[not_obs_vec]
          term3<-sum_Index_vec( Exp.b_Wa[c_Indices], starts_2_ends_notz)
          ti<-term1 - term2 + term3
          
          p[not_obs_vec] <- 1/(1+ exp(-ti) )
          #the group sizes are all different!
          p_tilde[c_Indices]<-rep( p[not_obs_vec], times=V_not_obs)
          P_tilde<-diag(c(p_tilde))
          
          #Calculation of a lower bound for H_alpha(q) and H_beta(q)
          #-----------------------------------------------------------------------------------
          H_alpha<-0 #initialize
          for (k in 1:comps.a)
          {
               c.k<-0 #initialize
               for ( j in 1:comps.a)
               {
                    sigma.j<-iter_gm.a$sigmas[[j]] + iter_gm.a$sigmas[[k]]
                    c.k<-c.k + iter_gm.a$weights[j]*dmvnorm(x=iter_gm.a$means[k,], 
                                                            mean=iter_gm.a$means[j,],
                                                            sigma=sigma.j) 
               }
               H_alpha<- H_alpha + iter_gm.a$weights[k]*log(c.k)
          }
          H_alpha<- -1*H_alpha
          
          H_beta<-0 #initialize
          for (k in 1:comps.b)
          {
               c.k<-0 #initialize
               for ( j in 1:comps.b)
               {
                    sigma.j<-iter_gm.b$sigmas[[j]] + iter_gm.b$sigmas[[k]]
                    c.k<-c.k + iter_gm.b$weights[j]*dmvnorm(x=iter_gm.b$means[k,], 
                                                            mean=iter_gm.b$means[j,],
                                                            sigma=sigma.j) 
               }
               H_alpha<- H_alpha + iter_gm.b$weights[k]*log(c.k)
          }
          H_beta<- -1*H_beta
          
          #Calculation of the lower bound of the marginal log likelihood
          #-----------------------------------------------------------------------------------
          #cat("\n dim(const4)=",dim(const4))
          #cat("\n dim(mean.alpha)=",dim(mean.alpha))
          #cat("\n trace(SigmaInv_alpha_0%*%E.alpha.alpha)=")
          #print(trace(SigmaInv_alpha_0%*%E.alpha.alpha))
          
          Log.alpha<- -0.5*trace(SigmaInv_alpha_0%*%E.alpha.alpha) + const4%*%mean.alpha
          Log.beta<- -0.5*trace(SigmaInv_beta_0%*%E.beta.beta) + const5%*%mean.beta
          #cat("\n Log.alpha=",Log.alpha)
          Logp<- Log.alpha + Log.beta + sum( p*(Exp.a_Xb - Exp.b_Xb) ) + sum( Exp.b_Xb ) + t(Y)%*%P_tilde%*%(Exp.a_Wa - Exp.b_Wa)
          
          E_log_pz <-sum(log(p[not_obs_vec]^p[not_obs_vec]) + log( (1-p[not_obs_vec])^(1-p[not_obs_vec]) ))
          
          lower.bound<- Logp + H_alpha + H_beta - E_log_pz
          #print(lower.bound)
          
          lower.bound.diff <- abs(lower.bound- lower.bound.old)
          lower.bound.old <- lower.bound
          #cat("\n Lower bound difference=", lower.bound.diff[1])
     }
     
     list(iter_gm.a=iter_gm.a, iter_gm.b=iter_gm.b, occup_p=c(p), Breakcounter=Breakcounter)
     
     } 

marg.gm2<-function(x, gm, index)
{
     #The following function takes the output from a function that fits
     #a Gaussian mixture distribution as the VB optimal distribution for
     #alpha and beta and then goes about plotting the marginal distributions
     #of a particular regression slope parameter.
     
     #x = domain of the parameter
     #gm = object from iterLap call that contain the fitted GM
     #index = indicate which parameter one is interested in
     
     #npoints<-10002
     #x<-seq(from=xmin, to=xmax, length.out=npoints)
     npoints=length(x)
     comps<-length(gm$weights)
     y<-NULL
     
     for (i in 1:npoints)
     {
          pdf.ht<-0
          for (j in 1:comps)
          {
               pdf.ht<-pdf.ht + gm$weights[j]*dnorm(x[i], mean=gm$means[j,index],
                                                    sd=sqrt(gm$sigmas[[j]][index,index]))
          }
          y[i]<-pdf.ht
     }
     list(x=x, y=y)
}

peroverlap2=function(x_mcmc,  gm, index, mu, sigma)
{
     #----------------------------------------------------------------------------
     #x_mcmc = the mcmc samples of the parameter
     #range of the parameter is assumed to be [xlo, xhi]
     #xhi=5
     #xlo=-5
     #gsize = the number of points used to perform the kernel density estimation
     #based on bkde from the KernSmooth package
     #npoints = the number of points used to subdivide [xlo,xhi] into
     #the variational posteriors are of the form sum( w[i]*N(mu[i], sigma[i]^2) )
     #We report the percentage overlap between the variational bayes posterior
     #(mixture of Gaussian distributions)
     #and the posterior distribution based on a kernel estimate using the mcmc
     #samples of the parameter
     #gm.a.b = object from iterLap call that contain the fitted GM
     #index = indicate which parameter one is interested in (for the GM)
     #index = 1 ==> the intercept parameter
     #index = 2 ==> the slope parameter
     #mu = the variational mean parameter using a Gaussian aprroximation
     #sd = the variational sd parameter using a Gaussian aprroximation
     #----------------------------------------------------------------------------
     
     require(KernSmooth)
     
     gsize=10000
     npoints=gsize+2
     k=bkde(x_mcmc,bandwidth=dpik(x_mcmc),gridsize=gsize)
     #lines(k$x, k$y, type="l")
     
     #Evaluate the two approximations at k$x
     #mu=vb1$alpha[1,1]
     #sigma=vb1$Sigma_alpha[1,1]^.5
     
     yy=dnorm(k$x, mean=mu, sd=sigma)
     #lines(k$x, yy, col="red")
     
     gm.marg<-marg.gm2(k$x, gm=gm, index=index)
     #lines(k$x, gm.marg$y, col="blue")
     
     po.gm<-1-0.5*sum(abs(gm.marg$y - k$y))*(k$x[2]-k$x[1])
     po.g<-1-0.5*sum(abs(yy - k$y))*(k$x[2]-k$x[1])
     
     #fun=approxfun(gm.marg$x, abs(gm.marg$y-zz),  yleft=0, yright=0, method="linear")
     #result <- integrate(fun, lower=xlo, upper=xhi, subdivisions=10000)
     #return(1-0.5*result$value)
     
     return(c(po.g, po.gm))
}


sink("occ_ss_probit.txt")
cat("
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
         probit(pz[i]) <- beta0 + beta1*X[i]  #covariate = AET_div_PET_s
         
         for (j in 1:nvisits[i]) #loop over the number of visits per site 
         {
              #observation process
              Y[i,j] ~ dbern(py[i,j])
              py[i,j] <- zb[i]*pd[i,j]
              probit(pd[i,j])<- alpha0+alpha1*W[i,j,2]  
         }#end loop over surveys
    }#end loop over sites
    
    alpha0 ~ dnorm(0, 1e-3)
    alpha1 ~ dnorm(0, 1e-3)
    
    beta0 ~ dnorm(0, 1e-3)
    beta1 ~ dnorm(0, 1e-3)
    
    #derived quantity
    occupied <-sum(zb[])
    }##model
    ",fill=TRUE)
sink()
#----------END OF FUNCTIONS TO BE USED -----------------------------------------

#-------------------------------------------------------------------------------
#1. Below we use MLE, VB and MCMC methods to fit a single season occupancy model
#-------------------------------------------------------------------------------

#set the working directory and load some packages
#setwd("add the path to the directory here")
library(MASS);library(msm); library(mvtnorm); require(coda); require(xtable)

#load the R data image. (raw data)
#load("sabap_data_Plos2016.RData")
load("GregDataAnalysis.RData")#load the large R data image. (raw data)

#Have a look at the names of the variables in the RData file
names(stori[[1]])

#'stori' is a list that contains the information for the five species used
#Each list item contains a dataframe with 17 variables. The variables are listed below
#"Cardno" = The observer card number                   
# "Pentad" =  Pentad identifier                   
# "QDGC" = The quarter degree grid cell identifier                      
# "Start_Date" = The date at the start of the observation period              
# "End_Date" = The date at the end of the observation period                    
# "ObserverNo" = Observer number               
# "Total_Spp" = Total number of species observed during the observation period               
# "Intensive_hours" = The number of intensive hours birded        
# "Total_hours" = The total number of hours spent birding              
# "Spp" = What species was observed? An NA indicates that a particular species was not observed.                      
# "Spp_On_Card" = An indicator variable indicating whether or not a particular species was observed              
# "Sequence_At_End_Intensive" = Sequence of birds at the end of your intense birding period
# "Spp_Sequence" = The sequence in which you saw the birds, in SABAP species numbers (i.e. hadedas are spp 84)            
# "Spp_In_Intensive_Hours" = Was the species observed during the intensive period of birding? (yes =1)      
# "Spp_in_hour" = Total number of birds seen per hour            
# "Spp_start_of_hour" = Total number of birds seen at the beginning of each hour        
# "Spp_end_of_hour" = Total number of birds seen at end of each hour

names(climate)

#Climatic variables obtained from 
#Huntley, B., Collingham, Y.C., Green, R.E., Hilton, G.M., Rahbek, C. & Willis, S.G. (2006) 
#Potential impacts of climatic change upon geographical distributions of birds. Ibis, 148, 8-28.
#The climate variables are as follows
# "SABAP_cell_label" = SABAP cell reference
# "centre_longitude" = longitude at the centre of the grid cell
# "centre_latitude" = latitude at the centre of the grid cell
# "GDD0" - annual thermal sum above 0 centigrade
# "GDD5" - annual thermal sum above 5 centigrade
# "MTCO" - mean temperature of the coldest month
# "MTWA"  - mean temperature of the warmest month
# "AET_PET"  - ratio of the annual integrals of actual and potential evapotranspiration
# "Dry_Intensity" - intensity of dry season
# "Wet_Intensity" - intensity of wet season

#Set up cluster to do caluclations!
#Windows version
#require(parallel)
#detectCores()
#require(doParallel)
#cl <- makeCluster(3) 
#registerDoParallel(cl)


sink.T=TRUE
use.species<- c(7, 47, 53, 200) #101:232
bayes.species<-NULL
sink1<-"Analysis2_3_to_5_2012_Bayes_quickrun_paper.txt"
what.year=112

if (sink.T==TRUE){sink(sink1)}
sink(sink1, append=TRUE)

#Do the analysis for each of the species in turn
for (ispecies in use.species)  
#output<- foreach(ispecies = use.species,.combine='rbind',
#                  .packages = c("msm","MASS","mvtnorm","spatstat","iterLap"),
#                  .inorder=T)%dopar%
{
     Timer1=proc.time()
     
     alldat<-stori[[ispecies]]
     
     #make some changes to the raw data
     alldat<-alldat[(!alldat$Sequence_At_End_Intensive>alldat$Total_Spp | is.na(alldat$Sequence_At_End_Intensive>alldat$Total_Spp)),]
     alldat<-alldat[(!alldat$Intensive_hours>alldat$Total_hours | is.na(alldat$Intensive_hours>alldat$Total_hours)),]
     alldat<-alldat[(!alldat$Spp_in_hour>alldat$Total_hours | is.na(alldat$Spp_in_hour>alldat$Total_hours)),]
     alldat<-alldat[(!alldat$Spp_end_of_hour>alldat$Total_Spp | is.na(alldat$Spp_end_of_hour>alldat$Total_Spp)),]
     alldat<-alldat[(!alldat$Spp_in_hour>alldat$Intensive_hours | is.na(alldat$Spp_in_hour>alldat$Intensive_hours)),]
     alldat<-alldat[(!alldat$Spp_Sequence>alldat$Total_Spp | is.na(alldat$Spp_Sequence>alldat$Total_Spp)),]
     alldat<-alldat[(alldat$Total_hours!=0 | is.na(alldat$Total_hours!=0)),]
     alldat<-alldat[(alldat$Intensive_hours>=2 | is.na(alldat$Intensive_hours>=2)),]
     alldat<-alldat[!(alldat$Spp_In_Intensive_Hours==0 & !is.na(alldat$Spp_in_hour)),] #this one does not run
     
     #create the Climate variables
     #------------------------------------------------------------------
     alldat$GDD0<-0
     alldat$GDD5<-0
     alldat$MTCO<-0
     alldat$MTWA<-0
     alldat$AET_PET<-0
     alldat$Wet_Intensity<-0
     alldat$Dry_Intensity<-0
     
     #create the climate vaiables
     unique.QDGC<-as.character(unique(alldat$QDGC))
     QDGC<-as.character(alldat$QDGC)
     l.unique.QDGC<-length(unique.QDGC)
     climate_SABAP_cell_label<-as.character(climate$SABAP_cell_label)
     
     for (iQDGC in 1:l.unique.QDGC)
     {
          index.in.climate_SABAP_cell_label<-which( climate_SABAP_cell_label==unique.QDGC[iQDGC]) 
          index.in.QDGC<-which( QDGC==unique.QDGC[iQDGC])
          
          alldat$GDD0[index.in.QDGC]<-climate$GDD0[index.in.climate_SABAP_cell_label]
          alldat$GDD5[index.in.QDGC]<-climate$GDD5[index.in.climate_SABAP_cell_label]
          alldat$MTCO[index.in.QDGC]<-climate$MTCO[index.in.climate_SABAP_cell_label]
          alldat$MTWA[index.in.QDGC]<-climate$MTWA[index.in.climate_SABAP_cell_label]
          alldat$AET_PET[index.in.QDGC]<-climate$AET_PET[index.in.climate_SABAP_cell_label]
          alldat$Wet_Intensity[index.in.QDGC]<-climate$Wet_Intensity[index.in.climate_SABAP_cell_label]
          alldat$Dry_Intensity[index.in.QDGC]<-climate$Dry_Intensity[index.in.climate_SABAP_cell_label]
     }
     #------------------------------------------------------------------
     
     Years<-strptime(as.character(alldat$Start_Date),format="%Y-%m-%d")$year
     
     #------------------------------------------------------------------
     #Select the correct subset of the data
     #we want to select a particular years data
     #and then restrict the analysis to those sites that has at least 3 visits and at most 5
     #this is done later however
     
     #select all of the rows with Year = 2012 (112)
     Rows<-which(Years==what.year)
     #cat("\n The number of rows in the subset data file = ", length(Rows))
     
     alldat<-alldat[Rows,] #select only the data for a particular year!!!
     #------------------------------------------------------------------
     #remove factors that is not in the subset
     alldat$Pentad<-factor(alldat$Pentad)
     alldat$QDGC<-factor(alldat$QDGC)
     
     #make some new variables
     alldat[,25:27]<-1 #dummy variables
     
     #this variable does change and is different for each particular visit
     alldat[, 28]<- (alldat[,7]-mean(alldat[,7]))/sd(alldat[,7]) #standardized total number of species
     colnames(alldat)[28]<-"nspp"
     
     #the huntley variables standardized across all of the data
     #these should not be used
     #this variable does not change with visits and is site specific
     alldat[,29:34]<-1 #dummy variables
     alldat[, 35]<- (alldat[,18]-mean(alldat[,18]))/sd(alldat[,18]) #standardized GDD0
     alldat[, 36]<- (alldat[,19]-mean(alldat[,19]))/sd(alldat[,19]) #standardized GDD5
     alldat[, 37]<- (alldat[,20]-mean(alldat[,20]))/sd(alldat[,20]) #standardized MTCO
     alldat[, 38]<- (alldat[,21]-mean(alldat[,21]))/sd(alldat[,21]) #standardized MTWA
     alldat[, 39]<- (alldat[,22]-mean(alldat[,22]))/sd(alldat[,22]) #standardized AET_div_PET
     alldat[, 40]<- (alldat[,23]-mean(alldat[,23]))/sd(alldat[,23]) #standardized Wet_Intensity
     alldat[, 41]<- (alldat[,24]-mean(alldat[,24]))/sd(alldat[,24]) #standardized Dry_Intensity
     alldat[,42:44]<-1 #dummy
     colnames(alldat)[35:41]<-c("GDD0_s","GDD5_s","MTCO_s",  "MTWA_s",  "AET_div_PET_s", "Wet_Intensity_s","Dry_Intensity_s")
     
     X<-alldat

     #detection covariates
     out.intensive <- unstack(X[ , c(8, 2)]) #number of intensive hours birding
     out.total <- unstack(X[ , c(9, 2)] ) # total number of hours birding
     out.nspp <- unstack(X[ , c(7, 2)])
     
     #occupancy covariates
     #The huntley variables 
     out.GDD0<-unstack(X[ , c(18, 2)])
     out.GDD5<-unstack(X[ , c(19, 2)])
     out.MTCO<-unstack(X[ , c(20, 2)])
     out.MTWA<-unstack(X[ , c(21, 2)])
     out.AET_div_PET<-unstack(X[ , c(22, 2)])
     out.Wet_Intensity<-unstack(X[ , c(23, 2)])
     out.Dry_Intensity<-unstack(X[ , c(23, 2)])

     #------------------------------------------------------------------
     #All of the sites and the detection-nondetection data
     out.y <- unstack(X[ , c(11, 2)]) #Spp_On_Card = alldat[ Rows, c(11, 2)]
     
     #No limit has been made of the number of surveys 
     len.y<-NULL
     max.y<-NULL
     ndetections<-NULL
     i=1
     for (i in 1:length(out.y))
     {  
          if( is.null( out.y[[i]] ) == FALSE)
          {
               len.y[i]<- length(out.y[[i]]) 
               max.y[i]<- max(out.y[[i]], na.rm=T) 
               ndetections[i]<- sum(out.y[[i]]>0, na.rm=T)
               i<-i+1
          }
     }
     maxSurveys<-max(len.y, na.rm=T)

     #create dummy data matrix
     DummyMatInit<-matrix(NA, nrow= length(unique(X$Pentad )), ncol = maxSurveys)
     y <- DummyMatInit # absence presence data
     
     j <- 1
     for(i in 1:length(out.y) )
     {
          if( is.null( out.y[[i]] ) == FALSE)
          {
               Len.y.i<-length( out.y[[i]])
               y[j, 1:Len.y.i] <- out.y[[i]]
               j <- j + 1
          }  
     }#now y is of dimension nsites by nvisits
     
     #Sub set the data and identify which rows we want to include in the analysis  
     Nvisits<-apply(y[,],1,function(x){length(na.omit(x))})
     
     table(Nvisits)

     #choose all in 2012
     nvisits.fit<-5#the number of site visits allowed in the analysis
     choose.rows<-which(Nvisits>=3 & Nvisits<=5)#which(Nvisits==3)
     
     #only use the chosen set of sites
     Y<-y[choose.rows, 1:nvisits.fit]
     DummyMatInit<-matrix(NA, nrow= length(choose.rows), ncol = maxSurveys)
     
     #detection covariates
     intensive <- DummyMatInit 
     nspp <- DummyMatInit 
     
     j <- 1
     for(i in choose.rows )
     {
          if( is.null( out.y[[i]] ) == FALSE)
          {
               #cat("\n i,j = ",i,j)
               Len.y.i<-length( out.y[[i]] )
               intensive[j, 1:Len.y.i] <- out.intensive[[i]]
               nspp[j, 1:Len.y.i] <- out.nspp[[i]]
               j <- j + 1
          }  
     }#now these variables are of dimension length(choose.rows) by maxSurveys
     
     #Standardise the detection covariates
     #we do this by standardizing across each of the visits
     for (ivisits in 1:nvisits.fit)
     {
          intensive[,ivisits] <- scale(intensive[,ivisits] ,T,T)
          nspp[,ivisits] <- scale(nspp[,ivisits] ,T,T)
     }

     #------------------------------------------------------------------
     # and now put the matrices into the correct format
     #------------------------------------------------------------------
     # detection covariates
     obs.covs <- list(intensive=intensive[,1:nvisits.fit], nspp=nspp[,1:nvisits.fit])
     #------------------------------------------------------------------
     # Clean up  the site specific variables
     #------------------------------------------------------------------
     #occupancy covariates are standardized
     site.dat<-unique(X[ , c(2, 35:41)])[choose.rows,]  #change later (include all of the X's Pentad is not a regressor!)
     site.dat$Pentad<-factor(site.dat$Pentad)
     for (ivars in 2:8){site.dat[,ivars]<-scale(site.dat[,ivars],T,T)}
     #OMdat <- unmarkedFrameOccu(y=Y, siteCovs = site.dat, obsCovs=obs.covs)
     #------------------------------------------------------------------
     
     Nvisits<-apply(y[choose.rows,],1,function(x){length(na.omit(x))})

     #------------------------------------------------------------------
     #The MLE fit - DONT WANT TO DO THE MLE!
     #------------------------------------------------------------------
     #f8<- occu( ~ nspp  ~AET_div_PET_s , data=OMdat,engine="C")
     #f8.summary<-try(summary(f8),T)
     #pvalues<-sum(c(sum(f8.summary[["state"]]["P(>|z|)"] <.05, na.rm=T), sum(f8.summary[["det"]]["P(>|z|)"]<0.05, na.rm=T) ))
     
     #if (pvalues==4) #if the MLE can be performed then continue
     #{
          bayes.species<-c(bayes.species,ispecies) #append
          cat("\n ....................................\n")
          #cat("\n ", as.character(sp[ispecies]), "- Speciesid = ", ispecies)
          cat("\n Speciesid = ", ispecies)
          cat("\n ....................................\n")
          
          cat("\n Number of sites = ", length(len.y[len.y!=0]))
          cat("\n Max surveys = ", maxSurveys)
          print( table(len.y[len.y!=0])[1:10] )
          cat("\n ....................................\n")
          
          #------------------------------------------------------------------
          #the VB fit
          #------------------------------------------------------------------

          ymat<-Y
          Xmat<-site.dat      
          Wmat <- list(intensive=intensive[,1:nvisits.fit], nspp=nspp[,1:nvisits.fit])
          design_mats<-vb_Designs(W=Wmat, X=Xmat, y=ymat)

          #Set the prior distributions used
          alpha_0 <- matrix(0, ncol=1, nrow=2) 
          beta_0 <- matrix(0, ncol=1, nrow=2)
          Sigma_beta_0 <- diag(2)*sqrt(1000)
          Sigma_alpha_0 <- diag(2)*sqrt(1000)
          
          #set the formula used
          form1<- ymat~ AET_div_PET_s  ~ nspp 
          
          
          Y.mat<-Y 
          if (  sum(apply(Y.mat, 1, function(x) max(x, na.rm = T)))/ nrow(Y.mat) > 0.45)
          {
               #VB fit - assuming that the posterior distributions can be approximated 
               #using a mixture of Gaussian distributions
               gm1<-vb_model2_gm(formula=form1, design_mats=design_mats, 
                                 alpha_0=alpha_0, beta_0=beta_0, 
                                 Sigma_alpha_0=Sigma_alpha_0, Sigma_beta_0=Sigma_beta_0, 
                                 epsilon=1e-3)
               
               #my code does not allow one to fit the simple model! i.e. 
               #an intercept only model cannot be fitted
               #The call to the Laplace algorithm VB code
               vb_a1<-vb_model2_prob_la(formula=form1, design_mats = design_mats,alpha_0 = alpha_0,beta_0 = beta_0,
                                      Sigma_alpha_0 = Sigma_alpha_0,Sigma_beta_0 = Sigma_beta_0, 
                                      epsilon=1e-5, Print=FALSE)
     
               
               #Store the results
               vb_g_outs<-rbind(cbind(c(vb_a1$beta), sqrt(diag(vb_a1$Sigma_beta))), cbind( c(vb_a1$alpha), sqrt(diag(vb_a1$Sigma_alpha)) ))
               vb_g_outs<-cbind(vb_g_outs,vb_g_outs[,1]/vb_g_outs[,2])
               colnames(vb_g_outs)<-c("Est","SD","t")
               
               cat("\n VB results \n")
               print(vb_g_outs)
               cat("\n ....................................\n")
               
               cat("\n VB GM results \n")
               print(gm1)
               cat("\n ....................................\n")
               
               #cat("\n MLE results \n")
               #print(f8)
               #cat("\n")       
               #mle.fit<-f8
     
               #------------------------------------------------------------------
               # The Bayesian fit
               #------------------------------------------------------------------
     
               #Y.mat<-Y
               X.use<-Xmat[,6]
               
               #set the number of draws to be used in simulation
               nsims = 1000
               #set the starting values as random values centered around the VB estimates
               a=runif(1,vb_g_outs[1,1]-vb_g_outs[1,2], vb_g_outs[1,1]+vb_g_outs[1,2]) 
               b=runif(1,vb_g_outs[2,1]-vb_g_outs[2,2], vb_g_outs[2,1]+vb_g_outs[2,2])
               e=runif(1,vb_g_outs[3,1]-vb_g_outs[3,2], vb_g_outs[3,1]+vb_g_outs[3,2]) 
               d=runif(1,vb_g_outs[4,1]-vb_g_outs[4,2], vb_g_outs[4,1]+vb_g_outs[4,2])
               
               n.use<-length(choose.rows)
               W = array(dim=c(n.use,5,2)) #number of observations ; number of visits; number of covariates + 1 (for the intercept)
               W[,,1] = 1 #intercept
               
               for (i in 1:5)
               {
                    W[,i,2]<-nspp[,i]
               }
               
               #identify which of the elements of W are missing using the observed data matrix
               W1=matrix(NA,nrow=n.use, ncol=5) #ncol=#visits
               W1[1:n.use,]<- W[1:n.use,,2]
               
               #here we ensure that the W1 matrix has 'NA' values whenever Y has 'NA' values
               for (i in 1:n.use)#sites
               {
                    for (j in 1:5)#site visits
                    {
                         #cat("\n", i, " ", j, " " , is.na(Y[i,j]), "\n")
                         W1[i,j]<-ifelse( is.na(Y.mat[i,j])==TRUE, NA, W1[i,j]) #the 1st detection covariates for the different site visits
                    }
               }
               W[,,2] = W1
               
               nvisits<-apply(Y.mat,1,function(x){length(na.omit(x))}) #the number of visits to each site
               
               
               #sample from posterior distribuion using DA algorithm
               X<-cbind(1, X.use)
     #check           
               thin.post=DA(X=X, W=W, Y=Y.mat, J=5, 
                            ndraws=nsims, alpha.start=c(e,d), beta.start=c(a,b),thin=1)
               
               #the output obtain from the DA function
               alpha0<-thin.post[,3]
               alpha1<-thin.post[,4]
               beta0<-thin.post[,1]
               beta1<-thin.post[,2]
     
               #the posterior samples for the occupancy prob at each site
               occupied = thin.post[,5]*dim(X)[1]#thin.post[, -c(1:5)]
               
               #posterior samples of number of occupied sites
               zb<-thin.post[,5]*dim(X)[1]
               
               #the MCMC posterior distribution of the parameters as well as the Vb posterior distribution
               #-------------------------------------------------------------------------------------------
               pdf(paste("SABAP_2013_Species_",ispecies,"_year_",what.year,".pdf",sep=""))
               #par(mfrow=c(2,2), mar=c(4,4,4,1))
               
               par(mfrow=c(2,2), mar=c(4,4,1,1))
               p1=density(alpha0)
               xx=seq(from=(min(p1$x)-.25), to=(max(p1$x)+.25), length.out=1000)
               d1=dnorm(xx, mean=vb_a1$alpha[1], sd=sqrt(vb_a1$Sigma_alpha[1,1]))
               hip1=max(p1$y)
               hip2=max(d1)         
               h<-hist(alpha0, plot=F, breaks=50)
               plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab=expression(alpha[0]), main="", ylab="Density")
               lines(xx, d1, col="red")
               legend("topleft", c("G","GM"), lwd=c(1,2), col=c("red","blue"),bty="n")
               p1gm=marg.gm2(x=seq(from=min(alpha0), to=max(alpha0), length.out=10002), gm=gm1$iter_gm.a, index=1) 
               lines(p1gm$x, p1gm$y, col="blue", lwd=2)
               
               #plots for alpha1
               p1=density(alpha1)
               xx=seq(from=min(p1$x), to=max(p1$x), length.out=1000)
               d1=dnorm(xx, mean=vb_a1$alpha[2], sd=sqrt(vb_a1$Sigma_alpha[2,2]))
               hip1=max(p1$y)
               hip2=max(d1)         
               h<-hist(alpha1, plot=F, breaks=50)
               plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab="nspp", main="", ylab="Density")
               lines(xx, d1, col="red")
               legend("topleft", c("G","GM"), lwd=c(1,2), col=c("red","blue"),bty="n")
               p1gm=marg.gm2(x=seq(from=min(alpha1), to=max(alpha1), length.out=10002), gm=gm1$iter_gm.a, index=2) 
               lines(p1gm$x, p1gm$y, col="blue", lwd=2) 
               
               #plots for beta0
               p1=density(beta0)
               xx=seq(from=min(p1$x), to=max(p1$x), length.out=1000)
               d1=dnorm(xx, mean=vb_a1$beta[1], sd=sqrt(vb_a1$Sigma_beta[1,1]))
               hip1=max(p1$y)
               hip2=max(d1)         
               h<-hist(beta0, plot=F, breaks=50)
               plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab=expression(beta[0]), main="", ylab="Density")
               lines(xx, d1, col="red")
               legend("topleft", c("G","GM"), lwd=c(1,2), col=c("red","blue"),bty="n")
               p1gm=marg.gm2(x=seq(from=min(beta0), to=max(beta0), length.out=10002), gm=gm1$iter_gm.b, index=1) 
               lines(p1gm$x, p1gm$y, col="blue", lwd=2) 
               
               #plots for beta1
               p1=density(beta1)
               xx=seq(from=min(p1$x), to=max(p1$x), length.out=1000)
               d1=dnorm(xx, mean=vb_a1$beta[2], sd=sqrt(vb_a1$Sigma_beta[2,2]))
               hip1=max(p1$y)
               hip2=max(d1)         
               h<-hist(beta1, plot=F, breaks=50)
               plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab="AETdivPETs", main="", ylab="Density")
               lines(xx, d1, col="red")
               legend("topleft", c("G", "GM"), lwd=c(1,2), col=c("red","blue"),bty="n")
               p1gm=marg.gm2(x=seq(from=min(beta1), to=max(beta1), length.out=10002), gm=gm1$iter_gm.b, index=2) 
               lines(p1gm$x, p1gm$y, col="blue", lwd=2)
               
               dev.off()
               
               #file.save<-paste("Analysis_3_to_5_2012_Bayes_Species_",ispecies,"_year_",what.year,".RData",sep="")
               #save.image(file.save) 
               
               Timer2=proc.time()
               cat("\n Timer = ", Timer2-Timer1, "\n")
               
               Timer2-Timer1
          }else
          {
               cat("\n observed detection prob < 0.45 \n")
          }
     #}
     #-----------------------------------------------------------------------
}
#stopCluster(cl)

cat("\n DONE RUNNING THE ANALYSIS! \n")

if (sink.T==T){sink()}
sink()
sink()
sink()


