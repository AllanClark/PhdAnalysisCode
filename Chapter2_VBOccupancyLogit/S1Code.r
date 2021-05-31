#----------FUNCTIONS TO BE USED -----------------------------------------------
peroverlap=function(x_mcmc, xlo, xhi, mu, sigma)
{
     #----------------------------------------------------------------------------
     #x_mcmc = the mcmc samples of the parameter
     #range of the parameter is assumed to be [xlo, xhi]
     #xhi=5
     #xlo=-5
     #gsize = the number of points used to perform the kernel density estimation
     #based on bkde from the KernSmooth package
     #npoints = the number of points used to subdivide [xlo,xhi] into
     #the variational posteriors are of the form N(mu, sigma^2)
     #We report the percentage overlap between the variational bayes posterior
     #and the posterior distribution based on a kernel estimate using the mcmc
     #samples of the parameter
     #----------------------------------------------------------------------------
     
     xlo=min(xlo, min(x_mcmc))
     xhi=max(xhi, max(x_mcmc))
     
     gsize=10000
     npoints=gsize+2
     k=bkde(x_mcmc,bandwidth=dpik(x_mcmc),gridsize=gsize)
     
     #number of point to add on the lower end
     nl=(npoints-gsize)/2
     nu=nl
     
     xl=seq(from=xlo, to=min(k$x), length.out=nl)
     yl=c(rep(0,nl-1), k$y[which(k$x==min(k$x))])
     
     xh=seq(from=max(k$x), xhi, length.out=nu)
     yh=c(k$y[which(k$x==max(k$x))],rep(0,nu-1))
     
     zz=c(yl,k$y,yh)
     
     #allowable range of parameter
     xx=c(xl, seq(from=min(k$x), to=max(k$x), length.out=gsize), xh)
     
     #y values of normal distribution
     yy=dnorm(xx, mean=mu, sd=sigma)
     
     fun=approxfun(xx, abs(yy-zz),  yleft=0, yright=0, method="linear")
     result <- integrate(fun, lower=xlo, upper=xhi, subdivisions=10000)
     
     return(1-0.5*result$value)
}

#VB: Laplace approximation
vb_model2_la<-function(formula, design_mats, alpha_0, beta_0, Sigma_alpha_0, 
                       Sigma_beta_0, LargeSample=TRUE, epsilon=1e-5, Print=FALSE)
{
     #Date: 20 October 2014
     #Allan E Clark
     #Multi-visit site occupancy model
     #Estimation is undertaken using variational bayes and Laplace approximation
     #----------------------------------------------------------------------------------
     
     #----------------------------------------------------------------------------------
     #Arguments
     #----------------------------------------------------------------------------------
     #formula<- y~ occupancy covariates  ~ site detection covariates
     #y <- n by J matrix of presence absence data
     #n <- the number of locations
     #J <- the number of visits to the sites
     #Assumed that each site is visited J times
     
     #X <- a dataframe that contains the covariates used to calculate site occupancy probabilities
     #W <- a named list that contains the site covariates used to calculate the site detection probabilities
     #W has the same form as an unmarkedFrameOccu object in the unmarked package
     
     #alpha_0 <- starting value of detection covariate coefficients
     #beta_0 <- starting value of occurence covariate coefficients
     #Sigma_alpha_0, Sigma_beta_0 - the variance covariance matrix of alpha_0 and beta_0
     #LargeSample<-TRUE - indicates that the number of sites is 'large' and that an
     #approximation to B(mu, sigma^2) is used instead of integrations
     #-----------------------------------------------------------------------------------
     
     #load the required functions
     bx<-function(x){log(1+exp(x))}
     
     b1<-function(x)
     {
          #the first deriv of bx
          1/(1+exp(-x))
     }
     
     b2<-function(x)
     {
          #the second deriv of bx
          exp(-x)/( (1+exp(-x))^2 )
     }
     
     Ebx <- function(x, mu, sigma2)
     {
          #bx(mu+sigma*x)*dnorm(x)
          argx<-mu+sqrt(sigma2)*x
          log(1+exp(argx))*dnorm(x)
     }
     
     B<-function(mu,sigma2)
     {
          #This does Brute force integrations
          #Note that the range is limited to -200 to 200
          integrate( Ebx, lower=-200, upper =200, mu=mu, sigma2=sigma2)$value
     }
     
     B2<-function(mu, sigma2)
     {
          #Approximation to B(mu,sigma2) for large sample sizes
          bx(mu)+b2(mu)*sigma2
     }
     
     B0.approx <- function(mu,sigma2)
     {
          #John Ormerod code
          sigma <- sqrt(sigma2)
          vB0 <- mu*pnorm(mu/sigma) + sigma2*dnorm(mu,sd<-sigma)
          vB0 <- vB0 + (0.6931472*exp(0.3298137*sigma2-0.8121745*mu))*pnorm( mu/sigma - 0.8121745*sigma)
          vB0 <- vB0 + (0.6931472*exp(0.3298137*sigma2+0.8121745*mu))*pnorm(-mu/sigma - 0.8121745*sigma)
          return(vB0)
     }
     
     B1.approx <- function(mu,sigma2)
     {
          #John Ormerod code
          sigma <- sqrt(sigma2)
          muOnSigma <- mu/sigma
          a1sigma <- 0.8121745*sigma
          a1mu <- 0.8121745*mu
          halfa1sqsigma2 <- 0.3298137*sigma2
          vB1 <- pnorm(muOnSigma)
          vB1 <- vB1 + 0.5629565*(exp(halfa1sqsigma2+a1mu)*pnorm(-muOnSigma - a1sigma) 
                                  - exp(halfa1sqsigma2-a1mu)*pnorm( muOnSigma - a1sigma))
          return(vB1)
     }
     
     B2.approx <- function(mu,sigma2)
     {
          #John Ormerod code
          sigma <- sqrt(sigma2)
          muOnSigma <- mu/sigma
          a1sigma <- 0.8121745*sigma
          a1mu <- 0.8121745*mu
          halfa1sqsigma2 <- 0.3298137*sigma2
          vB2 <- -0.1259130*dnorm(mu,sd<-sigma)
          vB2 <- vB2 + 0.4572189*((exp(halfa1sqsigma2-a1mu))*pnorm(muOnSigma - a1sigma) 
                                  + (exp(halfa1sqsigma2+a1mu))*pnorm(-muOnSigma - a1sigma))
          return(vB2)
     }
     
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
     
     logp<-function(par, W, X, Y, P_tilde, p_tilde, p, alpha_0, SigmaInv_alpha_0, beta_0, SigmaInv_beta_0, Exp=0)
     {
          #log(q(alpha, beta)) when 'Exp<-0'
          #the log prior distribution componenti is included in the calculations here
          
          ncol_W <- NCOL(W)
          
          alpha<- matrix(par[1:ncol_W], ncol=1)
          beta<- matrix(par[-c(1:ncol_W)] , ncol=1)
          
          alpha_x<-W%*%alpha
          alpha_diff <- alpha - alpha_0
          beta_x<-X%*%beta
          beta_diff <- beta - beta_0
          
          t1 <- crossprod(Y,P_tilde)%*%alpha_x - crossprod(p_tilde,bx(alpha_x) ) 
          + crossprod(p,beta_x) - sum( bx(beta_x) )
          t2 <- 0.5*(1-Exp)*( crossprod(alpha_diff,SigmaInv_alpha_0)%*%alpha_diff 
                              + crossprod(beta_diff,SigmaInv_beta_0)%*%beta_diff  )
          
          Logp <- t1[1]  - t2[1]
          return(Logp)
     }
     
     vb_Designs_check<-function(formula, Names)
     {
          #perform some checks on the design matrices
          #check that the names in the formula call are in the dataframes provided
          
          detVars<-all.vars(formula)
          
          if ( (sum(detVars[-1]%in% Names)==length(detVars[-1]))!= 1)
          {
               stop(print("\n \n CHECK YOUR FORMULA CALL. \n MISMATCH BETWEEN 
                          CALL AND DESIGN MATRICES. \n i.e. You included objects 
                          in the call: '~occupancy variables ~ detection variables' 
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
     
     #----------------------------------------------------------------------------------
     
     #Declare certain matrices and constants
     #---------------------------------------
     #create a matrix where the elements of y are stored one vector below the other
     #Y<-matrix(t(y), ncol=1)
     #X<-as.matrix(X) #X stored as a matrix
     
     #formula<- y~X1+X2~W1+W2+W3
     #design_mats<-vb_Designs(W=W.eg.l1, X=X.eg, y=y)
     
     #req_design_mats<-vb_ReqDesigns(formula=form1, design_mats)
     req_design_mats<-vb_ReqDesigns(formula, design_mats)
     W_vb<-req_design_mats$W
     X<-req_design_mats$X
     Y<-matrix(design_mats$Y$V1)
     
     n<-NROW(X)
     N<-length(Y)
     pn<-NCOL(X) #number of columns of X
     qn<-NCOL(W_vb) #number of columns of W
     
     const1 <- (pn+qn)*(log(2*pi)+1) 
     #+1 added since we require it for the calculation of the entropy of a mult norm
     
     const2<-determinant(Sigma_alpha_0, logarithm=T)$modulus[1] + 
          determinant(Sigma_beta_0, logarithm=T)$modulus[1]
     const3<- 0.5*(const1+const2)
     
     SigmaInv_alpha_0 <- chol2inv(chol(Sigma_alpha_0))
     SigmaInv_beta_0 <- chol2inv(chol(Sigma_beta_0))
     SigmaInv_times_alpha_0 <- SigmaInv_alpha_0%*%alpha_0
     SigmaInv_times_beta_0 <- SigmaInv_beta_0%*%beta_0
     #---------------------------------------
     
     #a matrix containing how many times each of the sites are visited
     #J<-dim(W)[2]
     V<-design_mats$nvisits #nvisits
     
     #create the W matrix
     #Note that here were are assuming that V(i) are all equal
     #W_vb <- cbind(1,c(t(W[,,qn])))
     
     #the prior mean and covariance matrix of alpha and beta
     #------------------------------------------------------
     #the starting value for the alpha and beta vector (and covariance matrices)
     #--------------------------------------------------------------------------
     
     alpha <- alpha_0*rnorm(1)*0  
     beta <- beta_0 *rnorm(1)*0 
     alpha_0<-alpha
     beta_0<-beta
     Sigma_alpha <- Sigma_alpha_0
     Sigma_beta <- Sigma_beta_0
     Sigma_alpha_inv <-solve(Sigma_alpha_0)
     Sigma_beta_inv <- solve(Sigma_beta_0)
     
     oldparams<- c(alpha, Sigma_alpha_inv, beta, Sigma_beta_inv)
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
     
     t1<-apply(X[not_obs_vec,],1, function(x,beta.=beta){crossprod(x,beta.)})
     
     starts_notz <- matrix(c(1+cumsum(V_not_obs)-V_not_obs) , ncol=1)
     ends_notz<- starts_notz+V_not_obs-1
     starts_2_ends_notz<-cbind(starts_notz, ends_notz)
     t2<-sum_Index_vec(bx(W_vb[c_Indices,]%*%alpha) , starts_2_ends_notz)
     p[not_obs_vec] <- 1/(1+ exp(-t1 + t2) )
     
     #the group sizes are all different!
     p_tilde[c_Indices]<-rep( p[not_obs_vec], times=V_not_obs)
     P_tilde<-diag(c(p_tilde))
     
     #Marginal likelihood approximation <- E(Log_p(y,z,alpha,beta))- E(q(alpha, beta, z))
     #-----------------------------------------------------------------------------------
     Logp <- logp(c(alpha,beta), W=W_vb, X=X, Y=Y, P_tilde=P_tilde, p_tilde=p_tilde, p=p, alpha_0=alpha_0, SigmaInv_alpha_0=SigmaInv_alpha_0, beta_0=beta_0, SigmaInv_beta_0=SigmaInv_beta_0,Exp=1)
     #Logp <- Logp -.5*( const1 + log(det(Sigma_alpha_0)) + log(det(Sigma_beta_0)) )  #E(Log_p(y,z,alpha,beta))

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
     
     #  cat("------------------------------------------------------------------------------------","\n")
     #  cat("STARTING VALUES","\n")
     #  cat("Log_ml <- ", Log_ml,", alpha <- ", round(alpha,digits<-3) , ", beta <- ", round(beta,digits<-3), "\n")
     #  cat("------------------------------------------------------------------------------------","\n")
     
     old_Log_ml<-Log_ml
     
     diff_Log_ml<-1
     its<-0
     
     #epsilon <- 1e-10
     Breakcounter<-0
     
     while (diff_Log_ml > epsilon)
     {
          its<-1+its
          #cat("Iteration #: ", its,"\n")
          #cat("------------------------ \n")
          
          if (its == 2000)
          {
               #print("Iterations large")
               Breakcounter<-1
               #list(alpha=alpha, beta=beta, Sigma_alpha=Sigma_alpha, Sigma_beta=Sigma_beta, occup_p=c(p), Log_mla=Log_ml, Breakcounter=Breakcounter)
               stop("Iterations large")
          }
          
          #Perform the Newton Rhaphson algortihm to calculate alpha and beta at iteration t
          #---------------------------------------------------------------------------------
          #crit <- .5 #used to assess convergence of the parameters - Rule 1
          crit<-0 #used for Rule 3
          
          #These are used in the Newton Rhaphson algortihm but don't change within the 'while' loop
          g_alpha_1<- crossprod(W_vb, P_tilde%*%Y)
          g_beta_1<- crossprod(X, p)
          
          iinside=0
          
          while (crit < nparams)
          {
               alpha_x<-W_vb%*%alpha
               g_alpha <- g_alpha_1 - crossprod(W_vb, p_tilde*b1(alpha_x)) - SigmaInv_alpha_0%*%alpha + SigmaInv_times_alpha_0
               Sigma_alpha_inv <-   XWX.fun2(W_vb, p_tilde*b2(alpha_x)) + SigmaInv_alpha_0 #Note Sigma_alpha is not stored!
               alpha <- alpha + solve(Sigma_alpha_inv,g_alpha)
               
               beta_x<-X%*%beta
               g_beta <- g_beta_1 -crossprod(X, b1(beta_x)) - SigmaInv_beta_0%*%beta + SigmaInv_times_beta_0
               Sigma_beta_inv <- XWX.fun2(X, b2(beta_x)) + SigmaInv_beta_0 #Note Sigma_beta is not stored!
               beta <- beta + solve(Sigma_beta_inv, g_beta)
               
               #Stopping method 3
               #-----------------
               newparams<-c(alpha, Sigma_alpha_inv, beta, Sigma_beta_inv)
               crit<-sum(abs(newparams-oldparams) <= epsilon) #STRONG CONDITION!
               oldparams<-newparams
               
               iinside<-iinside+1
               #cat("iinside ",iinside,"\n")
               
               if (iinside == 2000)
               {
                    #print("Newton Rhap iterations large")
                    Breakcounter<-1
                    #list(alpha=alpha, beta=beta, Sigma_alpha=Sigma_alpha, Sigma_beta=Sigma_beta, occup_p=c(p), Log_mla=Log_ml, Breakcounter=Breakcounter)
                    stop("Newton Rhap iterations large")
               }
          }
          #iinside=0
          #Now we calculate the covariance matrices
          Sigma_alpha<-solve(Sigma_alpha_inv)
          Sigma_beta<-solve(Sigma_beta_inv)
          
          #Calculate the occupancy probabilities at the different locations
          #-----------------------------------------------------------------
          
          t1<-apply(X[not_obs_vec,],1, function(x,beta.=beta){crossprod(x,beta.)}) #correct (not_obs_vec by 1 vector)
          ag1<-W_vb[c_Indices,]%*%alpha #'mu' for B2
          ag2<-X_row.Y_col(W_vb[c_Indices,],Sigma_alpha%*%t(W_vb[c_Indices,])) #'sigma2' for B2
          
          if (LargeSample== TRUE)
          {
               #the group sizes are all different! so this does not work nomore...
               #t2<-apply( matrix(B2(ag1, ag2), nrow=n_not_obs , byrow=T), 1, sum)
               t2<-sum_Index_vec(B2(ag1, ag2) , starts_2_ends_notz)
          }else
          {
               t2<- sum_Index_vec( matrix(apply(cbind(ag1,ag2),1, function(x){B(x[1],x[2])}), ncol=1), starts_2_ends_notz)
          }
          p[not_obs_vec] <- 1/(1+ exp(-t1 + t2) )
          p_tilde[c_Indices]<-rep( p[not_obs_vec], times=V_not_obs)
          P_tilde<-diag(c(p_tilde))
          
          #Marginal likelihood approximation <- E(Log_p(y,z,alpha,beta))- E(q(alpha, beta, z))
          #-----------------------------------------------------------------------------------
          Logp <- logp(c(alpha,beta), W=W_vb, X=X, Y=Y, P_tilde=P_tilde, p_tilde=p_tilde, p=p, alpha_0=alpha_0, SigmaInv_alpha_0=SigmaInv_alpha_0, beta_0=beta_0, SigmaInv_beta_0=SigmaInv_beta_0,Exp=1)
          #Logp <- Logp -const3  #E(Log_p(y,z,alpha,beta))
          #<ln( p(alpha, beta) )> when one is integrating wrt the optimal distribution!
          alpha_diff<-alpha - alpha_0
          beta_diff<-beta - beta_0
          f1<- trace(SigmaInv_alpha_0%*%Sigma_alpha) + crossprod(alpha_diff , SigmaInv_alpha_0)%*%alpha_diff
          f2<- trace(SigmaInv_beta_0%*%Sigma_beta) + crossprod(beta_diff , SigmaInv_beta_0)%*%beta_diff
          Logp<-Logp -0.5*( (pn+qn)*log(2*pi) - determinant(Sigma_alpha_0, logarithm=T)$modulus[1] - determinant(Sigma_beta_0, logarithm=T)$modulus[1] + f1 + f2)
          
          #the E(q(alpha, beta, z)) part
          #------------------------------
          E_log_pz <-sum(log(p[not_obs_vec]^p[not_obs_vec]) + log( (1-p[not_obs_vec])^(1-p[not_obs_vec]) ))
          E_log_p_alpha_beta <- -0.5*( const1 + determinant(Sigma_alpha, logarithm=T)$modulus[1] + determinant(Sigma_beta, logarithm=T)$modulus[1] )
          
          Log_ml <- Logp - E_log_p_alpha_beta[1] - E_log_pz
          
          diff_Log_ml <- abs(Log_ml-old_Log_ml)
          old_Log_ml <- Log_ml
          
          if (Print==TRUE)
          {cat("Iteration: ", its,", Log_ml <- ", round(Log_ml,digits<-6), ", alpha <- ", round(alpha,digits<-3) , ", beta <- ", round(beta,digits<-3), "\n")}
     }
     list(alpha=alpha, beta=beta, Sigma_alpha=Sigma_alpha, Sigma_beta=Sigma_beta, occup_p=c(p), Log_mla=Log_ml, Breakcounter=Breakcounter)
}

negLL_naive = function(param)
{
     #The following function calculates the negative loglikelihood using the 'naive' loglikelihood defined in Moreno and Lele 2010
     beta = param[1:dim(X)[2]]
     psi = as.vector(1/(1+exp(-X %*% beta))) #logistic link function used
     
     logL = rep(NA,n)
     yvec = apply(y,1,max)
     
     terms = psi^yvec * (1-psi)^(1-yvec)
     logL = log( prod(terms)   )
     (-1)*sum(logL)
}

#Penalized likelihood method
negLL_pen = function(param)
{
     #The following function does the penalized loglikelihood estimation of Moreno and Lele 2010
     
     beta = param[1:dim(X)[2]]
     psi = as.vector(1/(1+exp(-X %*% beta))) #logistic link function used
     alpha = param[(dim(X)[2]+1):(dim(X)[2]+dim(W)[3])]
     p = matrix(nrow=n, ncol=J)
     for (j in 1:J)
     {
          p[, j] = c(1/(1+exp(-W[,j,] %*% alpha)))
     }
     logL = rep(NA,n)
     for (i in 1:n) {
          yvec = y[i, jind[i,]]
          pvec = p[i, jind[i,]]
          terms = pvec^yvec * (1-pvec)^(1-yvec)
          logL[i] = log( psi[i] * prod(terms)  + ifelse(ysum[i]==0,1,0)*(1-psi[i])  )
     }
     
     penalty = lambda_0 * sum(abs(beta-beta_naive))
     logl_pen = sum(logL) - penalty
     return(-logl_pen)
}

vb_Designs<-function(W, X, y)
{
     #create the required 'response' and 'regressor matrices'
     #using all of the X and W data
     #the output is stored as a named list
     
     #create the Y matrix that will be used
     Y<-matrix(na.omit(matrix(t(y), ncol=1)))
     #col.names(Y)<-names(y)
     pres_abs <- apply(y,1,max,na.rm=T) #check if this will work for NA's
     
     #create the W matrix
     W.temp<-NULL
     nv<-length(W)
     
     for (i in 1:nv){W.temp<-cbind(W.temp, W[[i]])}
     #print(W.temp)
     
     nvisits<-apply(W[[1]],1,function(x){length(na.omit(x))})
     n<-length(nvisits)
     
     W.out<-NULL
     for (i in 1:n){W.out<-rbind(W.out, matrix( c(na.omit(W.temp[i,])), nrow=nvisits[i]) )}
     colnames(W.out)<-names(W)
     
     list(Y=as.data.frame(Y), X=as.data.frame(X), W=as.data.frame(W.out),
          Names=c( colnames(X), colnames(W.out)), nvisits=nvisits,
          pres_abs=pres_abs)
}

sink("occ.txt")
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
    logit(pz[i]) <- beta0 + beta1*X[i]  #covariate = AET_div_PET_s
    
    for (j in 1:nvisits[i]) #loop over the number of visits per site 
{
    #observation process
    Y[i,j] ~ dbern(py[i,j])
    py[i,j] <- zb[i]*pd[i,j]
    logit(pd[i,j])<- alpha0+alpha1*W[i,j,2]  
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
require(unmarked) #load the 'unmarked' package
require(jagsUI)
require(coda)
require(xtable)

#load("sabap_data_Plos2016.RData")#load the R data image. (raw data)
load("S1_Data.RData")#load the R data image. (raw data)


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

sink.T=TRUE
use.species<-1:5
bayes.species<-NULL
sink1<-"Analysis1_3_to_5_2012_Bayes2.txt"
what.year=112

if (sink.T==TRUE){sink(sink1)}
sink(sink1, append=TRUE)

#Do the analysis for each of the species in turn
for (ispecies in use.species)  
{
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
     OMdat <- unmarkedFrameOccu(y=Y, siteCovs = site.dat, obsCovs=obs.covs)
     #------------------------------------------------------------------
     
     Nvisits<-apply(y[choose.rows,],1,function(x){length(na.omit(x))})

     #------------------------------------------------------------------
     #The MLE fit
     #------------------------------------------------------------------
     f8<- occu( ~ nspp  ~AET_div_PET_s , data=OMdat,engine="C")
     f8.summary<-try(summary(f8),T)
     pvalues<-sum(c(sum(f8.summary[["state"]]["P(>|z|)"] <.05, na.rm=T), sum(f8.summary[["det"]]["P(>|z|)"]<0.05, na.rm=T) ))
     
     if (pvalues==4) #if the MLE can be performed then continue
     {
          bayes.species<-c(bayes.species,ispecies) #append
          cat("\n ....................................")
          cat("\n ", as.character(sp[ispecies]), "- Speciesid = ", ispecies)
          cat("\n ....................................")
          
          cat("\n Number of sites = ", length(len.y[len.y!=0]))
          cat("\n Max surveys = ", maxSurveys)
          print( table(len.y[len.y!=0])[1:10] )
               
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
          
          #my code does not allow one to fit the simple model! i.e. 
          #an intercept only model cannot be fitted
          #The call to the Laplace algorithm VB code
          vb_a1<-vb_model2_la(formula=form1, design_mats=design_mats, alpha_0=alpha_0, beta_0=beta_0, Sigma_alpha_0=Sigma_alpha_0, Sigma_beta_0=Sigma_beta_0,LargeSample=F, epsilon=1e-8, F)

          #Store the results
          vbouts<-rbind(cbind(c(vb_a1$beta), sqrt(diag(vb_a1$Sigma_beta))), cbind( c(vb_a1$alpha), sqrt(diag(vb_a1$Sigma_alpha)) ))
          vbouts<-cbind(vbouts,vbouts[,1]/vbouts[,2])
          colnames(vbouts)<-c("Est","SD","t")
          
          cat("\n VB results \n")
          print(vbouts)
          
          cat("\n MLE results \n")
          print(f8)
          cat("\n")       
          mle.fit<-f8

          #------------------------------------------------------------------
          # The Bayesian fit
          #------------------------------------------------------------------

          #inputs for the Baysian method
          niter <- 75000
          nburn <- 25000
          nthin <- 1
          nchains <- 3
          nkeep<-(niter-nburn)*nchains
          
          Y.mat<-Y
          X.use<-Xmat[,6]
          
          n.use<-length(choose.rows)
          
          W = array(dim=c(n.use,5,2)) #number of observations ; numbe of visits; number of covariates + 1 (for the intercept)
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
          
          #------------------------------------------------------------------
          # Input data for WinBUGS
          #------------------------------------------------------------------
          zinfo<-function(y)
          {
               #if the species was observed at a site then we set z==1
               z<-matrix(NA, ncol=1, nrow=NROW(y))
               zobs<-apply(y,1,function(x) max(x, na.rm=TRUE))
               z[which(zobs==1)]<-1
               return(z)
          }
          occ.data <- list(Y=Y.mat, X=X.use, W=W, n=dim(Y.mat)[1], zb=c(zinfo(Y.mat)) , nvisits=nvisits)
               
          #------------------------------------------------------------------
          # Initial values
          #------------------------------------------------------------------
          
          init.z<-function(y)
          {
               zb<-apply(y,1,function(x) max(x, na.rm=TRUE))
               
               if (length(which(zb==1))>0)
               {
                    zb[which(zb==1)]<-NA
               }
               return(zb)
          }
          
          inits <- function()
          {
               list("zb"=c(init.z(Y.mat)),"beta0"=rnorm(1),"beta1"=rnorm(1),
                    "alpha0"=rnorm(1),"alpha1"=rnorm(1))
          }
          
          inits1=inits(); inits2=inits(); inits3=inits()
          Inits=list(inits1,inits2,inits3)
          
          #------------------------------------------------------------------
          # Paramaters to be monitored       
          #------------------------------------------------------------------
          parameters <- c("beta0", "beta1", "alpha0", "alpha1", "occupied",  "zb", "pz") 
          
          #jagsUI often has initial value problems!!!
          occ.jags <- jags(occ.data, inits=Inits, parameters, "occ.txt",
                                n.iter = niter, n.burnin = nburn, n.thin = nthin, n.chains = nchains,
                                DIC = FALSE, parallel=TRUE)
          
          #occ.jags
          
          alpha0 = occ.jags$sims.list$alpha0
          alpha1 = occ.jags$sims.list$alpha1
          
          beta0 = occ.jags$sims.list$beta0
          beta1 = occ.jags$sims.list$beta1
          
          occupied = occ.jags$sims.list$occupied
          zb<-occ.jags$sims.list$zb
         
          RegParSims<-cbind(alpha0, alpha1, beta0,beta1)
          Means<-apply(RegParSims,2,mean)
          SDs<-apply(RegParSims,2,sd)
          
          occ.mcmc<-mcmc(data=RegParSims)
          #names(occ.jags)
          
          MCSE<-batchSE(occ.mcmc)
          mcmctab<-cbind(Means, SDs, MCSE)
          colnames(mcmctab)<-c("Mean","Std","MCSE")
          print(xtable(mcmctab, digits=3))
          
          r1<-cbind(coef(mle.fit,"det"),SE(mle.fit,"det"))
          r2<-cbind(coef(mle.fit,"state"), SE(mle.fit,"state"))
          r1r2<-rbind(r1,r2)
          #r1r2
          rownames(r1r2)<-c("Int (detection)","nspp","Int (Occupancy)","AET_div_PET_s")
          colnames(r1r2)<-c("Est", "SE")
          #r1r2
          v1<-cbind(vb_a1$alpha, sqrt(diag(vb_a1$Sigma_alpha)))
          v2<-cbind(vb_a1$beta, sqrt(diag(vb_a1$Sigma_beta)))
          v1v2<-cbind(rbind(v1,v2))
          #v1v2
          colnames(v1v2)<-c("Mean", "Std")
          RV<-cbind(r1r2,v1v2, mcmctab)
          #print(RV)
          print(xtable(RV, digits=3, align=c("|l|","c|","c|","c|","c|","c","c","c|")))
          
          # barplot of the distribution of the number of occupied sites
          #-----------------------------------------------------------------------
          
          occupiedmcmc<-occupied
          Occpla<-vb_a1$occup_p
          Nsim=1000000
          occ_n<-matrix(1,ncol=Nsim, nrow=length(Occpla))
          zmaybe<-which(apply(Y.mat,1,function(x)max(x,na.rm=T))==0) #only sample at these places
          
          for (i in zmaybe){occ_n[i,]<-rbinom(Nsim, p=Occpla[i], size=1)}
          
          occupiedla<-apply(occ_n, 2, sum)
          
          #pdf(paste("SABAP_2013_Species_",ispecies,"year",what.year,"_PAO.pdf",sep=""))
          par(mfrow=c(1,1), mar=c(4,3,2,2))
          t1<-table(occupiedmcmc)/nkeep
          t1.numeric<-as.numeric(dimnames(t1)$occupiedmcmc)
          t2<-table(occupiedla)/Nsim
          t2.numeric<-as.numeric(dimnames(t2)$occupiedla)
          min.labels<-min(t1.numeric, t2.numeric)
          max.labels<-max(t1.numeric, t2.numeric)
          range.labels<-min.labels:max.labels
          
          mcmc.prob=t1[match(range.labels, t1.numeric)]
          la.prob=t2[match(range.labels, t2.numeric)]
          la.prob[is.na(la.prob)]<-0
          data<-rbind(mcmc.prob,la.prob)
          colnames(data) <- paste(range.labels)
          barplot(data, beside=TRUE, space=c(0,1), las=1, col=c("black","white"), xlab="Number of occupied sites")
          legend("topright", c("MCMC","VB"), fill=c("black","white"), bty="n")
          title(main = "Posterior distribution of the \n number of occupied sites", font.main = 4)
          #dev.off()
          
          Summary.stats<-rbind(
               c(mean(occupiedmcmc),
                 sd(occupiedmcmc),
                 c(quantile(occupiedmcmc, 0.025), quantile(occupiedmcmc, 0.975) )),
               
               c(mean(occupiedla),
                 sd(occupiedla),
                 c(quantile(occupiedla, 0.025), quantile(occupiedla, 0.975) )) )
          
          colnames(Summary.stats)<-c("Mean","Std", "2.5%","97.5%")
          print(xtable(Summary.stats, digits=2))
          print(Summary.stats)
          
          #accuracy statistic for the example
          #-----------------------------------------------------------------------
          t1<-table(occupiedmcmc)/nkeep
          t1.numeric<-as.numeric(dimnames(t1)$occupiedmcmc)
          t2<-table(occupiedla)/Nsim
          t2.numeric<-as.numeric(dimnames(t2)$occupiedla)
          min.labels<-min(t1.numeric, t2.numeric)
          max.labels<-max(t1.numeric, t2.numeric)
          range.labels<-min.labels:max.labels
          
          mcmc.prob=t1[match(range.labels, t1.numeric)]
          la.prob=t2[match(range.labels, t2.numeric)]
          la.prob[is.na(la.prob)]<-0
          mcmc.prob[is.na(mcmc.prob)]<-0
          data<-rbind(mcmc.prob,la.prob)
          cat("\n The accuracy stat = ", 1-.5*sum(abs(data[1,]-data[2,])) )#the % overlap
          #-----------------------------------------------------------------------
          
          #the MCMC posterior distribution of the parameters as well as the Vb posterior distribution
          #-------------------------------------------------------------------------------------------
          #pdf(paste("SABAP_2013_Species_",ispecies,"_year_",what.year,".pdf",sep=""))
          par(mfrow=c(2,2), mar=c(4,4,4,1))
          
          #plots for alpha0
          p1=density(alpha0)
          xx=seq(from=(min(p1$x)-.25), to=(max(p1$x)+.25), length.out=1000)
          d1=dnorm(xx, mean=vb_a1$alpha[1], sd=sqrt(vb_a1$Sigma_alpha[1,1]))
          hip1=max(p1$y)
          hip2=max(d1)         
          h<-hist(alpha0, plot=F, breaks=50)
          plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab=expression(alpha[0]), main="", ylab="Density")
          lines(xx, d1, col="red")
          legend("topleft", c("VB"), lwd=c(1,2), col=c("red","blue"),bty="n")
          
          #plots for alpha1
          p1=density(alpha1)
          xx=seq(from=min(p1$x), to=max(p1$x), length.out=1000)
          d1=dnorm(xx, mean=vb_a1$alpha[2], sd=sqrt(vb_a1$Sigma_alpha[2,2]))
          hip1=max(p1$y)
          hip2=max(d1)         
          h<-hist(alpha1, plot=F, breaks=50)
          plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab="Number of observed Species", main="", ylab="Density")
          lines(xx, d1, col="red")
          legend("topleft", c("VB"), lwd=c(1,2), col=c("red","blue"),bty="n")
          
          #plots for beta0
          p1=density(beta0)
          xx=seq(from=min(p1$x), to=max(p1$x), length.out=1000)
          d1=dnorm(xx, mean=vb_a1$beta[1], sd=sqrt(vb_a1$Sigma_beta[1,1]))
          hip1=max(p1$y)
          hip2=max(d1)         
          h<-hist(beta0, plot=F, breaks=50)
          plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab=expression(beta[0]), main="", ylab="Density")
          lines(xx, d1, col="red")
          legend("topleft", c("VB"), lwd=c(1,2), col=c("red","blue"),bty="n")
          
          #plots for beta1
          p1=density(beta1)
          xx=seq(from=min(p1$x), to=max(p1$x), length.out=1000)
          d1=dnorm(xx, mean=vb_a1$beta[2], sd=sqrt(vb_a1$Sigma_beta[2,2]))
          hip1=max(p1$y)
          hip2=max(d1)         
          h<-hist(beta1, plot=F, breaks=50)
          plot(h$breaks,c(h$density,0),type="s",col="black", ylim=c(0,max(hip1,hip2)), xlab="AETdivPETs", main="", ylab="Density")
          lines(xx, d1, col="red")
          legend("topleft", c("VB"), lwd=c(1,2), col=c("red","blue"),bty="n")
          
          #dev.off()
          
          Compare<-cbind(
               rbind(c(mean(beta0),sd(beta0)),
                     c(mean(beta1),sd(beta1)),
                     c(mean(alpha0),sd(alpha0)),
                     c(mean(alpha1),sd(alpha1)) ) ,
               vbouts[,1:2])
          colnames(Compare)<-c("MCMC Est","MCMC sd","VB Est","VB sd")
          
          #sink("GregDataAnalysis_estimations_2012.txt")
          cat("\n Occupancy process \n")
          Compare[1:2,]
          cat("\n Detection process \n")
          Compare[3:4,]
          #sink()
          
          file.save<-paste("Analysis_3_to_5_2012_Bayes_Species_",ispecies,"_year_",what.year,".RData",sep="")
          save.image(file.save) 
     }
     #-----------------------------------------------------------------------
}
cat("\n DONE RUNNING THE ANALYSIS! \n")

if (sink.T==T){sink()}
