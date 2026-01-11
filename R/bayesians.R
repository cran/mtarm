#'
#' @title Out-of-sample predictive accuracy measures
#' @description Computes out-of-sample predictive accuracy measures for one or more fitted models
#' of the same class, based on their predictive distributions.
#' @param ... one or more fitted model objects of the same class.
#' @param newdata a data frame containing the future values of the output series required to evaluate
#' predictive performance.
#' @param n.ahead a positive integer specifying the number of forecast steps ahead to use in the
#' predictive performance evaluation.
#' @export out.of.sample
out.of.sample <- function(...,newdata,n.ahead) {
  UseMethod("out.of.sample")
}
#'
#' @title Deviance Information Criterion (DIC)
#' @description Computes the Deviance Information Criterion (DIC), an adjusted
#' within-sample measure of predictive accuracy, for models estimated using Bayesian methods.
#' @param ... one or more fitted model objects of the same class.
#' @return A numeric matrix containing the DIC values corresponding to each fitted object supplied in \code{...}.
#' @references Spiegelhalter D.J., Best N.G., Carlin B.P. and Van Der Linde A. (2002) Bayesian Measures of Model Complexity and Fit.
#'             Journal of the Royal Statistical Society Series B (Statistical Methodology), 64(4), 583–639.
#' @references Spiegelhalter D.J., Best N.G., Carlin B.P. and Van der Linde A. (2014). The deviance information criterion:
#'             12 years on. Journal of the Royal Statistical Society Series B (Statistical Methodology), 76(3), 485–493.
#' @export DIC
DIC <- function(...) {
  UseMethod("DIC")
}
#' @title Watanabe-Akaike or Widely Available Information Criterion (WAIC)
#' @description Computes Watanabe-Akaike or Widely Available Information Criterion (WAIC), an adjusted
#' within-sample measure of predictive accuracy, for models estimated using Bayesian methods.
#' @param ... one or more fitted model objects of the same class.
#' @return A numeric matrix containing the WAIC values corresponding to each fitted object supplied in \code{...}.
#' @references Watanabe S. (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in
#'             Singular Learning Theory. The Journal of Machine Learning Research, 11, 3571–3594.
#' @export WAIC
WAIC <- function(...) {
  UseMethod("WAIC")
}
#'
#' @title Deviance Information Criterion (DIC) for objects of class \code{mtar}
#' @description This function computes the Deviance Information Criterion (DIC) for objects of class \code{mtar}.
#' @param ...	one or several objects of the class \emph{mtar}.
#' @return A numeric matrix containing the DIC values corresponding to each \emph{mtar} object in the input.
#' @method DIC mtar
#' @export
#' @seealso \link{WAIC}
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1a <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'               subset={Date<="2016-03-14"}, dist="Student-t",
#'               ars=ars(nregim=3,p=c(1,1,2)), n.burnin=2000, n.sim=3000,
#'               n.thin=2)
#' fit1b <- update(fit1a,dist="Slash")
#' fit1c <- update(fit1a,dist="Laplace")
#' DIC(fit1a,fit1b,fit1c)
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2a <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'               subset={Date<="2009-04-04"}, dist="Laplace",
#'               ars=ars(nregim=3,p=5), n.burnin=2000, n.sim=3000, n.thin=2)
#' fit2b <- update(fit2a,dist="Slash")
#' fit2c <- update(fit2a,dist="Student-t")
#' DIC(fit2a,fit2b,fit2c)
#'
#' ###### Example 3: Temperature, precipitation, and two river flows in Iceland
#' data(iceland.rf)
#' fit3a <- mtar(~ Jokulsa + Vatnsdalsa | Temperature | Precipitation,
#'               data=iceland.rf, subset={Date<="1974-12-21"}, row.names=Date,
#'               ars=ars(nregim=2,p=15,q=4,d=2), n.burnin=2000, n.sim=3000,
#'               n.thin=2, dist="Slash")
#' fit3b <- update(fit3a,dist="Laplace")
#' fit3c <- update(fit3a,dist="Student-t")
#' DIC(fit3a,fit3b,fit3c)
#' }
#'
DIC.mtar <- function(...){
  # Internal function to compute minus twice the log-likelihood
  mtlogLik <- function(dist,y,X,beta,Sigma,delta,log,nu){
      # Compute residuals
      resu <- y-X%*%beta
      # Case of skewed distributions
      if(dist %in% c("Skew-normal","Skew-Student-t")){
         # Auxiliary matrices for skew distributions
         A <- chol2inv(chol(Sigma + as.vector(delta^2)*diag(k)))
         muv <- resu%*%(A*matrix(delta,k,k,byrow=TRUE))
         Sigmav <- diag(k) - matrix(delta,k,k)*A*matrix(delta,k,k,byrow=TRUE)
         out <- colSums(t(resu)*tcrossprod(A,resu))
         # Skew-normal likelihood contribution
         if(dist=="Skew-normal"){
            if(k > 1) sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pmvnorm(lower=-x,upper=rep(Inf,k),sigma=Sigmav)))
            else sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pnorm(-x/sqrt(Sigmav),lower.tail=FALSE)))
            out <- 2^k*exp(-out/2)*(det(A)^(1/2))*sum0/(2*pi)^(k/2)
         }
         # Skew-Student-t likelihood contribution
         if(dist=="Skew-Student-t"){
            muv <- muv*matrix(sqrt((nu + k)/(nu + out)),nrow(muv),k)
            if(k > 1) sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pmvt(lower=-x,upper=rep(Inf,k),sigma=Sigmav,df=round(nu,0)+k)))
            else sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pt(-x/sqrt(Sigmav),df=nu,lower.tail=FALSE)))
            out <- 2^k*(1 + out/nu)^(-(nu+k)/2)*(det(A)^(1/2))*gamma((nu+k)/2)*sum0/((nu*pi)^(k/2)*gamma(nu/2))
         }
      }else{# Case of symmetric distributions
           out <- colSums(t(resu)*tcrossprod(chol2inv(chol(Sigma)),resu))
           # Distribution-specific likelihood kernels
           out <- switch(dist,
                         "Gaussian"={exp(-out/2)},
                         "Laplace"={besselK(sqrt(out)/2,(2-k)/2)*out^((2-k)/4)/(2^((k+2)/2))},
                         "Student-t"={(1 + out/nu)^(-(nu+k)/2)*gamma((nu+k)/2)/((nu/2)^(k/2)*gamma(nu/2))},
                         "Contaminated normal"={nu[1]*exp(-out*nu[2]/2)*nu[2]^(k/2) + (1-nu[1])*exp(-out/2)},
                         "Slash"={ifelse(out==0,nu/(nu+k),gamma((k+nu)/2)*(nu/2)*pgamma(1,shape=(k+nu)/2,rate=out/2)/((out/2)^((k+nu)/2)))},
                         "Hyperbolic"={besselK(nu*sqrt(1+out),(2-k)/2)*(1+out)^((2-k)/4)*nu^(k/2)/besselK(nu,1)})
           # Normalizing constant
           out <- out/((2*pi)^(k/2)*det(Sigma)^(1/2))
      }
      # Convert to deviance (-2 log-likelihood)
      out <- -2*sum(log(out))
      # Optional Jacobian term if log-transformation is used
      if(log) out <- out + 2*sum(y)
      return(out)
  }
  # Collect fitted mtar objects
  another <- list(...)
  call. <- match.call()
  # Initialize output container
  out <- matrix(0,length(another),1)
  outnames <- vector()
  # Loop over each fitted model
  for(l in 1:(length(another))){
      n.sim <- another[[l]]$n.sim
      dist <- another[[l]]$dist
      Dbar <- vector()
      k <- ncol(another[[l]]$data[[1]]$y)
      Dbarv <- vector()
      # Indices for threshold series if multiple regimes
      if(another[[l]]$regim > 1) lims <- (another[[l]]$ps+1):(length(another[[l]]$threshold.series))
      # Loop over MCMC iterations
      for(j in 1:n.sim){
          Dbar <- 0
          # Determine regimes at iteration j
          if(another[[l]]$regim > 1){
             Z <- another[[l]]$threshold.series[lims-another[[l]]$chains$h[j]]
             regs <- cut(Z,breaks=c(-Inf,sort(another[[l]]$chains$thresholds[,j]),Inf),labels=FALSE)
          }else regs <- matrix(1,nrow(another[[l]]$data[[1]]$y),1)
          # Loop over regimes
          for(i in 1:another[[l]]$regim){
              # Variable selection indicator (if SSVS is used)
              if(another[[l]]$ssvs) zetai <- another[[l]]$chains[[i]]$zeta[,j] else zetai <- rep(1,ncol(another[[l]]$data[[i]]$X))
              # Regression coefficients and covariance matrix
              betai <- matrix(another[[l]]$chains[[i]]$location[zetai==1,((j-1)*k + 1):(j*k)],sum(zetai),k)
              Sigmai <- matrix(another[[l]]$chains[[i]]$scale[,((j-1)*k + 1):(j*k)],k,k)
              # Observations belonging to the current regime
              places <- regs==i
              yy <- matrix(another[[l]]$data[[i]]$y[places,],sum(places),k)
              XX <- matrix(another[[l]]$data[[i]]$X[places,zetai==1],sum(places),sum(zetai==1))
              # Accumulate deviance
              Dbar <- switch(another[[l]]$dist,
                             "Skew-Student-t"={Dbar + mtlogLik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])},
                             "Skew-normal"={Dbar + mtlogLik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log)},
                             "Student-t"=,"Hyperbolic"=,"Slash"=,"Contaminated normal"={Dbar + mtlogLik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])},
                             "Gaussian"=,"Laplace"={Dbar + mtlogLik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,log=another[[l]]$log)})
          }
          Dbarv <- c(Dbarv,Dbar)
      }
      # Compute deviance at posterior mean parameters
      Dhat <- 0
      # Posterior means of extra parameters
      if(another[[l]]$dist %in% c("Skew-Student-t","Student-t","Hyperbolic","Slash","Contaminated normal"))
         extra <- rowMeans(matrix(another[[l]]$chains$extra,ncol=n.sim))
      # Posterior means of skewness parameters
      if(another[[l]]$dist %in% c("Skew-Student-t","Skew-normal"))
         delta <- rowMeans(matrix(another[[l]]$chains$delta,ncol=n.sim))
      # Regime classification using posterior mean thresholds
      if(another[[l]]$regim > 1){
         h <- as.integer(round(mean(another[[l]]$chains$h)))
         thresholds <- rowMeans(matrix(another[[l]]$chains$thresholds,ncol=n.sim))
         Z <- another[[l]]$threshold.series[lims-h]
         regs <- cut(Z,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
      }else regs <- matrix(1,nrow(another[[l]]$data[[1]]$y),1)
      # Loop over regimes using posterior mean parameters
      for(i in 1:another[[l]]$regim){
          location <- matrix(0,nrow(another[[l]]$chains[[i]]$location),k)
          scale <- matrix(0,k,k)
          for(j in 1:k){
              location[,j] <- rowMeans(matrix(another[[l]]$chains[[i]]$location[,seq(j,n.sim*k,by=k)],nrow(location),n.sim))
              scale[,j] <- rowMeans(matrix(another[[l]]$chains[[i]]$scale[,seq(j,n.sim*k,by=k)],k,n.sim))
          }
          # Apply variable selection at posterior mean (if SSVS is used)
          if(another[[l]]$ssvs){
             zeta <- rowMeans(another[[l]]$chains[[i]]$zeta)
             location <- matrix(location[zeta>0.5,],sum(zeta>0.5),k)
          }else zeta <- rep(1,ncol(another[[l]]$data[[i]]$X))
          # Data for current regime
          places <- regs==i
          yy <- matrix(another[[l]]$data[[i]]$y[places,],sum(places),k)
          XX <- matrix(another[[l]]$data[[i]]$X[places,zeta>0.5],sum(places),sum(zeta>0.5))
          # Accumulate deviance at posterior mean
          switch(another[[l]]$dist,
                 "Skew-Student-t"={Dhat <- Dhat + mtlogLik(dist=dist,y=yy,X=XX,beta=location,Sigma=scale,delta=delta,log=another[[l]]$log,nu=extra)},
                 "Skew-normal"={Dhat <- Dhat + mtlogLik(dist=dist,y=yy,X=XX,beta=location,Sigma=scale,delta=delta,log=another[[l]]$log)},
                 "Student-t"=,"Hyperbolic"=,"Slash"=,"Contaminated normal"={Dhat <- Dhat + mtlogLik(dist=dist,y=yy,X=XX,beta=location,Sigma=scale,log=another[[l]]$log,nu=extra)},
                 "Gaussian"=,"Laplace"={Dhat <- Dhat + mtlogLik(dist=dist,y=yy,X=XX,beta=location,Sigma=scale,log=another[[l]]$log)})
      }
      # Deviance Information Criterion
      out[l] <- Dhat + 2*(mean(Dbarv) - Dhat)
      outnames[l] <- as.character(call.[l+1])
  }
  # Assign names and return DIC values
  rownames(out) <- outnames
  colnames(out) <- "DIC"
  return(out)
}
#'
#' @title Watanabe-Akaike or Widely Available Information Criterion (WAIC) for objects of class \code{mtar}
#' @description This function computes the Watanabe-Akaike or Widely Available Information Criterion (WAIC),
#' for objects of class \code{mtar}.
#' @param ...	one or several objects of the class \emph{mtar}.
#' @return A numeric matrix containing the WAIC values corresponding to each \emph{mtar} object in the input.
#' @method WAIC mtar
#' @export
#' @seealso \link{DIC}
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1a <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'               subset={Date<="2016-03-14"}, dist="Student-t",
#'               ars=ars(nregim=3,p=c(1,1,2)), n.burnin=2000, n.sim=3000,
#'               n.thin=2)
#' fit1b <- update(fit1a,dist="Slash")
#' fit1c <- update(fit1a,dist="Laplace")
#' WAIC(fit1a,fit1b,fit1c)
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2a <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'               subset={Date<="2009-04-04"}, dist="Laplace",
#'               ars=ars(nregim=3,p=5), n.burnin=2000, n.sim=3000, n.thin=2)
#' fit2b <- update(fit2a,dist="Slash")
#' fit2c <- update(fit2a,dist="Student-t")
#' WAIC(fit2a,fit2b,fit2c)
#'
#' ###### Example 3: Temperature, precipitation, and two river flows in Iceland
#' data(iceland.rf)
#' fit3a <- mtar(~ Jokulsa + Vatnsdalsa | Temperature | Precipitation,
#'               data=iceland.rf, subset={Date<="1974-12-21"}, row.names=Date,
#'               ars=ars(nregim=2,p=15,q=4,d=2), n.burnin=2000, n.sim=3000,
#'               n.thin=2, dist="Slash")
#' fit3b <- update(fit3a,dist="Laplace")
#' fit3c <- update(fit3a,dist="Student-t")
#' WAIC(fit3a,fit3b,fit3c)
#' }
#'
WAIC.mtar <- function(...){
  # Internal function to compute the likelihood contribution
  Lik <- function(dist,y,X,beta,Sigma,delta,log,nu){
      # Compute residuals
      resu <- y-X%*%beta
      # Case of skewed distributions
      if(dist %in% c("Skew-normal","Skew-Student-t")){
         A <- chol2inv(chol(Sigma + as.vector(delta^2)*diag(k)))
         muv <- resu%*%(A*matrix(delta,k,k,byrow=TRUE))
         Sigmav <- diag(k) - matrix(delta,k,k)*A*matrix(delta,k,k,byrow=TRUE)
         out <- colSums(t(resu)*tcrossprod(A,resu))
         # Skew-normal likelihood
         if(dist=="Skew-normal"){
            if(k > 1) sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pmvnorm(lower=-x,upper=rep(Inf,k),sigma=Sigmav)))
            else sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pnorm(-x/sqrt(Sigmav),lower.tail=FALSE)))
            out <- 2^k*exp(-out/2)*(det(A)^(1/2))*sum0/(2*pi)^(k/2)
         }
         # Skew-Student-t likelihood
         if(dist=="Skew-Student-t"){
            muv <- muv*matrix(sqrt((nu + k)/(nu + out)),nrow(muv),k)
            if(k > 1) sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pmvt(lower=-x,upper=rep(Inf,k),sigma=Sigmav,df=round(nu,0)+k)))
            else sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pt(-x/sqrt(Sigmav),df=nu,lower.tail=FALSE)))
            out <- 2^k*(1 + out/nu)^(-(nu+k)/2)*(det(A)^(1/2))*gamma((nu+k)/2)*sum0/((nu*pi)^(k/2)*gamma(nu/2))
         }
      }else{# Case of symmetric distributions
           out <- colSums(t(resu)*tcrossprod(chol2inv(chol(Sigma)),resu))
           # Distribution-specific likelihood kernels
           out <- switch(dist,
                         "Gaussian"={exp(-out/2)},
                         "Laplace"={besselK(sqrt(out)/2,(2-k)/2)*out^((2-k)/4)/(2^((k+2)/2))},
                         "Student-t"={(1 + out/nu)^(-(nu+k)/2)*gamma((nu+k)/2)/((nu/2)^(k/2)*gamma(nu/2))},
                         "Contaminated normal"={nu[1]*exp(-out*nu[2]/2)*nu[2]^(k/2) + (1-nu[1])*exp(-out/2)},
                         "Slash"={ifelse(out==0,nu/(nu+k),gamma((k+nu)/2)*(nu/2)*pgamma(1,shape=(k+nu)/2,rate=out/2)/((out/2)^((k+nu)/2)))},
                         "Hyperbolic"={besselK(nu*sqrt(1+out),(2-k)/2)*(1+out)^((2-k)/4)*nu^(k/2)/besselK(nu,1)})
           # Normalizing constant
           out <- out/((2*pi)^(k/2)*det(Sigma)^(1/2))
      }
      # Adjustment for log-transformed responses
      if(log) out <- out*apply(matrix(y,nrow(X),k),1,function(x) prod(exp(-x)))
      return(out)
  }
  # Collect fitted mtar objects
  another <- list(...)
  call. <- match.call()
  # Initialize output container
  out <- matrix(0,length(another),1)
  outnames <- vector()
  # Loop over each fitted model
  for(l in 1:(length(another))){
      n.sim <- another[[l]]$n.sim
      dist <- another[[l]]$dist
      # Containers for pointwise likelihood averages
      Dbar <- matrix(0,nrow(another[[l]]$data[[1]]$y),1)
      Dbarlog <- matrix(0,nrow(another[[l]]$data[[1]]$y),1)
      # Number of response variables
      k <- ncol(another[[l]]$data[[1]]$y)
      # Indices for threshold series if multiple regimes
      if(another[[l]]$regim > 1) lims <- (another[[l]]$ps+1):(length(another[[l]]$threshold.series))
      # Loop over MCMC iterations
      for(j in 1:n.sim){
          # Loop over regimes
          for(i in 1:another[[l]]$regim){
              # Variable selection indicator (if SSVS is used)
              if(another[[l]]$ssvs) zetai <- another[[l]]$chains[[i]]$zeta[,j] else zetai <- rep(1,ncol(another[[l]]$data[[i]]$X))
              # Regression coefficients and covariance matrix
              betai <- matrix(another[[l]]$chains[[i]]$location[zetai==1,((j-1)*k + 1):(j*k)],sum(zetai),k)
              Sigmai <- matrix(another[[l]]$chains[[i]]$scale[,((j-1)*k + 1):(j*k)],k,k)
              # Determine regimes at iteration j
              if(another[[l]]$regim > 1){
                 Z <- another[[l]]$threshold.series[lims-another[[l]]$chains$h[j]]
                 regs <- cut(Z,breaks=c(-Inf,sort(another[[l]]$chains$thresholds[,j]),Inf),labels=FALSE)
              }else regs <- matrix(1,nrow(another[[l]]$data[[i]]$y),1)
              # Observations belonging to the current regime
              places <- regs==i
              yy <- matrix(another[[l]]$data[[i]]$y[places,],sum(places),k)
              XX <- matrix(another[[l]]$data[[i]]$X[places,zetai==1],sum(places),sum(zetai))
              # Likelihood contribution for current iteration and regime
              tempi <- switch(dist,
                              "Skew-Student-t"={Lik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])},
                              "Skew-normal"={Lik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log)},
                              "Student-t"=,"Hyperbolic"=,"Slash"=,"Contaminated normal"={Lik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])},
                              "Gaussian"=,"Laplace"={Lik(dist=dist,y=yy,X=XX,beta=betai,Sigma=Sigmai,log=another[[l]]$log)})
              # Accumulate pointwise likelihoods
              Dbar[places] <- Dbar[places] + tempi/n.sim
              Dbarlog[places] <- Dbarlog[places] + log(tempi)/n.sim
          }
      }
      # Compute WAIC components
      a <- sum(unlist(lapply(Dbar,function(x) sum(log(x)))))
      b <- sum(unlist(lapply(Dbarlog,sum)))
      # Widely Applicable Information Criterion
      out[l] <- -2*a + 4*(a-b)
      outnames[l] <- as.character(call.[l+1])
  }
  # Assign names and return WAIC values
  rownames(out) <- outnames
  colnames(out) <- "WAIC"
  return(out)
}
#'
#' @title Coercion of \code{mtar} objects to \code{mcmc} objects
#' @description This method converts an object of class \code{mtar} into a list of
#' \code{mcmc} objects, each corresponding to a Markov chain produced during
#' Bayesian estimation.
#' @param x an object of class \code{mtar} obtained from a call to \code{mtar()}.
#' @param ... additional arguments passed to specific coercion methods.
#' @return A list of \code{mcmc} objects containing the posterior simulation draws
#' generated by the \code{mtar()} routine.
#' @method as.mcmc mtar
#' @seealso \code{\link[coda]{as.mcmc}}
#' @export
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1 <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'              subset={Date<="2016-03-14"}, dist="Student-t",
#'              ars=ars(nregim=3,p=c(1,1,2)), n.burnin=2000, n.sim=3000,
#'              n.thin=2)
#' fit1.mcmc <- coda::as.mcmc(fit1)
#' summary(fit1.mcmc)
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2 <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'              subset={Date<="2009-04-04"}, dist="Laplace",
#'              ars=ars(nregim=3,p=5), n.burnin=2000, n.sim=3000, n.thin=2)
#' fit2.mcmc <- coda::as.mcmc(fit2)
#' summary(fit2.mcmc)
#'
#' ###### Example 3: Temperature, precipitation, and two river flows in Iceland
#' data(iceland.rf)
#' fit3 <- mtar(~ Jokulsa + Vatnsdalsa | Temperature | Precipitation,
#'              data=iceland.rf, subset={Date<="1974-12-21"}, row.names=Date,
#'              ars=ars(nregim=2,p=15,q=4,d=2), n.burnin=2000, n.sim=3000,
#'              n.thin=2, dist="Slash")
#' fit3.mcmc <- coda::as.mcmc(fit3)
#' summary(fit3.mcmc)
#' }
as.mcmc.mtar <- function(x,...){
  # Number of response variables
  k <- ncol(x$data[[1]]$y)
  # Names of the components of the output series
  yn <- colnames(x$data[[1]]$y)
  # Number of MCMC simulations
  n.sim <- x$n.sim
  # Initialize output list
  out <- list()
  # MCMC metadata
  out$thin <- x$n.thin
  out$start <- x$n.burnin+1
  out$end <- max(seq(x$n.burnin + 1,length.out=x$n.sim,by=x$n.thin))
  out$sz <- x$n.sim
  out$regim <- x$ars$nregim
  # Containers for chains of location and scale parameters
  out$location <- list()
  out$scale <- list()
  if(x$ssvs) out$zeta <- list()
  # Loop over regimes
  for(r in 1:x$regim){
      # Posterior inclusion probabilities for SSVS if present
      if(x$ssvs) zeta <- rowMeans(x$chains[[r]]$zeta) else zeta <- rep(1,nrow(x$chains[[r]]$location))
      # Extract and reshape regression coefficients
      datos <- matrix(x$chains[[r]]$location[zeta>0.5,],nrow=sum(zeta>0.5),ncol=n.sim*k)
      out_ <- vector()
      # Organize coefficients by response variable
      for(i in 1:k){
          temp <- t(matrix(datos[,seq(i,n.sim*k,k)],nrow=nrow(datos),ncol=n.sim))
          colnames(temp) <- paste0(yn[i],":",rownames(x$chains[[r]]$location)[zeta>0.5])
          out_ <- cbind(out_,temp)
      }
      # Store location parameters as an mcmc object
      out$location[[r]] <- mcmc(out_,thin=x$n.thin,start=x$n.burnin+1)
      # Extract and reshape scale parameters
      datos <-  matrix(x$chains[[r]]$scale,nrow=k,ncol=n.sim*k)
      out_ <- vector()
      # Store upper triangular elements of scale parameters
      for(i in 1:k){
          for(j in i:k){
              temp <- t(matrix(datos[i,seq(j,n.sim*k,k)],nrow=1,ncol=n.sim))
              colnames(temp) <- paste0(yn[i],".",yn[j])
              out_ <- cbind(out_,temp)
          }
      }
      # Store scale parameters as an mcmc object
      out$scale[[r]] <- mcmc(out_,thin=x$n.thin,start=x$n.burnin+1)
      # Process SSVS indicators if present
      if(x$ssvs){
         # Indices of deterministic and lagged components
         quien <- x$deterministic + seq(k,x$ars$p[r]*k,k)
         namezeta <- paste0("OS.lag(",1:x$ars$p[r],")")
         # Add exogenous series lags if present
         if(x$ars$q[r]>0){
            quien <- c(quien,max(quien)+seq(x$r,x$ars$q[r]*x$r,x$r))
            namezeta <- c(namezeta,paste0("ES.lag(",1:x$ars$q[r],")"))
         }
         # Add threshold series lags if present
         if(x$ars$d[r]>0){
            quien <- c(quien,max(quien)+seq(1,x$ars$d[r],1))
            namezeta <- c(namezeta,paste0("TS.lag(",1:x$ars$d[r],")"))
         }
         temp <- matrix(x$chains[[r]]$zeta[quien,],length(quien),n.sim)
         rownames(temp) <- namezeta
         # Store SSVS indicators as an mcmc object
         out$zeta[[r]] <- mcmc(t(temp),thin=x$n.thin,start=x$n.burnin+1)
      }
  }
  # Store extra parameters as an mcmc object if present
  if(!(x$dist %in% c("Skew-normal","Gaussian","Laplace"))){
     datos <-  matrix(x$chains$extra,nrow=nrow(x$chains$extra),ncol=n.sim)
     out_ <- t(datos)
     if(nrow(datos)==1) colnames(out_) <- "nu" else colnames(out_) <- paste("nu",1:nrow(datos))
     out$extra <- mcmc(out_,thin=x$n.thin,start=x$n.burnin+1)
  }
  # Store skewness parameters as an mcmc object if present
  if(x$dist %in% c("Skew-normal","Skew-Student-t")){
     datos <-  matrix(x$chains$delta,nrow=nrow(x$chains$delta),ncol=n.sim)
     out_ <- t(datos)
     if(nrow(datos)==1) colnames(out_) <- "delta" else colnames(out_) <- paste("delta",1:nrow(datos))
     out$skewness <- mcmc(out_,thin=x$n.thin,start=x$n.burnin+1)
  }
  # Store thresholds and delay parameters as mcmc objects if present
  if(x$regim > 1){
     datos <-  matrix(x$chains$thresholds,nrow=nrow(x$chains$thresholds),ncol=n.sim)
     out_ <- t(datos)
     if(nrow(datos)==1) colnames(out_) <- "threshold" else colnames(out_) <- paste0("Threshold.",1:nrow(datos))
     out$thresholds <- mcmc(out_,thin=x$n.thin,start=x$n.burnin+1)
     out_ <- matrix(x$chains$h,nrow=n.sim,ncol=1)
     colnames(out_) <- "delay"
     out$delay <- mcmc(out_,thin=x$n.thin,start=x$n.burnin+1)
  }
  # Assign class and return
  class(out) <- "mtarmcmc"
  return(out)
}
#'
#'
#' @method print summtarmcmc
#' @export
print.summtarmcmc <- function(x,...,digits=max(3,getOption("digits")-2)){
  # Print basic MCMC information: iteration range, thinning, and sample size
  cat("\n\n Iterations = ",paste0(x$star,":",x$end))
  cat("\n Thinning interval = ",x$thin)
  cat("\n Sample size per chain = ",x$sz,"\n\n")
  # Print threshold parameters if the model has multiple regimes
  if(!is.null(x$thresholds)){
     cat("\nThresholds:\n")
     print(x$thresholds,digits=digits)
  }
  # Loop over regimes and print regime-specific parameters
  for(j in 1:x$regim){
      cat("\n\nRegime ",j,"\n")
      # Print location parameters
      cat("\n\nAutoregressive coefficients:\n")
      print(x$location[[j]],digits=digits)
      # Print scale parameters
      cat("\n\nScale parameter:\n")
      print(x$scale[[j]],digits=digits)
  }
  # Print skewness parameters if present
  if(!is.null(x$skewness)){
     cat("\n\nSkewness parameter:\n")
     print(x$skewness,digits=digits)
  }
  # Print extra parameters if present
  if(!is.null(x$extra)){
     cat("\n\nExtra parameter:\n")
     print(x$extra,digits=digits)
  }
}
#'
#'
#' @method summary mtarmcmc
#' @export
summary.mtarmcmc <- function(object,...,quantiles=c(0.025,0.25,0.5,0.75,0.975)){
  out <- list()
  # Store basic MCMC information
  out$start <- object$start
  out$end <- object$end
  out$thin <- object$thin
  out$sz <- object$sz
  out$regim <- object$regim
  # Number of requested quantiles
  q <- length(quantiles)
  # Helper function to compute summary statistics and quantiles
  myfunc <- function(x){
    # Extract mean, standard deviation, and Monte Carlo error, then append quantiles
    x2 <- c(summary(as.mcmc(x))$statistics[-3],quantile(x,probs=quantiles))
    names(x2) <- c("Mean","Sd","Sd(Mean)",paste0(round(quantiles*100,digits=1),"%"))
    return(x2)
  }
  # Summarize threshold parameters if present
  if(!is.null(object$thresholds)){
     out$thresholds <- t(apply(object$thresholds,2,myfunc))
     rownames(out$thresholds) <- colnames(object$thresholds)
  }
  # Initialize lists for regime-specific location and scale summaries
  out$location <- list()
  out$scale <- list()
  # Loop over regimes and summarize location and scale parameters
  for(j in 1:object$regim){
      out$location[[j]] <- t(apply(object$location[[j]],2,myfunc))
      rownames(out$location[[j]]) <- colnames(object$location[[j]])
      out$scale[[j]] <- t(apply(object$scale[[j]],2,myfunc))
      rownames(out$scale[[j]]) <- colnames(object$scale[[j]])
  }
  # Summarize skewness parameters if present
  if(!is.null(object$skewness)){
     out$skewness <- t(apply(object$skewness,2,myfunc))
     rownames(out$skewness) <- colnames(object$skewness)
  }
  # Summarize extra distributional parameters if present
  if(!is.null(object$extra)){
     out$extra <- t(apply(object$extra,2,myfunc))
     rownames(out$extra) <- colnames(object$extra)
  }
  # Assign class to the summary object and return it
  class(out) <- "summtarmcmc"
  return(out)
}
#'
#'
#' @method plot mtarmcmc
#' @export
plot.mtarmcmc <- function(x, trace=TRUE, density=TRUE, smooth=FALSE, bwf, auto.layout=TRUE, ask=dev.interactive(), ...){
  # Loop over regimes to plot regime-specific parameters
  for(j in 1:x$regim){
      # Open a new graphics device for location parameters
      dev.new()
      cat(paste("\nAutoregressive coefficients Regime",j,"\n"))
      plot(x$location[[j]], trace=trace, density=density, smooth=smooth, bwf,
           auto.layout=auto.layout, ask=ask, ...)
      # Open a new graphics device for scale parameters
      dev.new()
      cat(paste("\nScale parameter Regime",j,"\n"))
      plot(x$scale[[j]], trace=trace, density=density, smooth=smooth, bwf,
           auto.layout=auto.layout, ask=ask, ...)
  }
  # Plot skewness parameters if present
  if(!is.null(x$skewness)){
     dev.new()
     cat("\nSkewness parameter\n")
     plot(x$skewness, trace=trace, density=density, smooth=smooth, bwf,
          auto.layout=auto.layout, ask=ask, ...)
  }
  # Plot extra parameters if present
  if(!is.null(x$extra)){
     dev.new()
     cat("\nExtra parameter\n")
     plot(x$extra, trace=trace, density=density, smooth=smooth, bwf,
          auto.layout=auto.layout, ask=ask, ...)
  }
}
#'
#' @title Highest Posterior Density intervals for objects of class \code{mtar}
#' @param obj an object of class \code{mtar} generated by a call to the function \code{mtar()}.
#' @param prob a numeric scalar in the interval \eqn{(0,1)} giving the target probability content of
#' the intervals. By default, \code{prob} is set to \code{0.95}.
#' @param ... Optional additional arguments for methods. None are used at present.
#' @method HPDinterval mtar
#' @seealso \code{\link[coda]{HPDinterval}}
#' @export
HPDinterval.mtar <- function(obj, prob=0.95, ...){
  # Convert the mtar object to an mcmc object
  obj2 <- as.mcmc(obj)
  # Remove components not needed for HPD interval output
  obj2$regim <- obj2$sz <- obj2$end <- obj2$start <- obj2$thin <- obj2$delay <- NULL
  # Compute HPD intervals for thresholds if they are present
  if(!is.null(obj2$thresholds)){
     obj2$thresholds <- HPDinterval(obj2$thresholds,prob=prob, ...)
     attr(obj2$thresholds,"Probability") <- NULL
  }
  # Initialize temporary objects to merge location parameters across regimes
  temp <- data.frame(ns=NA)
  temp2 <- matrix(NA,ncol(obj2$scale[[1]]),1)
  rownames(temp2) <- colnames(obj2$scale[[1]])
  # Helper vectors used to reorder and relabel coefficients
  cc <- c(colnames(obj$data[[1]]$y),":Time",":Season")
  cc2 <- c(paste0(1:ncol(obj$data[[1]]$y),colnames(obj$data[[1]]$y)),":.1Time",":.2Season")
  # Loop over regimes to compute HPD intervals for location and scale parameters
  for(j in 1:obj$ars$nregim){
      # HPD intervals for location parameters
      obj2$location[[j]] <- HPDinterval(obj2$location[[j]],prob=prob, ...)
      obj2$location[[j]] <- data.frame(ns=rownames(obj2$location[[j]]),ps=obj2$location[[j]])
      # Merge parameters across regimes
      temp <- merge(temp,obj2$location[[j]],by.x="ns",by.y="ns",all.x=TRUE,all.y=TRUE)
      # HPD intervals for scale parameters
      obj2$scale[[j]] <- HPDinterval(obj2$scale[[j]],prob=prob, ...)
      temp2 <- cbind(temp2,obj2$scale[[j]])
  }
  # Remove placeholder rows
  temp <- temp[!is.na(temp[,1]),]
  for(jj in 1:(length(cc))) temp[,1] <- gsub(cc[jj],cc2[jj],temp[,1])
  temp <- temp[sort(temp[,1],index=TRUE)$ix,]
  # Temporary renaming to ensure correct ordering of location parameters
  for(jj in 1:(length(cc))) temp[,1] <- gsub(cc2[jj],cc[jj],temp[,1])
  # Store location HPD intervals as a matrix
  obj2$location <- as.matrix(temp[,-1])
  rownames(obj2$location) <- temp[,1]
  # Store scale HPD intervals as a matrix
  obj2$scale <- matrix(temp2[,-1],nrow(temp2),2*obj$ars$nregim)
  rownames(obj2$scale) <- rownames(temp2)
  # Set column names for lower and upper bounds by regime
  colnames(obj2$location) <- colnames(obj2$scale) <- rep(c("lower","upper"),obj$ars$nregim)
  colnames(obj2$location) <- colnames(obj2$scale) <- sort(c(paste0("Regime ",1:obj$ars$nregim,":lower"),paste0("Regime ",1:obj$ars$nregim,":upper")))
  # Compute HPD intervals for skewness parameters if present
  if(!is.null(obj2$skewness)){
     obj2$skewness <- HPDinterval(obj2$skewness,prob=prob, ...)
     attr(obj2$skewness,"Probability") <- NULL
  }
  # Compute HPD intervals for extra parameters if present
  if(!is.null(obj2$extra)){
     obj2$extra <- HPDinterval(obj2$extra,prob=prob, ...)
     attr(obj2$extra,"Probability") <- NULL
  }
  # Store the probability level used for the HPD intervals
  obj2$prob <- prob
  # Set class for the HPD interval object
  class(obj2) <- "HDPmtar"
  return(obj2)
}

#' @method print HDPmtar
#' @export
print.HDPmtar <- function(x, digits=max(3, getOption("digits") - 2), ...){
  # Display the probability level used for the HPD intervals
  cat("\nProbability = ",x$prob,"\n\n")
  # Print HPD intervals for threshold parameters if they are present
  if(!is.null(x$thresholds)){
     cat("Thresholds:\n")
     print(x$thresholds,digits=digits,na.print="")
  }
  # Print HPD intervals for location parameters
  cat("\n\nAutoregressive coefficients:\n")
  print(x$location,digits=digits,na.print="")
  # Print HPD intervals for scale parameters
  cat("\n\nScale parameter:\n")
  print(x$scale,digits=digits,na.print="")
  # Print HPD intervals for skewness parameters if they are present
  if(!is.null(x$skewness)){
     cat("\n\nSkewness parameter:\n")
     print(x$skewness,digits=digits,na.print="")
  }
  # Print HPD intervals for extra parameters if they are present
  if(!is.null(x$extra)){
     cat("\n\nExtra parameter:\n")
     print(x$extra,digits=digits,na.print="")
  }
}

#' @method DIC listmtar
#' @export
DIC.listmtar <- function(x,...){
  # Initialize an empty vector to store DIC values
  out <- vector()
  # Loop over each element in the list and compute its DIC
  for(i in 1:length(x)) out <- c(out,DIC(x[[i]]))
  # Convert the vector of DIC values into a column matrix
  out <- matrix(out,length(out),1)
  # Assign column and row names to the output matrix
  colnames(out) <- "DIC"
  rownames(out) <- names(x)
  # Return the matrix of DIC values
  return(out)
}

#' @method WAIC listmtar
#' @export
WAIC.listmtar <- function(x,...){
  # Initialize an empty vector to store WAIC values
  out <- vector()
  # Loop over each element in the list and compute its WAIC
  for(i in 1:length(x)) out <- c(out,WAIC(x[[i]]))
  # Convert the vector of WAIC values into a column matrix
  out <- matrix(out,length(out),1)
  # Assign column and row names to the output matrix
  colnames(out) <- "WAIC"
  rownames(out) <- names(x)
  # Return the matrix of WAIC values
  return(out)
}
#'
