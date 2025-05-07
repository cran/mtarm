#' @importFrom grDevices dev.new
#' @importFrom graphics abline hist lines par
#' @importFrom methods missingArg
#' @importFrom stats pnorm pt density median model.extract model.matrix pgamma qgamma qnorm quantile rbeta rbinom rexp rgamma rnorm runif sd terms na.omit acf pacf qqnorm residuals
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Formula Formula model.part
#' @importFrom mvtnorm pmvnorm pmvt
#' @importFrom GIGrvg rgig
#' @importFrom coda mcmc
#' @title Returns of the closing prices of three financial indexes
#'
#' @description These data correspond to the returns on closing prices of the Colcap,
#' Bovespa, and S&P 500 indexes from February 10, 2010 to March 31, 2016 (1505 time
#' points). Colcap is a leading indicator of the price dynamics of the 20 most liquid
#' shares on the Colombian stock market. Bovespa is the Brazilian stock market index,
#' the world's thirteenth largest and most important stock exchange, and the first in
#' Latin America. Finally, the Standard & Poor's 500 (S&P 500) index is a stock index
#' based on the 500 largest companies in the United States.
#'
#' @docType data
#'
#' @usage data(returns)
#'
#' @format A data frame with 1505 rows and 4 variables:
#' \describe{
#'   \item{Date}{a vector that indicates the date each measurement was performed.}
#'   \item{COLCAP}{a numerical vector indicating the returns on closing prices of COLCAP.}
#'   \item{SP500}{a numerical vector indicating the returns on closing prices of SP500.}
#'   \item{BOVESPA}{a numerical vector indicating the returns on closing prices of BOVESPA.}
#' }
#' @keywords datasets
#' @references Romero, L.V. and Calderon, S.A. (2021) Bayesian estimation of a multivariate TAR model when the noise
#'             process follows a Student-t distribution. Communications in Statistics - Theory and Methods, 50, 2508-2530.
#' @examples
#' data(returns)
#' dev.new()
#' plot(ts(as.matrix(returns[,-1])), main="Returns")
#'
"returns"
#'
#' @title Rainfall and two river flows in Colombia
#'
#' @description The data represent daily rainfall (in mm) and two river flows (in \eqn{m^3}/s)
#' in southern Colombia. A meteorological station located at an altitude of 2400 meters was
#' used to measure rainfall. The El Trebol hydrological station was used to measure the flow
#' in the Bedon river at an altitude of 1720 meters. The Villalosada hydrological station
#' measured the flow in the La Plata river at an altitude of 1300 meters. Geographically, the
#' stations are located near the equator. The last characteristic allows for control over
#' hydrological and meteorological factors that might distort the dynamic relationship between
#' rainfall and river flows. January 1, 2006, to April 14, 2009, was the sample period.
#' @docType data
#'
#' @usage data(riverflows)
#'
#' @format A data frame with 1200 rows and 4 variables:
#' \describe{
#'   \item{Date}{a vector that indicates the date each measurement was performed.}
#'   \item{Bedon}{a numerical vector indicating the Bedon river flow.}
#'   \item{LaPlata}{a numerical vector indicating the La Plata river flow.}
#'   \item{Rainfall}{a numerical vector indicating the rainfall.}
#' }
#' @keywords datasets
#' @references Calderon, S.A. and Nieto, F.H. (2017) Bayesian analysis of multivariate threshold autoregressive models
#'             with missing data. Communications in Statistics - Theory and Methods, 46, 296-318.
#' @examples
#' data(riverflows)
#' dev.new()
#' plot(ts(as.matrix(riverflows[,-1])), main="Rainfall and river flows")
#'
"riverflows"
#'
#'
#' @title Deviance information criterion (DIC)
#' @description This function computes the Deviance Information Criterion (DIC) for objects of class \code{mtar}.
#' @param ...	one or several objects of the class \emph{mtar}.
#' @param verbose an (optional) logical switch indicating if should the report of results be printed. By default,\code{verbose} is set to TRUE.
#' @param digits an (optional) integer indicating the number of digits to print. By default,\code{digits} is set to \code{max(3, getOption("digits") - 2)}.
#' @return A \code{data.frame} with the values of the DIC for each \emph{mtar} object in the input.
#' @references Spiegelhalter D.J., Best N.G., Carlin B.P. and Van Der Linde A. (2002) Bayesian Measures of Model Complexity and Fit.
#'             Journal of the Royal Statistical Society Series B (Statistical Methodology), 64(4), 583–639.
#' @references Spiegelhalter D.J., Best N.G., Carlin B.P. and Van der Linde A. (2014). The deviance information criterion:
#'             12 years on. Journal of the Royal Statistical Society Series B (Statistical Methodology), 76(3), 485–493.
#' @export DIC
#' @seealso \link{WAIC}
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1a <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'               dist="Gaussian", ars=list(p=c(1,1,2)), n.burnin=2000,
#'               n.sim=3000, n.thin=2)
#' fit1b <- update(fit1a,dist="Slash")
#' fit1c <- update(fit1a,dist="Student-t")
#' DIC(fit1a,fit1b,fit1c)
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2a <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'               dist="Gaussian", ars=list(p=c(5,5,5)), n.burnin=2000,
#'               n.sim=3000, n.thin=2)
#' fit2b <- update(fit2a,dist="Slash")
#' fit2c <- update(fit2a,dist="Student-t")
#' DIC(fit2a,fit2b,fit2c)
#' }
#'
DIC <- function(...,verbose=TRUE,digits=max(3, getOption("digits") - 2)){
  mtlogLik <- function(dist,y,X,beta,Sigma,delta,log,nu){
    resu <- y-X%*%beta
    if(dist %in% c("Skew-normal","Skew-Student-t")){
      A <- chol2inv(chol(Sigma + as.vector(delta^2)*diag(k)))
      muv <- resu%*%(A*matrix(delta,k,k,byrow=TRUE))
      Sigmav <- diag(k) - matrix(delta,k,k)*A*matrix(delta,k,k,byrow=TRUE)
      out <- colSums(t(resu)*tcrossprod(A,resu))
      if(dist=="Skew-normal"){
        if(k > 1) sum0 <- apply(muv,1,function(x) pmvnorm(lower=-x,upper=rep(Inf,k),sigma=Sigmav))
        else sum0 <- apply(muv,1,function(x) pnorm(-x/sqrt(Sigmav),lower.tail=FALSE))
        out <- 2^k*exp(-out/2)*(det(A)^(1/2))*sum0/(2*pi)^(k/2)
      }
      if(dist=="Skew-Student-t"){
        muv <- muv*matrix(sqrt((nu + k)/(nu + out)),nrow(muv),k)
        if(k > 1) sum0 <- apply(muv,1,function(x) pmvt(lower=-x,upper=rep(Inf,k),sigma=Sigmav,df=round(nu,0)+k))
        else sum0 <- apply(muv,1,function(x) pt(-x/sqrt(Sigmav),df=nu,lower.tail=FALSE))
        out <- 2^k*(1 + out/nu)^(-(nu+k)/2)*(det(A)^(1/2))*gamma((nu+k)/2)*sum0/((nu*pi)^(k/2)*gamma(nu/2))
      }
    }else{
      out <- colSums(t(resu)*tcrossprod(chol2inv(chol(Sigma)),resu))
      if(dist=="Gaussian") out <- exp(-out/2)
      if(dist=="Laplace")
        out <- besselK(sqrt(out)/2,(2-k)/2)*out^((2-k)/4)/(2^((k+2)/2))
      if(dist=="Student-t")
        out <- (1 + out/nu)^(-(nu+k)/2)*gamma((nu+k)/2)/((nu/2)^(k/2)*gamma(nu/2))
      if(dist=="Contaminated normal")
        out <- nu[1]*exp(-out*nu[2]/2)*nu[2]^(k/2) + (1-nu[1])*exp(-out/2)
      if(dist=="Slash")
        out <- ifelse(out==0,nu/(nu+k),gamma((k+nu)/2)*(nu/2)*pgamma(1,shape=(k+nu)/2,rate=out/2)/((out/2)^((k+nu)/2)))
      if(dist=="Hyperbolic")
        out <- besselK(nu*sqrt(1+out),(2-k)/2)*(1+out)^((2-k)/4)*nu^(k/2)/besselK(nu,1)
      out <- out/((2*pi)^(k/2)*det(Sigma)^(1/2))
    }
    out <- -2*sum(log(out))
    if(log) out <- out + 2*sum(y)
    return(out)
  }
  another <- list(...)
  if(any(lapply(another,function(xx) class(xx)[1])!="mtar"))
    stop("Only mtar-type objects are supported!!",call.=FALSE)
  call. <- match.call()
  out <- matrix(0,length(another),1)
  outnames <- vector()

  for(l in 1:(length(another))){
    n.sim <- another[[l]]$n.sim
    dist <- another[[l]]$dist
    Dbar <- vector()
    k <- ncol(another[[l]]$data[[1]]$y)
    Dbarv <- vector()
    if(another[[l]]$regim > 1) lims <- (another[[l]]$ps+1):(length(another[[l]]$threshold.series))
    for(j in 1:n.sim){
      Dbar <- 0
      if(another[[l]]$regim > 1){
        Z <- another[[l]]$threshold.series[lims-another[[l]]$chains$h[j]]
        regs <- cut(Z,breaks=c(-Inf,sort(another[[l]]$chains$thresholds[,j]),Inf),labels=FALSE)
      }else regs <- matrix(1,nrow(another[[l]]$data[[1]]$y),1)
      for(i in 1:another[[l]]$regim){
        betai <- matrix(another[[l]]$chains[[i]]$location[,((j-1)*k + 1):(j*k)],ncol(another[[l]]$data[[i]]$X),k)
        Sigmai <- matrix(another[[l]]$chains[[i]]$scale[,((j-1)*k + 1):(j*k)],k,k)
        places <- regs==i
        if(another[[l]]$dist %in% c("Skew-Student-t"))
          Dbar <- Dbar + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])
        if(another[[l]]$dist %in% c("Skew-normal"))
          Dbar <- Dbar + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log)
        if(another[[l]]$dist %in% c("Student-t","Hyperbolic","Slash","Contaminated normal"))
          Dbar <- Dbar + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])
        if(another[[l]]$dist %in% c("Gaussian","Laplace"))
          Dbar <- Dbar + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,log=another[[l]]$log)
      }
      Dbarv <- c(Dbarv,Dbar)
    }
    Dhat <- 0
    if(another[[l]]$dist %in% c("Skew-Student-t","Student-t","Hyperbolic","Slash","Contaminated normal"))
      extra <- rowMeans(matrix(another[[l]]$chains$extra,ncol=n.sim))
    if(another[[l]]$dist %in% c("Skew-Student-t","Skew-normal"))
      delta <- rowMeans(matrix(another[[l]]$chains$delta,ncol=n.sim))
    if(another[[l]]$regim > 1){
      h <- as.integer(round(mean(another[[l]]$chains$h)))
      thresholds <- rowMeans(matrix(another[[l]]$chains$thresholds,ncol=n.sim))
      Z <- another[[l]]$threshold.series[lims-h]
      regs <- cut(Z,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
    }else regs <- matrix(1,nrow(another[[l]]$data[[1]]$y),1)
    for(i in 1:another[[l]]$regim){
      location <- matrix(0,nrow(another[[l]]$chains[[i]]$location),k)
      scale <- matrix(0,k,k)
      for(j in 1:k){
        location[,j] <- rowMeans(matrix(another[[l]]$chains[[i]]$location[,seq(j,n.sim*k,by=k)],nrow(location),n.sim))
        scale[,j] <- rowMeans(matrix(another[[l]]$chains[[i]]$scale[,seq(j,n.sim*k,by=k)],k,n.sim))
      }
      places <- regs==i
      if(another[[l]]$dist %in% c("Skew-Student-t"))
        Dhat <- Dhat + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=location,Sigma=scale,delta=delta,log=another[[l]]$log,nu=extra)
      if(another[[l]]$dist %in% c("Skew-normal"))
        Dhat <- Dhat + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=location,Sigma=scale,delta=delta,log=another[[l]]$log)
      if(another[[l]]$dist %in% c("Student-t","Hyperbolic","Slash","Contaminated normal"))
        Dhat <- Dhat + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=location,Sigma=scale,log=another[[l]]$log,nu=extra)
      if(another[[l]]$dist %in% c("Gaussian","Laplace"))
        Dhat <- Dhat + mtlogLik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=location,Sigma=scale,log=another[[l]]$log)
    }
    out[l] <- Dhat + 2*(mean(Dbarv) - Dhat)
    outnames[l] <- as.character(call.[l+1])
  }
  rownames(out) <- outnames
  colnames(out) <- "DIC"
  out <- round(out,digits=digits)
  if(verbose) print(out)
  return(invisible(out))
}
#'
#' @title Watanabe-Akaike or Widely Available Information Criterion (WAIC)
#' @description This function computes the Watanabe-Akaike or Widely Available Information Criterion (WAIC) for objects of class \code{mtar}.
#' @param ...	one or several objects of the class \emph{mtar}.
#' @param verbose an (optional) logical switch indicating if should the report of results be printed. By default,\code{verbose} is set to TRUE.
#' @param digits an (optional) integer indicating the number of digits to print. By default,\code{digits} is set to \code{max(3, getOption("digits") - 2)}.
#' @return A \code{data.frame} with the values of the WAIC for each \emph{mtar} object in the input.
#' @references Watanabe S. (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in
#'             Singular Learning Theory. The Journal of Machine Learning Research, 11, 3571–3594.
#' @export WAIC
#' @seealso \link{DIC}
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1a <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'               dist="Gaussian", ars=list(p=c(1,1,2)), n.burnin=100,
#'               n.sim=3000, n.thin=2)
#' fit1b <- update(fit1a,dist="Slash")
#' fit1c <- update(fit1a,dist="Student-t")
#' WAIC(fit1a,fit1b,fit1c)
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2a <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'               dist="Gaussian", ars=list(p=c(5,5,5)), n.burnin=100,
#'               n.sim=3000, n.thin=2)
#' fit2b <- update(fit2a,dist="Slash")
#' fit2c <- update(fit2a,dist="Student-t")
#' WAIC(fit2a,fit2b,fit2c)
#' }
#'
WAIC <- function(...,verbose=TRUE,digits=max(3, getOption("digits") - 2)){
  Lik <- function(dist,y,X,beta,Sigma,delta,log,nu){
    resu <- y-X%*%beta
    if(dist %in% c("Skew-normal","Skew-Student-t")){
      A <- chol2inv(chol(Sigma + as.vector(delta^2)*diag(k)))
      muv <- resu%*%(A*matrix(delta,k,k,byrow=TRUE))
      Sigmav <- diag(k) - matrix(delta,k,k)*A*matrix(delta,k,k,byrow=TRUE)
      out <- colSums(t(resu)*tcrossprod(A,resu))
      if(dist=="Skew-normal"){
        if(k > 1) sum0 <- apply(muv,1,function(x) pmvnorm(lower=-x,upper=rep(Inf,k),sigma=Sigmav))
        else sum0 <- apply(muv,1,function(x) pnorm(-x/sqrt(Sigmav),lower.tail=FALSE))
        out <- 2^k*exp(-out/2)*(det(A)^(1/2))*sum0/(2*pi)^(k/2)
      }
      if(dist=="Skew-Student-t"){
        muv <- muv*matrix(sqrt((nu + k)/(nu + out)),nrow(muv),k)
        if(k > 1) sum0 <- apply(muv,1,function(x) pmvt(lower=-x,upper=rep(Inf,k),sigma=Sigmav,df=round(nu,0)+k))
        else sum0 <- apply(muv,1,function(x) pt(-x/sqrt(Sigmav),df=nu,lower.tail=FALSE))
        out <- 2^k*(1 + out/nu)^(-(nu+k)/2)*(det(A)^(1/2))*gamma((nu+k)/2)*sum0/((nu*pi)^(k/2)*gamma(nu/2))
      }
    }else{
      out <- colSums(t(resu)*tcrossprod(chol2inv(chol(Sigma)),resu))
      if(dist=="Gaussian") out <- exp(-out/2)
      if(dist=="Laplace")
        out <- besselK(sqrt(out)/2,(2-k)/2)*out^((2-k)/4)/(2^((k+2)/2))
      if(dist=="Student-t")
        out <- (1 + out/nu)^(-(nu+k)/2)*gamma((nu+k)/2)/((nu/2)^(k/2)*gamma(nu/2))
      if(dist=="Contaminated normal")
        out <- nu[1]*exp(-out*nu[2]/2)*nu[2]^(k/2) + (1-nu[1])*exp(-out/2)
      if(dist=="Slash")
        out <- ifelse(out==0,nu/(nu+k),gamma((k+nu)/2)*(nu/2)*pgamma(1,shape=(k+nu)/2,rate=out/2)/((out/2)^((k+nu)/2)))
      if(dist=="Hyperbolic")
        out <- besselK(nu*sqrt(1+out),(2-k)/2)*(1+out)^((2-k)/4)*nu^(k/2)/besselK(nu,1)
      out <- out/((2*pi)^(k/2)*det(Sigma)^(1/2))
    }
    if(log) out <- out*apply(matrix(y,nrow(X),k),1,function(x) prod(exp(-x)))
    return(out)
  }
  another <- list(...)
  if(any(lapply(another,function(xx) class(xx)[1])!="mtar"))
    stop("Only mtar-type objects are supported!!",call.=FALSE)
  call. <- match.call()
  out <- matrix(0,length(another),1)
  outnames <- vector()

  for(l in 1:(length(another))){
    n.sim <- another[[l]]$n.sim
    dist <- another[[l]]$dist
    Dbar <- matrix(0,nrow(another[[l]]$data[[1]]$y),1)
    Dbarlog <- matrix(0,nrow(another[[l]]$data[[1]]$y),1)
    k <- ncol(another[[l]]$data[[1]]$y)
    if(another[[l]]$regim > 1) lims <- (another[[l]]$ps+1):(length(another[[l]]$threshold.series))
    for(j in 1:n.sim){
      for(i in 1:another[[l]]$regim){
        betai <- matrix(another[[l]]$chains[[i]]$location[,((j-1)*k + 1):(j*k)],ncol(another[[l]]$data[[i]]$X),k)
        Sigmai <- matrix(another[[l]]$chains[[i]]$scale[,((j-1)*k + 1):(j*k)],k,k)
        if(another[[l]]$regim > 1){
          Z <- another[[l]]$threshold.series[lims-another[[l]]$chains$h[j]]
          regs <- cut(Z,breaks=c(-Inf,sort(another[[l]]$chains$thresholds[,j]),Inf),labels=FALSE)
        }else regs <- matrix(1,nrow(another[[l]]$data[[i]]$y),1)
        places <- regs==i
        if(dist=="Skew-Student-t")
          tempi <- Lik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])
        if(dist=="Skew-normal")
          tempi <- Lik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,delta=another[[l]]$chains$delta[,j],log=another[[l]]$log)
        if(dist %in% c("Student-t","Hyperbolic","Slash","Contaminated normal"))
          tempi <- Lik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,log=another[[l]]$log,nu=another[[l]]$chains$extra[,j])
        if(dist %in% c("Gaussian","Laplace"))
          tempi <- Lik(dist=dist,y=another[[l]]$data[[i]]$y[places,],X=another[[l]]$data[[i]]$X[places,],beta=betai,Sigma=Sigmai,log=another[[l]]$log)
        Dbar[places] <- Dbar[places] + tempi/n.sim
        Dbarlog[places] <- Dbarlog[places] + log(tempi)/n.sim
      }
    }
    a <- sum(unlist(lapply(Dbar,function(x) sum(log(x)))))
    b <- sum(unlist(lapply(Dbarlog,sum)))
    out[l] <- -2*a + 4*(a-b)
    outnames[l] <- as.character(call.[l+1])
  }
  rownames(out) <- outnames
  colnames(out) <- "WAIC"
  out <- round(out,digits=digits)
  if(verbose) print(out)
  return(invisible(out))
}
#'
#' @title Bayesian estimation of a multivariate Threshold Autoregressive (TAR) model.
#' @description This function uses Gibbs sampling to generate a sample from the posterior
#'              distribution of the parameters of a multivariate TAR model when the noise
#'              process follows Gaussian, Student-\eqn{t}, Slash, Symmetric Hyperbolic,
#'              Contaminated normal, Laplace, Skew-normal or skew-\eqn{t} distribution.
#' @param formula A three-part expression of type \code{Formula} describing the TAR model
#'                to be fitted to the data. In the first part, the variables in the
#'                multivariate output series are listed; in the second part, the threshold
#'                series is specified, and in the third part, the variables in the
#'                multivariate exogenous series are specified.
#' @param ars a list composed of three objects, namely: \code{p}, \code{q} and \code{d},
#'            each of which corresponds to a vector of non-negative integers with as many
#'            elements as there are regimes in the TAR model.
#' @param Intercept an (optional) logical variable. If \code{TRUE}, then the model
#'                  includes an intercept.
#' @param data an (optional) data frame, list or environment (or object coercible by
#'             \link{as.data.frame} to a data frame) containing the variables in the model.
#'             If not found in data, the variables are taken from \code{environment(formula)},
#'             typically the environment from which \code{mtar} is called.
#' @param subset an (optional) vector specifying a subset of observations to be used in the
#'               fitting process.
#' @param dist an (optional) character string that allows the user to specify the
#'             multivariate distribution to be used to describe the noise process
#'             behavior. The available options are: Gaussian ("Gaussian"), Student-\eqn{t}
#'             ("Student-t"), Slash ("Slash"), Symmetric Hyperbolic ("Hyperbolic"),
#'             Laplace ("Laplace"), Contaminated normal ("Contaminated normal"),
#'             Skew-normal ("Skew-normal") and Skew-Student-\eqn{t} ("Skew-Student-t").
#'             By default, \code{dist} is set to "Gaussian".
#'
#' @param n.sim an (optional) positive integer specifying the required number of iterations
#'              for the simulation after the burn-in period. By default, \code{n.sim} is set
#'              to 500.
#' @param n.burnin an (optional) positive integer specifying the required number of burn-in
#'                 iterations for the simulation. By default, \code{n.burnin} is set to 100.
#' @param n.thin an (optional) positive integer specifying the required thinning interval
#'               for the simulation. By default, \code{n.thin} is set to 1.
#' @param row.names an (optional) vector that allows the user to name the time point to
#'                  which each row in the data set corresponds.
#' @param prior an (optional) list that allows the user to specify the values of the
#'              hyperparameters, that is, allows to specify the values of the parameters
#'              of the prior distributions.
#' @param log an (optional) logical variable. If \code{TRUE}, then the behaviour of the output
#'            series is described using the exponentiated version of \code{dist}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return an object of class \emph{mtar} in which the main results of the model fitted to the data are stored, i.e., a
#' list with components including
#' \tabular{ll}{
#' \code{chains}   \tab list with several arrays, which store the values of each model parameter in each iteration of the simulation,\cr
#' \tab \cr
#' \code{n.sim}    \tab number of iterations of the simulation after the burn-in period,\cr
#' \tab \cr
#' \code{n.burnin} \tab number of burn-in iterations in the simulation,\cr
#' \tab \cr
#' \code{n.thin}   \tab thinning interval in the simulation,\cr
#' \tab \cr
#' \code{regim}    \tab number of regimes, \cr
#' \tab \cr
#' \code{ars}      \tab list composed of three objects, namely: \code{p}, \code{q} and \code{d},
#'                      each of which corresponds to a vector of non-negative integers with as
#'                      many elements as there are regimes in the TAR model,\cr
#' \tab \cr
#' \code{dist}     \tab name of the multivariate distribution used to describe the behavior of
#'                      the noise process,\cr
#' \tab \cr
#' \code{threshold.series}  \tab vector with the values of the threshold series,\cr
#' \tab \cr
#' \code{response.series}   \tab matrix with the values of the output series,\cr
#' \tab \cr
#' \code{covariable.series} \tab matrix with the values of the exogenous series,\cr
#' \tab \cr
#' \code{Intercept}    \tab If \code{TRUE}, then the model included an intercept term,\cr
#' \tab \cr
#' \code{formula}      \tab the formula,\cr
#' \tab \cr
#' \code{call}         \tab the original function call.\cr
#' }
#' @export mtar
#' @seealso \link{DIC}, \link{WAIC}
#' @references Nieto, F.H. (2005) Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data.
#'             Communications in Statistics - Theory and Methods, 34, 905-930.
#' @references Romero, L.V. and Calderon, S.A. (2021) Bayesian estimation of a multivariate TAR model when the noise
#'             process follows a Student-t distribution. Communications in Statistics - Theory and Methods, 50, 2508-2530.
#' @references Calderon, S.A. and Nieto, F.H. (2017) Bayesian analysis of multivariate threshold autoregressive models
#'             with missing data. Communications in Statistics - Theory and Methods, 46, 296-318.
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1 <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'              dist="Gaussian", ars=list(p=c(1,1,2)), n.burnin=100,
#'              n.sim=3000, n.thin=2)
#' summary(fit1)
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2 <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'              dist="Gaussian", ars=list(p=c(5,5,5)), n.burnin=2000,
#'              n.sim=3000, n.thin=2)
#' summary(fit2)
#' }
#'
mtar <- function(formula, data, subset, Intercept=TRUE, ars, row.names, dist="Gaussian", prior=list(), n.sim=500, n.burnin=100, n.thin=1, log=FALSE, ...){
  if(!(dist %in% c("Gaussian","Student-t","Hyperbolic","Laplace","Slash","Contaminated normal","Skew-Student-t","Skew-normal")))
    stop("Only 'Gaussian', 'Student-t', 'Hyperbolic', 'Laplace', 'Slash', 'Contaminated normal', 'Skew-normal' and 'Skew-Student-t' distributions are supported!",call.=FALSE)
  if(missing(data)) data <- environment(formula)
  ars$p <- ceiling(abs(ars$p))
  regim <- length(ars$p)
  mmf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "row.names"), names(mmf), 0)
  mmf <- mmf[c(1,m)]
  mmf$drop.unused.levels <- TRUE
  mmf[[1]] <- as.name("model.frame")
  mmf$formula <- Formula(formula)
  mmf <- eval(mmf, parent.frame())
  if(!missingArg(row.names)) row.names <- as.vector(as.character(model.extract(mmf,row.names)))

  mx <- model.part(Formula(formula), data = mmf, rhs = 1, terms = TRUE)
  D <- model.matrix(mx, data = mmf)
  if(attr(terms(mx),"intercept")){
    Dnames <- colnames(D)
    D <- matrix(D[,-1],ncol=length(Dnames)-1)
    colnames(D) <- Dnames[-1]
  }
  if(log){
    if(any(D<0)) stop(paste0("There are non-positive values in the output series, so it cannot be described using the log-",tolower(dist)," distribution."),call.=FALSE)
    D <- log(D)
  }
  k <- ncol(D)
  if(regim > 1){
    mz <- model.part(Formula(formula), data = mmf, rhs = 2, terms = TRUE)
    Z <- model.matrix(mz, data = mmf)
    if(attr(terms(mz),"intercept")){
      Znames <- colnames(Z)
      Z <- as.matrix(Z[,-1])
      colnames(Z) <- Znames[-1]
    }
    ta <- vector()
  }
  if(!is.null(ars$q)){
    ars$q <- ceiling(abs(ars$q))
    if(max(ars$q) > 0){
      mx2 <- model.part(Formula(formula), data = mmf, rhs = 3, terms = TRUE)
      X2 <- model.matrix(mx2, data = mmf)
      if(attr(terms(mx2),"intercept")){
        X2names <- colnames(X2)
        X2 <- as.matrix(X2[,-1])
        colnames(X2) <- X2names[-1]
      }
      r <- ncol(X2)
    }
  }else ars$q <- rep(0,regim)
  if(is.null(ars$d)) ars$d <- rep(0,regim) else ars$d <- ceiling(abs(ars$d))
  if(is.null(prior$theta0)) prior$theta0 <- 0
  if(is.null(prior$delta0)) prior$delta0 <- 1000000000
  if(is.null(prior$omega0)) prior$omega0 <- 1/1000000000
  if(is.null(prior$phi0)) prior$phi0 <- 1000000000
  Omega0 <- diag(k)*prior$omega0
  Phi0 <- diag(k)*prior$phi0
  if(is.null(prior$tau0)) prior$tau0 <- k
  if(is.null(prior$hmin)) prior$hmin <- 0 else prior$hmin <- ceiling(abs(prior$hmin))
  if(is.null(prior$hmax)) prior$hmax <- 3 else prior$hmax <- ceiling(abs(prior$hmax))
  prior$hmin <- min(prior$hmin,prior$hmax)
  prior$hmax <- max(prior$hmin,prior$hmax)
  data <- list()
  chains <- list()
  n.sim <- ceiling(abs(n.sim))
  n.thin <- ceiling(abs(n.thin))
  n.burnin <- ceiling(abs(n.burnin))
  rep <- n.sim*n.thin + n.burnin
  ids <- matrix(seq(1,n.sim*n.thin*k,k))
  ids <- ids[seq(1,n.sim*n.thin,n.thin)]
  ids <- n.burnin*k + as.vector(apply(matrix(ids),1,function(x) x + c(0:(k-1))))
  Sigmanew2 <- list()
  name <- list()
  tn1 <- function(mu,var){
    temp <- pnorm(0,mean=mu,sd=sqrt(var))
    temp <- temp + runif(1)*(1 - temp)
    return(ifelse(temp>0.9999999999999999,-mu*0.001,qnorm(temp)*sqrt(var) + mu))
  }
  myf1 <- function(x){
    mu <- ais[,x]
    var <- Ais2/usi[x]
    tn1(mu=mu,var=var)
  }
  myf2 <- function(x){
    mu <- ais[,x]
    var <- Ais2/usi[x]
    ind <- TRUE
    indc <- 1
    while(ind & indc<=30){
      temp <- crossprod(Ais3/sqrt(usi[x]),matrix(rnorm(k*50),k,50))
      temp2 <- colSums(mu + temp > 0)==k
      indc <- indc + 1
      ind <- !any(temp2)
    }
    if(ind){
      myfout <- vector()
      for(indc in 1:k) myfout <- c(myfout,tn1(mu=mu[indc],var=var[indc,indc]))
    }else myfout <- (mu + temp)[,min(c(1:50)[temp2])]
    myfout
  }
  if(regim > 1) ps <- max(ars$p,ars$q,ars$d,prior$hmax) else ps <- max(ars$p,ars$q,ars$d)
  if(regim > 1){
    Zs <- Z[(ps+1-prior$hmin):(nrow(D)-prior$hmin),]
    t1 <- max(Zs);t0 <- min(Zs)
    thresholds <- quantile(Zs,probs=c(1:(regim-1))/regim)
    regs <- cut(Zs,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
    thresholds.chains <- matrix(thresholds,regim-1,1)
    hs.chains <- matrix(prior$hmin,1,1)
  }else regs <- matrix(1,nrow(D)-ps,1)
  for(i in 1:regim){
    y <- matrix(D[(ps+1):nrow(D),1:k],ncol=k)
    X <- matrix(1,nrow(y),1)
    for(j in 1:ars$p[i]) X <- cbind(X,D[((ps+1)-j):(nrow(D)-j),])
    name[[i]] <-  c("(Intercept)",paste0(rep(colnames(D)[1:k],ars$p[i]),sort(paste0(".lag(",rep(1:ars$p[i],k))),")"))
    if(!Intercept){
      X <- matrix(X[,-1],nrow(X),ncol(X)-1)
      name[[i]] <- name[[i]][-1]
    }
    if(ars$q[i]!=0){
      for(j in 1:ars$q[i]) X <- cbind(X,X2[((ps+1)-j):(nrow(D)-j),])
      name[[i]] <- c(name[[i]],paste0(rep(colnames(X2),ars$q[i]),sort(paste0(".lag(",rep(1:ars$q[i],r))),")"))
    }
    if(ars$d[i]!=0){
      for(j in 1:ars$d[i]) X <- cbind(X,Z[((ps+1)-j):(nrow(D)-j),])
      name[[i]] <- c(name[[i]],paste0(rep(colnames(Z),ars$d[i]),sort(paste0(".lag(",1:ars$d[i])),")"))
    }
    colnames(X) <- name[[i]]
    if(!missingArg(row.names)){
      row.names2 <- row.names[(ps+1):nrow(D)]
      rownames(X) <- rownames(y) <- row.names2
    }
    colnames(y) <- colnames(D)
    data[[i]] <- list()
    data[[i]]$y <- y
    data[[i]]$X <- X
    places <- regs == i
    X <- matrix(X[places,],sum(places),ncol(X))
    y <- matrix(y[places,],nrow(X),ncol(y))
    b0 <- chol2inv(chol(crossprod(X)))
    betanew <- crossprod(b0,crossprod(X,y))
    Sigmanew <- crossprod(y-X%*%betanew)/nrow(X)
    Sigmanew2[[i]] <- chol2inv(chol(Sigmanew))
    chains[[i]] <- list()
    chains[[i]]$location <- betanew
    chains[[i]]$scale <- Sigmanew
  }
  tol <- unlist(lapply(data,function(x) ncol(x$X)))*k*1.2
  us <- matrix(1,nrow(data[[1]]$X),1)
  ss <- matrix(0,nrow(data[[1]]$X),1)
  zs <- matrix(0,nrow(data[[1]]$y),k)
  chains$delta <- matrix(0,k,1)
  if(dist=="Hyperbolic"){
    if(is.null(prior$gamma0)) prior$gamma0 <- 0.1
    if(is.null(prior$eta0)) prior$eta0 <- 4
    chains$extra <- matrix((prior$gamma0+prior$eta0)/2,1,1)
    nus <- matrix(seq(prior$gamma0,prior$eta0,length=10001),10001,1)
    num1 <- nus[2:length(nus)] - nus[1:(length(nus)-1)]
    num2 <- (nus[2:length(nus)] + nus[1:(length(nus)-1)])/2
    resto <- log(nus) - log(besselK(nus,nu=1))
  }
  if(dist %in% c("Student-t","Skew-Student-t")){
    if(is.null(prior$gamma0)) prior$gamma0 <- 2
    if(is.null(prior$eta0)) prior$eta0 <- 100
    chains$extra <- matrix((prior$gamma0+prior$eta0)/2,1,1)
    nus <- matrix(seq(prior$gamma0,prior$eta0,length=10001),10001,1)
    num1 <- nus[2:length(nus)] - nus[1:(length(nus)-1)]
    num2 <- (nus[2:length(nus)] + nus[1:(length(nus)-1)])/2
    resto <- (nus/2)*log(nus/2) - lgamma(nus/2)
  }
  if(dist=="Slash"){
    if(is.null(prior$gamma0)) prior$gamma0 <- 1/1000000000
    if(is.null(prior$eta0)) prior$eta0 <- 1/1000000000
    chains$extra <- matrix(100,1,1)
  }
  if(dist=="Contaminated normal"){
    chains$extra <- matrix(c(0.01,0.99),2,1)
    if(is.null(prior$gamma01)) prior$gamma01 <- 1/1000000000
    if(is.null(prior$eta01)) prior$eta01 <- 1/1000000000
    if(is.null(prior$gamma02)) prior$gamma02 <- 1/1000000000
    if(is.null(prior$eta02)) prior$eta02 <- 1/1000000000
  }
  Loglik <- function(h,thresholds){
    regs <- cut(Z[(ps+1-h):(nrow(D)-h),],breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
    result <- 0
    for(i in 1:regim){
      places <- regs == i
      X <- matrix(data[[i]]$X[places,],sum(places),ncol(data[[i]]$X))
      y <- matrix(data[[i]]$y[places,] - matrix(chains$delta[,ncol(chains$delta)],nrow(X),k,byrow=TRUE)*zs[places,],nrow(X),k)
      usi <- us[places]
      Xu <- matrix(sqrt(usi),nrow(X),ncol(X))*X
      yu <- matrix(sqrt(usi),nrow(X),k)*y
      resu <- yu-Xu%*%chains[[i]]$location[,(ncol(chains[[i]]$location)-k+1):ncol(chains[[i]]$location)]
      ds <- colSums(t(resu)*tcrossprod(Sigmanew2[[i]],resu))
      result <- result -0.5*sum(ds - log(det(Sigmanew2[[i]])) - k*log(usi) + k*log(2*pi))
    }
    return(result)
  }
  bar <- txtProgressBar(min=0, max=rep, initial=0, width=min(50,rep), char="+", style=3)
  for(j in 1:rep){
    if(dist %in% c("Skew-normal","Skew-Student-t")){
      bis <- matrix(0,k,1)
      Bis <- matrix(0,k,k)
    }
    if(regim > 1){
      hs <- hs.chains[,j]
      thresholds <- thresholds.chains[,j]
      regs <- cut(Z[(ps+1-hs):(nrow(D)-hs),],breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
    }
    for(i in 1:regim){
      places <- regs == i
      X <- matrix(data[[i]]$X[places,],sum(places),ncol(data[[i]]$X))
      n <- nrow(X); s <- ncol(X)
      y2z <- matrix(data[[i]]$y[places,],n,k)
      y <- matrix(y2z - matrix(chains$delta[,j],n,k,byrow=TRUE)*zs[places,],n,k)
      mu0s <- matrix(prior$theta0,s,k)
      Sigmarinv <- diag(s)/prior$delta0
      if(dist %in% c("Gaussian","Skew-normal")) us[places] <- rep(1,n)
      else{
        resu <- y-X%*%chains[[i]]$location[,(1+(j-1)*k):(j*k)]
        usi <- colSums(t(resu)*tcrossprod(Sigmanew2[[i]],resu))
      }
      if(dist=="Laplace")
        us[places] <- apply(matrix(usi,length(usi),1),1,function(x) 1/rgig(n=1,lambda=(2-k)/2,chi=x,psi=1/4))
      if(dist=="Hyperbolic")
        us[places] <- apply(matrix(usi,length(usi),1),1,function(x) 1/rgig(n=1,lambda=(2-k)/2,chi=x+1,psi=chains$extra[,j]^2))
      if(dist=="Student-t")
        us[places] <- rgamma(n,shape=(chains$extra[,j]+k)/2,scale=2/(chains$extra[,j]+usi))
      if(dist=="Skew-Student-t")
        us[places] <- rgamma(n,shape=(chains$extra[,j]+2*k)/2,scale=2/(chains$extra[,j]+usi+rowSums(matrix(zs[places,]^2,n,k))))
      if(dist=="Slash"){
        u0 <- pgamma(q=1,shape=(chains$extra[,j]+k)/2,scale=2/usi)
        us[places] <- qgamma(p=runif(n)*u0,shape=(chains$extra[,j]+k)/2,scale=2/usi)
        us[places] <- ifelse(us[places]<.Machine$double.xmin,.Machine$double.xmin,us[places])
      }
      if(dist=="Contaminated normal"){
        a <- chains$extra[1,j]*chains$extra[2,j]^(k/2)*exp(-chains$extra[2,j]*usi/2)
        b <- (1-chains$extra[1,j])*exp(-usi/2)
        us[places] <- ifelse(runif(n)<=a/(a+b),chains$extra[2,j],1)
      }
      usi <- us[places]
      zsi <- matrix(0,n,k)
      if(dist %in% c("Skew-normal","Skew-Student-t")){
        ais <- matrix(chains$delta[,j],k,n)*tcrossprod(Sigmanew2[[i]],y2z-X%*%chains[[i]]$location[,(1+(j-1)*k):(j*k)])
        Ais <- diag(k) + matrix(chains$delta[,j],k,k)*Sigmanew2[[i]]*matrix(chains$delta[,j],k,k,byrow=TRUE)
        Ais2 <- chol2inv(chol(Ais))
        if(k>1) Ais3 <- chol(Ais2)
        ais <- crossprod(Ais2,ais)
        if(k==1) zsi <- matrix(unlist(lapply(1:n,myf1)),n,k,byrow=TRUE)
        else zsi <- matrix(unlist(lapply(1:n,myf2)),n,k,byrow=TRUE)
        zs[places,] <- zsi
        y <- matrix(data[[i]]$y[places,] - matrix(chains$delta[,j],n,k,byrow=TRUE)*zsi,n,k)
      }
      Xu <- matrix(sqrt(usi),n,s)*X
      yu <- matrix(sqrt(usi),n,k)*y
      A <- chol2inv(chol(Sigmarinv + crossprod(Xu)))
      M <- crossprod(A,(crossprod(Xu,yu) + crossprod(Sigmarinv,mu0s)))
      betanew <- M + crossprod(chol(A),matrix(rnorm(s*k),s,k))%*%chol(chains[[i]]$scale[,(1+(j-1)*k):(j*k)])
      chains[[i]]$location <- cbind(chains[[i]]$location,betanew)
      Omega <- Omega0 + crossprod(betanew-mu0s,Sigmarinv)%*%(betanew-mu0s) + crossprod(yu-Xu%*%betanew)
      Omegachol <- chol(chol2inv(chol(Omega)))
      Sigmanew2[[i]] <- tcrossprod(crossprod(Omegachol,matrix(rnorm(k*(prior$tau0+s+n)),k,prior$tau0+s+n)))
      if(j < n.burnin) Sigmanew <- chol2inv(chol(Sigmanew2[[i]])) else Sigmanew <- Sigmanew2[[i]]
      chains[[i]]$scale <- cbind(chains[[i]]$scale,Sigmanew)
      if(dist=="Contaminated normal"){
        resu <- y-X%*%betanew
        ss[places] <- colSums(t(resu)*tcrossprod(Sigmanew2[[i]],resu))
      }
      if(dist %in% c("Skew-normal","Skew-Student-t")){
        bis <- bis + rowSums(matrix(matrix(usi,k,n,byrow=TRUE)*t(zsi)*tcrossprod(Sigmanew2[[i]],(y2z-X%*%betanew)),k,n))
        Bis <- Bis + Reduce(`+`, lapply(1:n,function(x) usi[x]*(matrix(zsi[x,],k,k)*Sigmanew2[[i]]*matrix(zsi[x,],k,k,byrow=TRUE))))
      }
    }
    if(dist=="Hyperbolic"){
      etanew <- length(us)*(resto - (1/2)*(nus^2*mean(1/us)))
      etanew <- exp(etanew - max(etanew))
      probs <- num1*(etanew[2:length(nus)] + etanew[1:(length(nus)-1)])/2
      etanew <- sample(x=num2,size=1,prob=probs/sum(probs))
      chains$extra <- cbind(chains$extra,matrix(etanew,1,1))
    }
    if(dist %in% c("Student-t","Skew-Student-t")){
      etanew <- length(us)*resto + (nus/2)*sum(log(us)-us)
      etanew <- exp(etanew - max(etanew))
      probs <- num1*(etanew[2:length(nus)] + etanew[1:(length(nus)-1)])/2
      etanew <- sample(x=num2,size=1,prob=probs/sum(probs))
      chains$extra <- cbind(chains$extra,matrix(etanew,1,1))
    }
    if(dist=="Slash"){
      etanew <- rgamma(1,shape=prior$gamma0+length(us),scale=2/(prior$eta0-sum(log(us))))
      chains$extra <- cbind(chains$extra,matrix(etanew,1,1))
    }
    if(dist=="Contaminated normal"){
      a <- us==chains$extra[2,j]
      etanew1 <- max(0.01,rbeta(1,shape1=prior$gamma01+sum(a),shape2=prior$eta01+length(a)-sum(a)))
      u0 <- pgamma(q=1,shape=(sum(a)*k + prior$gamma02)/2,scale=2/(prior$eta02 + sum(ss*a)))
      etanew2 <- max(0.01,qgamma(p=runif(1)*u0,shape=(sum(a)*k + prior$gamma02)/2,scale=2/(prior$eta02 + sum(ss*a))))
      chains$extra <- cbind(chains$extra,matrix(c(etanew1,etanew2),2,1))
    }
    if(dist %in% c("Skew-normal","Skew-Student-t")){
      Bis2 <- chol2inv(chol(Bis + diag(k)/prior$phi0))
      chains$delta <- cbind(chains$delta,Bis2%*%bis + crossprod(chol(Bis2),matrix(rnorm(k),k,1)))
    }else chains$delta <- cbind(chains$delta,matrix(0,k,1))
    if(regim > 1){
      a0 <- (thresholds.chains[,j] - t0)/(t1 - t0)
      if(length(a0) > 1) a0 <- c(a0[1],diff(a0))
      tp <- 100000
      a0 <- tp*c(a0,1-sum(a0))
      ind <- TRUE
      while(ind){
        r0 <- rgamma(length(a0),shape=a0,scale=rep(1,length(a0)))
        r0 <- r0/sum(r0)
        thresholds.new <- t0 + cumsum(r0[-length(a0)])*(t1 - t0)
        if(length(unique(thresholds.new))==(length(r0)-1)){
          indl <- table(cut(Z[(ps+1-hs):(nrow(D)-hs),],breaks=c(-Inf,sort(thresholds.new),Inf),labels=FALSE))
          if(length(indl) == regim) ind <- any(indl < tol)
        }
      }
      a <- min(1,exp(Loglik(hs,thresholds.new) - Loglik(hs,thresholds.chains[,j]) +
              sum((tp*r0-1)*log(a0/tp) - lgamma(tp*r0)) + lgamma(sum(tp*r0)) - sum((a0-1)*log(r0) - lgamma(a0)) - lgamma(sum(a0))))
      if(runif(1) > a){
        thresholds.new <- thresholds.chains[,j]
        ta <- c(ta,0)
      }
      else ta <- c(ta,1)
      thresholds.chains <- cbind(thresholds.chains,thresholds.new)
      resul <- vector()
      for(h in prior$hmin:prior$hmax) resul <- c(resul,Loglik(h,thresholds.new))
      resul <- resul-max(resul)
      resul <- exp(resul)/sum(exp(resul))
      hs.chains <- cbind(hs.chains,sample(prior$hmin:prior$hmax,size=1,prob=resul))
    }
    setTxtProgressBar(bar,j)
  }
  for(i in 1:regim){
    chains[[i]]$location <- matrix(chains[[i]]$location[,ids],nrow=nrow(chains[[i]]$location),ncol=n.sim*k)
    rownames(chains[[i]]$location) <- name[[i]]
    chains[[i]]$scale <- matrix(chains[[i]]$scale[,ids],nrow=k,ncol=k*n.sim)
    colnames(chains[[i]]$scale) <- rep(colnames(D),n.sim)
    rownames(chains[[i]]$scale) <- colnames(D)
  }
  if(dist %in% c("Skew-Student-t","Student-t","Hyperbolic","Slash","Contaminated normal"))
    chains$extra <- matrix(chains$extra[,n.burnin + seq(1,n.sim*n.thin,n.thin)],nrow=ifelse(dist=="Contaminated normal",2,1),ncol=n.sim)
  if(dist %in% c("Skew-normal","Skew-Student-t"))
    chains$delta <- matrix(chains$delta[,n.burnin + seq(1,n.sim*n.thin,n.thin)],nrow=k,ncol=n.sim)
  out_ <- list(data=data,chains=chains,n.sim=n.sim,regim=regim,name=name,dist=dist,ps=ps,ars=ars,formula=Formula(formula),Intercept=Intercept,call=match.call(),log=log,response.series=D)
  if(regim > 1){
    out_$chains$thresholds <- matrix(thresholds.chains[,n.burnin + seq(1,n.sim*n.thin,n.thin)],ncol=n.sim)
    out_$chains$h <- hs.chains[,n.burnin + seq(1,n.sim*n.thin,n.thin)]
    out_$threshold.series=Z
    out_$ts=paste0(colnames(Z),".lag(",mean(out_$chains$h),")")
    out_$ta=100*mean(ta[n.burnin+seq(1,n.sim*n.thin)])
  }
  if(max(ars$q) > 0) out_$covariable.series=X2
  class(out_) <- "mtar"
  return(out_)
}
#' @method coef mtar
#' @export
coef.mtar <- function(object,...,FUN=mean){
  k <- ncol(object$data[[1]]$y)
  n.sim <- object$n.sim
  out_ <- list()
  for(i in 1:object$regim){
    out_[[i]] <- list()
    out <- outs <- vector()
    for(j in 1:k){
      temp <- object$chains[[i]]$location[,seq(j,n.sim*k,k)]
      temps <- matrix(matrix(object$chains[[i]]$scale,k,n.sim*k)[,seq(j,n.sim*k,k)],nrow=k)
      out <- cbind(out,apply(temp,1,FUN))
      outs <- cbind(outs,apply(temps,1,FUN))
    }
    out_[[i]]$location <- out
    out_[[i]]$scale <- outs
    colnames(out_[[i]]$location) <- colnames(out_[[i]]$scale) <- rownames(out_[[i]]$scale) <- colnames(object$data[[1]]$y)
  }
  if(object$regim > 1){
    out_$delay <- as.integer(round(apply(matrix(object$chains$h,1,n.sim),1,FUN)))
    out_$thresholds <- matrix(apply(object$chains$thresholds,1,FUN),object$regim-1,1)
    rownames(out_$thresholds) <- paste0("Threshold",1:(object$regim-1))
    colnames(out_$thresholds) <- ""
  }
  if(object$dist %in% c("Skew-normal","Skew-Student-t")){
    out_$skewness <- matrix(apply(object$chains$delta,1,FUN),k,1)
    rownames(out_$skewness) <- paste0("delta",1:k)
    colnames(out_$skewness) <- ""
  }
  if(object$dist %in% c("Slash","Contaminated normal","Student-t","Hyperbolic","Skew-Student-t")){
    out_$extra <- matrix(apply(object$chains$extra,1,FUN),nrow(object$chains$extra),1)
    rownames(out_$extra) <- paste0("nu",1:nrow(object$chains$extra))
    colnames(out_$extra) <- ""
  }
  return(out_)
}

#' @method summary mtar
#' @export
summary.mtar <- function(object, credible=0.95, digits=max(3, getOption("digits") - 2),...){
  k <- ncol(object$data[[1]]$y)
  n.sim <- object$n.sim
  out_ <- list()
  resumen <- function(x){
    x <- matrix(x,ifelse(is.null(nrow(x)),1,nrow(x)),ifelse(is.null(ncol(x)),length(x),ncol(x)))
    y <- matrix(0,nrow(x),4)
    y[,1] <- rowMeans(x)
    y[,2] <- apply(x,1,function(x) min(mean(sign(median(x))*x > 0),1-1/200000))
    y[,2] <- 2*(1 - y[,2])
    ks <- seq(credible,1,length=n.sim*(1-credible))
    lis <- t(apply(x,1,quantile,probs=ks-credible))
    lss <- t(apply(x,1,quantile,probs=ks))
    dif <- apply(abs(lss-lis),1,which.min)
    y[,3] <- lis[cbind(1:nrow(x),dif)]
    y[,4] <- lss[cbind(1:nrow(x),dif)]
    colnames(y) <- c("   Mean"," 2(1-PD) ","HDI_low","HDI_high")
    return(y)
  }
  if(object$regim > 1){
    thresholds <- matrix(round(resumen(matrix(object$chains$thresholds,nrow=object$regim-1)),digits=digits)[,c(1,3,4)],ncol=3)
    h <- as.integer(round(mean(object$chains$h)))
    thresholds1 <- paste0(c("(-Inf",paste0("(",round(thresholds[,1],digits=digits))),",",c(paste0(round(thresholds[,1],digits=digits),"]"),"Inf)"))
    thresholds2 <- paste0(c("(-Inf",paste0("(",round(thresholds[,2],digits=digits))),",",c(paste0(round(thresholds[,2],digits=digits),"]"),"Inf)"))
    thresholds3 <- paste0(c("(-Inf",paste0("(",round(thresholds[,3],digits=digits))),",",c(paste0(round(thresholds[,3],digits=digits),"]"),"Inf)"))
    d <- data.frame(cbind(thresholds1,thresholds2,thresholds3))
    rownames(d) <- paste("Regime",1:nrow(d))
    colnames(d) <- rep(" ",3)
  }
  cat("\nResponse          :",ifelse(length(colnames(object$data[[1]]$y))==1,colnames(object$data[[1]]$y),paste(colnames(object$data[[1]]$y),collapse="    |    ")))
  if(object$regim > 1) cat("\nThreshold series  :",object$ts,"(Mean)")
  cat("\nError distribution:",object$dist)
  cat("\n\n")
  if(object$regim > 1){
    cat("\nThresholds (Mean, HDI_low, HDI_high)\n")
    print(d)
    out_$thresholds <- d
    out_$location <- list()
    out_$scale <- list()
  }
  for(i in 1:object$regim){
    out <- outs <- vector()
    for(j in 1:k){
      temp <- object$chains[[i]]$location[,seq(j,n.sim*k,k)]
      temps <- matrix(matrix(object$chains[[i]]$scale,k,n.sim*k)[,seq(j,n.sim*k,k)],nrow=k)
      if(j > 1){
        out <- cbind(out,matrix(0,nrow(out),1),round(resumen(temp),digits=digits))
        outs <- cbind(outs,round(resumen(temps)[,c(1,3,4)],digits=digits+1))
      }
      else{
        out <- round(resumen(temp),digits=digits)
        outs <- round(resumen(temps)[,c(1,3,4)],digits=digits+1)
      }
    }
    outs <- matrix(outs,k,3*k)
    rownames(out) <- object$name[[i]]
    outs <- matrix(cbind(outs,0),k,3*k+1)
    outs <- outs[,c(seq(1,3*k,3),3*k+1,seq(2,3*k,3),3*k+1,seq(3,3*k,3))]
    outs <- matrix(outs,k,length(outs)/k)
    rownames(outs) <- colnames(object$data[[1]]$y)
    colnames(outs) <- c(rownames(outs),"",rownames(outs),"",rownames(outs))
    cat("\n\nRegime",i,":")
    cat("\nAutoregressive coefficients\n")
    print(format(out, justify = "right", format = "+/-", zero.print="   |   "), quote=FALSE)
    cat("\nScale parameter (Mean, HDI_low, HDI_high)\n")
    print(format(outs, justify = "right", format = "+/-", zero.print="   ."), quote=FALSE)
    out_$location[[i]] <- out
    out_$scale[[i]] <- outs
  }
  if(object$dist %in% c("Skew-normal","Skew-Student-t")){
    out <- round(resumen(object$chains$delta),digits=digits+1)
    rownames(out) <- paste0(paste0("delta",1:k),paste0(rep("",max(nchar(object$name[[1]]))-6),collapse=" "))
    cat("\n\nSkewness parameter","\n")
    print(format(out, justify = "right", flag="+", zero.print="   .   "), quote=FALSE)
    out_$delta <- out
  }
  if(object$dist %in% c("Slash","Contaminated normal","Student-t","Hyperbolic","Skew-Student-t")){
    out <- round(resumen(object$chains$extra),digits=digits+1)
    out[,2] <- 0
    if(object$dist %in% c("Slash","Student-t","Hyperbolic","Skew-Student-t")) rownames(out)[nrow(out)] <- paste0("nu",paste0(rep("",max(nchar(object$name[[1]]))-1),collapse=" "))
    else rownames(out)[nrow(out):(nrow(out)-1)] <- paste0(c("nu2","nu1"),paste0(rep("",max(nchar(object$name[[1]]))-2),collapse=" "))
    cat("\n\nExtra parameter","\n")
    print(format(out, justify = "right", flag="+", zero.print="   .   "), quote=FALSE)
    out_$extra <- out
  }
  cat("\n\n")
  return(invisible(out_))
}
#'
#' @method plot mtar
#' @export
plot.mtar <- function(x,...,identify){
  n.sim <- x$n.sim
  dist <- x$dist
  k <- ncol(x$data[[1]]$y)
  out <- matrix(NA,1,k)
  idsi <- vector()
  out_2 <- matrix(NA,1,2*k)
  if(x$regim > 1) lims <- (x$ps+1):(length(x$threshold.series))
  if(x$dist %in% c("Student-t","Hyperbolic","Slash","Contaminated normal"))
    extra <- rowMeans(matrix(x$chains$extra,ncol=n.sim))
  ids <- 1:nrow(x$data[[1]]$X)
  for(i in 1:x$regim){
    location <- matrix(0,nrow(x$chains[[i]]$location),k)
    scale <- matrix(0,k,k)
    for(j in 1:k){
      location[,j] <- rowMeans(matrix(x$chains[[i]]$location[,seq(j,n.sim,by=k)],nrow(location),n.sim))
      scale[,j] <- rowMeans(matrix(x$chains[[i]]$scale[,seq(j,n.sim,by=k)],k,n.sim))
    }
    if(x$regim > 1){
      h <- as.integer(round(mean(x$chains$h)))
      thresholds <- rowMeans(matrix(x$chains$thresholds,ncol=n.sim))
      Z <- x$threshold.series[lims-h]
      regs <- cut(Z,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
    }else regs <- matrix(1,nrow(x$data[[i]]$y),1)
    places <- regs==i
    y <- x$data[[i]]$y[places,]
    X <- x$data[[i]]$X[places,]
    resu <- y-X%*%location
    ei <- eigen(scale,symmetric=TRUE)
    resu <- resu%*%(ei$vectors%*%diag(1/sqrt(as.vector(ei$values)))%*%t(ei$vectors))
    out <- rbind(out,resu)
    out_2 <- rbind(out_2,cbind(y,X%*%location))
    idsi <- c(idsi,ids[places])
  }
  out <- matrix(out,nrow(x$data[[1]]$X),k)
  TT <- 50000
  sim <- matrix(rnorm(TT*k),TT,k)
  if(dist=="Gaussian") u <- 1
  if(dist=="Student-t") u <- 1/rgamma(TT,shape=extra/2,rate=extra/2)
  if(dist=="Slash") u <- 1/rbeta(TT,shape1=extra/2,shape2=1)
  if(dist=="Contaminated normal") u <- 1/(1 - (1-extra[2])*rbinom(TT,1,extra[1]))
  if(dist=="Hyperbolic") u <- rgig(n=TT,lambda=1,chi=1,psi=extra^2)
  if(dist=="Laplace") u <- rexp(TT,rate=1/8)
  sim2 <- sort(rowSums(sim^2)*u)
  out2 <- apply(matrix(rowSums(out*out),nrow(out),1),1,function(x) mean(sim2<=x))
  out3 <- ifelse(out2>=0.5,1-out2,out2)
  out3 <- qnorm(ifelse(.Machine$double.xmin>=out3,.Machine$double.xmin,out3))*ifelse(out2>0.5,-1,1)
  sim <- sim*matrix(sqrt(u),TT,k)
  for(i in 1:k){
     temp <- apply(matrix(out[,i],nrow(out),1),1,function(x) mean(sim[,i]<=x))
     temp2 <- ifelse(temp>=0.5,1-temp,temp)
     out[,i] <- qnorm(ifelse(.Machine$double.xmin>=temp2,.Machine$double.xmin,temp2))*ifelse(temp>0.5,-1,1)
  }
  idsx <- sort(idsi,index=TRUE)$ix
  out3 <- matrix(cbind(out,out3)[idsx,],length(out3),k+1)
  out_2 <- na.omit(out_2);out_2 <- matrix(out_2[idsx,],nrow(out_2),ncol(out_2))
  row.names(out3) <- row.names(x$data[[1]]$y)
  row.names(out_2) <- row.names(x$data[[1]]$y)
  nano_ <- list(...)
  nano_$y <- out3[,k+1]
  nano_$type <- "p"
  if(is.null(nano_$pch)) nano_$pch <- 20
  dev.new()
  outm <- do.call("qqnorm",nano_)
  abline(0,1,lty=3)
  if(!missingArg(identify)) identify(outm$x,outm$y,n=max(1,floor(abs(identify))),labels=row.names(out3))
  return(invisible(list(plot=out3,output=out_2)))
}

#' @title Forecasting of a multivariate TAR model.
#' @description This function computes forecasting from a fitted multivariate TAR model.
#' @param object an object of the class \emph{mtar}.
#' @param data an (optional) data frame, list or environment (or object coercible by
#'             \link{as.data.frame} to a data frame) containing the future values of the threshold
#'             series as well as the exogenous series in the model.
#'             If not found in data, the variables are taken from \code{environment(formula)},
#'             typically the environment from which \code{mtar} is called.
#' @param credible an (optional) value for the level of the credible intervals. By default, \code{credible} is set to 0.95.
#' @param row.names an vector that allows the user to name the time point to
#'                  which each row in the data set \code{data} corresponds.
#' @param out.of.sample an (optional) logical variable. If \code{TRUE}, then the log-score is computed, which is a measurement to assess density forecasts. Therefore, the data.frame specified in the argument \code{data} must to include the true values of the output series.
#' @param setar an (optional) positive integer indicating the component of the output series which is
#'              the threshold variable. By default,\code{setar} is set to \code{NULL}, which indicates that the fitted model is not a SETAR.
#' @return a list with the following component
#' \tabular{ll}{
#' \code{ypred}   \tab a matrix with the results of the forecasting,\cr
#' \tab \cr
#' \code{summary} \tab a matrix with the mean and credible intervals of the forecasting,\cr
#' }
#' @export forecasting
#' @references Nieto, F.H. (2005) Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data.
#'             Communications in Statistics - Theory and Methods, 34, 905-930.
#' @references Romero, L.V. and Calderon, S.A. (2021) Bayesian estimation of a multivariate TAR model when the noise
#'             process follows a Student-t distribution. Communications in Statistics - Theory and Methods, 50, 2508-2530.
#' @references Calderon, S.A. and Nieto, F.H. (2017) Bayesian analysis of multivariate threshold autoregressive models
#'             with missing data. Communications in Statistics - Theory and Methods, 46, 296-318.
#' @references Karlsson, S. (2013) Chapter 15-Forecasting with Bayesian Vector Autoregression. In Elliott, G. and
#'             Timmermann, A. Handbook of Economic Forecasting, Volume 2, 791–89, Elsevier.
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1 <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'              dist="Gaussian", ars=list(p=c(1,1,2)), n.burnin=100,
#'              n.sim=3000, n.thin=2)
#' out1 <- forecasting(fit1,data=subset(returns,Date >= "2016-03-20"),row.names=Date)
#' out1$summary
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2 <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'              dist="Gaussian", ars=list(p=c(5,5,5)), n.burnin=2000,
#'              n.sim=3000, n.thin=2)
#' out2 <- forecasting(fit2,data=subset(riverflows,Date >= "2009-04-09"),row.names=Date)
#' out2$summary
#' }
#'
forecasting <- function(object,data,out.of.sample=FALSE,credible=0.95,row.names,setar=NULL){
  if(class(object)[1]!="mtar")
    stop("Only mtar-type objects are supported!!",call.=FALSE)
  regim <- object$regim
  gendist <- function(dist,Sigma,extra,delta){
    Sigma <- as.matrix(Sigma)
    nor <- matrix(rnorm(ncol(Sigma)),1,ncol(Sigma))%*%chol(Sigma)
    if(dist %in% c("Gaussian","Skew-normal")) u <- 1
    if(dist %in% c("Student-t","Skew-Student-t")) u <- 1/rgamma(1,shape=extra/2,rate=extra/2)
    if(dist=="Slash") u <- 1/rbeta(1,shape1=extra/2,shape2=1)
    if(dist=="Contaminated normal") u <- 1/(1 - (1-extra[2])*rbinom(1,1,extra[1]))
    if(dist=="Hyperbolic") u <- rgig(n=1,lambda=1,chi=1,psi=extra^2)
    if(dist=="Laplace") u <- rexp(1,rate=1/8)
    if(dist %in% c("Skew-normal","Skew-Student-t"))
      out_ <- sqrt(u)*matrix(delta*qnorm(0.5 + runif(ncol(Sigma))*0.5) + nor,1,ncol(Sigma)) else out_ <- nor*sqrt(u)
    return(out_)
  }
  if(out.of.sample){
    Lik <- function(dist,y,M,Sigma,delta,log,nu){
      resu <- matrix(y-M,1,k)
      if(dist %in% c("Skew-normal","Skew-Student-t")){
        D <- as.vector(delta)*diag(k); A <- chol2inv(chol(Sigma + D**2))
        out0 <- resu%*%tcrossprod(A,resu)
        muv <- resu%*%(A*matrix(delta,k,k,byrow=TRUE))
        Sigmav <- diag(k) - matrix(delta,k,k)*A*matrix(delta,k,k,byrow=TRUE)
        if(dist=="Skew-normal"){
          out <- out0 + k*log(2*pi)
          if(k > 1) out <- out - 2*log(pmvnorm(lower=-muv,upper=rep(Inf,k),sigma=Sigmav))
          else out <- out - 2*log(pnorm(-muv/sqrt(Sigmav),lower.tail=FALSE))
        }
        if(dist=="Skew-Student-t"){
          out <- (nu+k)*log(1 + out0/nu) - 2*lgamma((nu+k)/2) + 2*lgamma(nu/2) + k*log(nu*pi)
          muv <- muv*sqrt((nu + k)/(nu + out0))
          if(k > 1) out <- out - 2*log(pmvt(lower=-muv,upper=rep(Inf,k),sigma=Sigmav,df=round(nu,0)+k))
          else out <- out - 2*log(pt(-muv/sqrt(Sigmav),df=nu,lower.tail=FALSE))
        }
      }else{
        out <- resu%*%tcrossprod(chol2inv(chol(Sigma)),resu)
        if(dist=="Laplace")
          out <- -2*log(besselK(sqrt(out)/2,(2-k)/2)) - ((2-k)/2)*log(out)	+ (k+2)*log(2)
        if(dist=="Student-t")
          out <- (nu+k)*log(1 + out/nu) - 2*lgamma((nu+k)/2) + 2*lgamma(nu/2) + k*log(nu/2)
        if(dist=="Contaminated normal")
          out <- -2*log(nu[1]*exp(-out*nu[2]/2)*nu[2]^(k/2) + (1-nu[1])*exp(-out/2))
        if(dist=="Slash")
          out <- ifelse(out==0,-2*log(nu/(nu+k)),-2*lgamma((k+nu)/2) + (k+nu)*log(out/2) -2*log((nu/2)*pgamma(1,shape=(k+nu)/2,rate=out/2)))
        if(dist=="Hyperbolic")
          out <- -2*log(besselK(nu*sqrt(1+out),(2-k)/2)) -((2-k)/2)*log(1+out) - k*log(nu) +2*log(besselK(nu,1))
        out <- out + k*log(2*pi)
      }
      out <- -(out + log(det(Sigma)))/2
      if(log) out <- out + sum(y)
      return(exp(out))
    }
  }
  mmf <- match.call(expand.dots = FALSE)
  m <- match(c("data", "row.names"), names(mmf), 0)
  mmf <- mmf[c(1,m)]
  mmf$drop.unused.levels <- TRUE
  mmf[[1]] <- as.name("model.frame")
  mmf <- eval(mmf, parent.frame())
  row.names <- as.vector(as.character(model.extract(mmf,row.names)))
  if(out.of.sample){
    mx <- model.part(Formula(object$formula), data = mmf, rhs = 1, terms = TRUE)
    D <- model.matrix(mx, data = mmf)
    if(attr(terms(mx),"intercept")){
      Dnames <- colnames(D)
      D <- matrix(D[,-1],ncol=length(Dnames)-1)
      colnames(D) <- Dnames[-1]
    }
    ytrue <- D
  }
  if(regim > 1 & is.null(setar)){
    mz <- model.part(Formula(object$formula), data = mmf, rhs = 2, terms = TRUE)
    Z <- model.matrix(mz, data = mmf)
    if(attr(terms(mz),"intercept")){
      Znames <- colnames(Z)
      Z <- as.matrix(Z[,-1])
    }
    Z <- rbind(object$threshold.series,Z)
  }
  pasos <- length(row.names)
  if(max(object$ars$q) > 0){
    mx2 <- model.part(Formula(object$formula), data = mmf, rhs = 3, terms = TRUE)
    X2 <- model.matrix(mx2, data = mmf)
    if(attr(terms(mx2),"intercept")){
      X2names <- colnames(X2)
      X2 <- as.matrix(X2[,-1])
      colnames(X2) <- X2names[-1]
    }
    X2 <- rbind(object$covariable.series,X2)
  }
  y <- object$response.series
  n <- nrow(y)
  k <- ncol(y)
  ysim <- rbind(matrix(y,n,k*object$n.sim),matrix(0,pasos,k*object$n.sim))
  if(out.of.sample){
    prek <- matrix(0,pasos,object$n.sim)
    ytrue <- rbind(matrix(y,n,k),matrix(ytrue,nrow(ytrue),k))
  }
  for(i in 1:object$n.sim){
    h <- object$chains$h[i]
    thresholds <- object$chains$thresholds[,i]
    for(j in (n+1):(n+pasos)){
      if(regim > 1){
        if(!is.null(setar)) iz <- ysim[j-h,(i-1)*k + setar] else iz <- Z[j-h]
        regs <- cut(iz,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
      }else regs <- 1
      X <- 1;for(l in 1:object$ars$p[regs]) X <- c(X,ysim[j-l,((i-1)*k+1):(i*k)])
      if(!object$Intercept) X <- X[-1]
      if(object$ars$q[regs] > 0) for(l in 1:object$ars$q[regs]) X <- c(X,X2[j-l,])
      if(object$ars$d[regs] > 0) for(l in 1:object$ars$d[regs]) if(!is.null(setar)) X <- c(X,ysim[j-l,(i-1)*k + setar]) else X <- c(X,Z[j-l,])
      M <- matrix(X,1,length(X))%*%object$chains[[regs]]$location[,((i-1)*k+1):(i*k)]
      if(object$dist %in% c("Gaussian","Laplace")) extrai <- 0 else extrai <- object$chains$extra[,i]
      if(object$dist %in% c("Skew-normal","Skew-Student-t")) deltai <- object$chains$delta[,i]
      else deltai <- rep(0,k)
      ysim[j,((i-1)*k+1):(i*k)] <- M + gendist(object$dist,object$chains[[regs]]$scale[,((i-1)*k+1):(i*k)],extrai,deltai)
      if(out.of.sample){
        Xous <- 1;for(l in 1:object$ars$p[regs]) Xous <- c(Xous,ytrue[j-l,])
        if(!object$Intercept) Xous <- Xous[-1]
        if(object$ars$q[regs] > 0) for(l in 1:object$ars$q[regs]) Xous <- c(Xous,X2[j-l,])
        if(object$ars$d[regs] > 0) for(l in 1:object$ars$d[regs]) if(!is.null(setar)) Xous <- c(Xous,ytrue[j-l,setar]) else Xous <- c(Xous,Z[j-l,])
        M <- matrix(Xous,1,length(Xous))%*%object$chains[[regs]]$location[,((i-1)*k+1):(i*k)]
        prek[j-n,i] <- Lik(object$dist,ytrue[j,],M,object$chains[[regs]]$scale[,((i-1)*k+1):(i*k)],deltai,object$log,extrai)
      }
    }
  }
  ysim <- matrix(ysim[-c(1:n),],pasos,k*object$n.sim)
  colnames(ysim) <- rep(colnames(y),object$n.sim)
  if(object$log) ysim <- exp(ysim)
  out_ <- vector()
  predi <- function(x){
    x <- matrix(x,1,length(x))
    ks <- seq(credible,1,length=object$n.sim*(1-credible))
    lis <- t(apply(x,1,quantile,probs=ks-credible))
    lss <- t(apply(x,1,quantile,probs=ks))
    dif <- apply(abs(lss-lis),1,which.min)
    out_ <- c(mean(x),lis[dif],lss[dif])
    names(out_) <- c("Mean","HDI_Low","HDI_high")
    return(out_)
  }
  for(i in 1:k) out_ <- cbind(out_,t(apply(matrix(ysim[,seq(i,k*object$n.sim,k)],pasos,object$n.sim),1,predi)))
  rownames(out_) <- row.names
  rownames(ysim) <- row.names
  out__ <- list(ypred=ysim,summary=out_,y=y)
  if(out.of.sample) out__$log.score <- -log(rowMeans(prek))
  class(out__) <- "forecastmtar"
  return(out__)
}
#' @method plot forecastmtar
#' @export
plot.forecastmtar <- function(x,...,n,historical=list(),forecasts=list(),forecasts.PI=list(),main=NULL){
  k <- ncol(x$y)
  out_ <- x$summary
  pasos <- nrow(out_)
  y <- matrix(x$y,nrow(x$y),ncol(x$y))
  if(missingArg(n)) n <- nrow(y) else n <- ceiling(abs(n))
  y2 <- rbind(matrix(y[(nrow(y)-n+1):nrow(y),],ncol=k),matrix(out_[,seq(1,3*k,3)],ncol=k))
  forecasts$xlim <- historical$xlim <- c(1,n+pasos)
  forecasts$x <- (n+1):(n+pasos)
  forecasts.PI$ylab <- forecasts.PI$xlab <- historical$xlab <- historical$ylab <- ""
  if(is.null(historical$type)) historical$type <- "l"
  if(is.null(historical$lty)) historical$lty <- 1
  if(is.null(historical$col)) historical$col <- "black"
  forecasts$xlim <- c(1,n+pasos)
  if(is.null(forecasts$type)) forecasts$type <- "l"
  if(is.null(forecasts$lty)) forecasts$lty <- 1
  if(is.null(forecasts$col)) forecasts$col <- "blue"
  if(is.null(forecasts$xlab)) forecasts$xlab <- ""
  if(is.null(forecasts$ylab)) forecasts$ylab <- ""
  if(is.null(forecasts.PI$density)) forecasts.PI$density <- NA
  if(is.null(forecasts.PI$col)) forecasts.PI$col <- "light gray"
  for(i in 1:k){
    dev.new()
    historical$ylim <- forecasts.PI$ylim <- forecasts$ylim <- range(y2[,i],out_[,1:3+3*(i-1)])
    historical$x <- y[(nrow(y)-n+1):nrow(y),i]
    if(is.null(main)) historical$main <- colnames(x$y)[i]
    else historical$main <- main[i]
    do.call("plot",historical)
    xs <- c((n+1):(n+pasos),(n+pasos):(n+1))
    ys <- c(out_[,2+3*(i-1)],out_[nrow(out_):1,3+3*(i-1)])
    par(new=TRUE)
    forecasts.PI$x <- xs
    forecasts.PI$y <- ys
    do.call("polygon",forecasts.PI)
    forecasts$y <- out_[,1+3*(i-1)]
    forecasts$main <- colnames(y)[i]
    par(new=TRUE)
    do.call("plot",forecasts)
  }
}
#'
#' @title Converts chains from the Bayesian estimation of a multivariate TAR model to a mcmc object.
#' @description This function converts the chains obtained from the Bayesian estimation of a multivariate TAR model to a \code{mcmc} object to be analyzed with the \pkg{coda} package.
#' @param object an object of the class \emph{mtar}.
#'
#' @return a \code{mcmc}-type object.
#' @export convert
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1 <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'              dist="Gaussian", ars=list(p=c(1,1,2)), n.burnin=100,
#'              n.sim=3000, n.thin=2)
#' chains1 <- convert(fit1)
#' summary(chains1$location[[1]])
#' plot(chains1$location[[1]])
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2 <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'              dist="Gaussian", ars=list(p=c(5,5,5)), n.burnin=2000,
#'              n.sim=3000, n.thin=2)
#' chains2 <- convert(fit2)
#' summary(chains2$location[[2]])
#' plot(chains2$location[[2]])
#' }
#'
convert <- function(object){
  if(class(object)[1]!="mtar")
    stop("Only mtar-type objects are supported!!",call.=FALSE)
  k <- ncol(object$data[[1]]$y)
  yn <- colnames(object$data[[1]]$y)
  n.sim <- object$n.sim
  out <- list()
  out$location <- list()
  out$scale <- list()
  for(r in 1:object$regim){
    datos <- matrix(object$chains[[r]]$location,nrow=nrow(object$chains[[r]]$location),ncol=n.sim*k)
    out_ <- vector()
    for(i in 1:k){
      temp <- t(matrix(datos[,seq(i,n.sim*k,k)],nrow=nrow(datos),ncol=n.sim))
      colnames(temp) <- paste0(yn[i],":",rownames(object$chains[[r]]$location))
      out_ <- cbind(out_,temp)
    }
    out$location[[r]] <- mcmc(out_)
    datos <-  matrix(object$chains[[r]]$scale,nrow=k,ncol=n.sim*k)
    out_ <- vector()
    for(i in 1:k){
      for(j in i:k){
        temp <- t(matrix(datos[i,seq(j,n.sim*k,k)],nrow=1,ncol=n.sim))
        colnames(temp) <- paste0(yn[i],".",yn[j])
        out_ <- cbind(out_,temp)
      }
    }
    out$scale[[r]] <- mcmc(out_)
  }
  if(!(object$dist %in% c("Gaussian","Laplace"))){
    datos <-  matrix(object$chains$extra,nrow=nrow(object$chains$extra),ncol=n.sim)
    out_ <- t(datos)
    if(nrow(datos)==1) colnames(out_) <- "nu" else colnames(out_) <- paste("nu",1:nrow(datos))
    out$extra <- mcmc(out_)
  }
  if(object$dist %in% c("Skew-normal","Skew-Student-t")){
    datos <-  matrix(object$chains$delta,nrow=nrow(object$chains$delta),ncol=n.sim)
    out_ <- t(datos)
    if(nrow(datos)==1) colnames(out_) <- "delta" else colnames(out_) <- paste("delta",1:nrow(datos))
    out$skewness <- mcmc(out_)
  }
  if(object$regim > 1){
    datos <-  matrix(object$chains$thresholds,nrow=nrow(object$chains$thresholds),ncol=n.sim)
    out_ <- t(datos)
    if(nrow(datos)==1) colnames(out_) <- "threshold" else colnames(out_) <- paste("threshold",1:nrow(datos))
    out$thresholds <- mcmc(out_)
    out_ <- matrix(object$chains$h,nrow=n.sim,ncol=1)
    colnames(out_) <- "delay"
    out$delay <- mcmc(out_)
  }
  return(invisible(out))
}
#'
#' @title Simulation of multivariate time series according to a TAR model
#' @description This function simulates multivariate time series according to a user-specified TAR model.
#' @param n a positive integer value indicating the length of the desired output series.
#' @param k a positive integer value indicating the dimension of the desired output series.
#' @param ars a list composed of three objects, namely: \code{p}, \code{q} and \code{d},
#'            each of which corresponds to a vector of \eqn{l} non-negative integers, where \eqn{l} represents the number of regimes in the TAR model.
#' @param Intercept an (optional) logical variable. If \code{TRUE}, then the model includes an intercept.
#' @param delay an (optional) non-negative integer value indicating the delay in the threshold series.
#' @param thresholds a vector with \eqn{l-1} real values sorted ascendingly.
#' @param parms a list with as many sublists as regimes in the user-specified TAR model. Each sublist is composed of two matrices. The first corresponds
#'              to location parameters, while the second corresponds to scale parameters.
#' @param t.series a matrix with the values of the threshold series.
#' @param ex.series a matrix with the values of the multivariate exogenous series.
#' @param dist an (optional) character string which allows the user to specify the multivariate
#'             distribution to be used to describe the behavior of the noise process. The
#'             available options are: Gaussian ("Gaussian"), Student-\eqn{t} ("Student-t"),
#'             Slash ("Slash"), Symmetric Hyperbolic ("Hyperbolic"), Laplace ("Laplace"), and
#'             contaminated normal ("Contaminated normal"). By default,\code{dist} is set to
#'             "Gaussian".
#' @param delta an (optional) vector with the values of the skewness parameters. By default,\code{delta} is set to \code{NULL}.
#' @param extra a value indicating the value of the extra parameter of the noise process distribution, if any.
#' @param setar an (optional) positive integer indicating the component of the output series which should
#'              be the threshold variable. By default,\code{setar} is set to \code{NULL}.
#'
#' @return a \code{data.frame} containing the output series, threshold series (if any), and multivariate exogenous series (if any).
#' @export simtar
#' @examples
#' \donttest{
#' ###### Simulation of a trivariate TAR model with two regimes
#' n <- 2000
#' k <- 3
#' ars <- list(p=c(1,2))
#' Z <- as.matrix(arima.sim(n=n+max(ars$p),list(ar=c(0.5))))
#' Intercept <- TRUE
#' parms <- list()
#' for(j in 1:length(ars$p)){
#'    np <- Intercept + ars$p[j]*k
#'    parms[[j]] <- list()
#'    parms[[j]]$location <- c(ifelse(runif(np*k)<=0.5,1,-1)*rbeta(np*k,shape1=4,shape2=16))
#'    parms[[j]]$location <- matrix(parms[[j]]$location,np,k)
#'    parms[[j]]$scale <- rgamma(k,shape=1,scale=1)*diag(k)
#' }
#' thresholds <- quantile(Z,probs=(0.85 + runif(1)*0.3)*seq(1,length(ars$p)-1)/length(ars$p))
#' out1 <- simtar(n=n,k=k,ars=ars,Intercept=Intercept,parms=parms,
#'                thresholds=thresholds,t.series=Z,dist="Student-t",extra=6)
#' str(out1)
#'
#' fit1 <- mtar(~ Y1 + Y2 + Y3 | t.series, data=out1, ars=ars, dist="Student-t",
#'              nsim=3000, n.burn=2000, n.thin=2)
#' summary(fit1)
#'
#' ###### Simulation of a trivariate VAR model
#' n <- 2000
#' k <- 3
#' ars <- list(p=2)
#' Intercept <- TRUE
#' parms <- list()
#' for(j in 1:length(ars$p)){
#'    np <- Intercept + ars$p[j]*k
#'    parms[[j]] <- list()
#'    parms[[j]]$location <- c(ifelse(runif(np*k)<=0.5,1,-1)*rbeta(np*k,shape1=4,shape2=16))
#'    parms[[j]]$location <- matrix(parms[[j]]$location,np,k)
#'    parms[[j]]$scale <- rgamma(k,shape=1,scale=1)*diag(k)
#' }
#' out2 <- simtar(n=n,k=k,ars=ars,Intercept=Intercept,parms=parms,
#'                dist="Slash",extra=2)
#' str(out2)
#'
#' fit2 <- mtar(~ Y1 + Y2 + Y3, data=out2, ars=ars, dist="Slash", nsim=3000,
#'              n.burn=2000, n.thin=2)
#' summary(fit2)
#' toy <- data.frame(Date=rep(0,10))
#' f2 <- forecasting(fit2, data=toy, row.names=Date)
#' f2$summary
#' plot(f2, n=100)
#'
###### Simulation of a trivariate SETAR model with two regimes
#' n <- 5010
#' k <- 3
#' ars <- list(p=c(1,2))
#' Intercept <- TRUE
#' parms <- list()
#' for(j in 1:length(ars$p)){
#'   np <- Intercept + ars$p[j]*k
#'   parms[[j]] <- list()
#'   parms[[j]]$location <- c(ifelse(runif(np*k)<=0.5,1,-1)*rbeta(np*k,shape1=4,shape2=16))
#'   parms[[j]]$location <- matrix(parms[[j]]$location,np,k)
#'   parms[[j]]$scale <- rgamma(k,shape=1,scale=1)*diag(k)
#' }
#' out3 <- simtar(n=n, k=k, ars=ars, Intercept=Intercept, parms=parms, delay=2,
#'                thresholds=-1, dist="Laplace", setar=2)
#' str(out3)
#'
#' fit3 <- mtar(~ Y1 + Y2 + Y3 | Y2, data=out3, ars=ars, dist="Laplace",
#'              nsim=3000,n.burn=2000, n.thin=2, prior=list(hmin=1))
#' summary(fit3)
#' toy <- data.frame(Date=rep(0,10))
#' f3 <- forecasting(fit3, setar=2, data=toy, row.names=Date)
#' f3$summary
#' plot(f3, n=100)
#' }
simtar <- function(n,k=2,ars=list(p=1),Intercept=TRUE,parms,delay=0,thresholds=0,t.series,ex.series,dist="Gaussian",delta=NULL,extra,setar=NULL){
  n <- ceiling(abs(n))
  k <- ceiling(abs(k))
  if(k<=0 | k!=floor(k)) stop("The argument 'k' must be a positive integer!",call.=FALSE)
  if(is.null(ars$p)) stop("The argument 'ars=list(p=)' must be a positive integer!",call.=FALSE)
  regim <- length(ars$p)
  if(is.null(ars$d) | regim==1) ars$d <- rep(0,length(ars$p))
  if(is.null(ars$q) | missingArg(ex.series)) ars$q <- rep(0,length(ars$p))
  if(regim > 1){
    if(delay<0 | delay!=floor(delay)) stop("The argument 'delay' must be a non-negative integer!",call.=FALSE)
    if(delay==0 & !is.null(setar)) stop("In SETAR models the argument 'delay' must be a positive integer!",call.=FALSE)
    if(missingArg(t.series) & is.null(setar)) stop("The argument 't.series' is required!",call.=FALSE)
    if(missingArg(thresholds)) stop("The argument 'thresholds' is required!",call.=FALSE) else thresholds <- sort(thresholds)
    if(length(thresholds) != regim-1)
      stop(paste0("The length of argument 'thresholds' must be ",regim-1,", and the length of the argument supplied by the user is ",length(thresholds)),call.=FALSE)
  }else{
    delay <- 0
    ars$d <- rep(0,length(ars$p))
    ars$q <- rep(0,length(ars$p))
  }
  ps <- max(ars$p,ars$q,ars$d,delay)
  if(regim > 1 & is.null(setar)){
    if(nrow(t.series)!=n+ps)
      stop(paste0("The number of rows of the argument 't.series' must be ",n+ps,", and the number of rows of the argument supplied by the user is ",nrow(t.series)),call.=FALSE)
  }
  if(max(ars$q)>0){
    if(missingArg(ex.series)) stop("An exogenous series is required!",call.=FALSE)
    if(nrow(ex.series)!=n+ps)
      stop(paste0("The number of rows of the argument 'ex.series' must be ",n+ps,", and the number of rows of the argument supplied by the user is ",nrow(ex.series)),call.=FALSE)
    r <- ncol(ex.series)
  }else r <- 0
  if(!(dist %in% c("Gaussian","Student-t","Hyperbolic","Laplace","Slash","Contaminated normal","Skew-normal","Skew-Student-t")))
    stop("Only 'Gaussian', 'Student-t', 'Hyperbolic', 'Laplace', 'Slash', 'Contaminated normal', 'Skew-normal' and 'Skew-Student-t' distributions are supported!",call.=FALSE)
  if(dist %in% c("Student-t","Hyperbolic","Slash","Contaminated normal","Skew-Student-t")){
    if(is.null(extra))
      stop("For 'Student-t', 'Hyperbolic', 'Slash', 'Contaminated normal' and 'Skew-Student-t' distributions an extra parameter value must be specified!",call.=FALSE)
    if((dist %in% c("Student-t","Skew-Student-t") & extra<=0) | (dist=="Hyperbolic" & extra<=0) | (dist=="Slash" & extra<=0))
      stop("For 'Student-t', 'Skew-Student-t', 'Hyperbolic' and 'Slash' distributions the extra parameter value must be positive!",call.=FALSE)
    if(dist=="Contaminated normal" & (any(extra<=0) | any(extra>=1)))
      stop("For 'Contaminated normal' distribution the extra parameter must be a bidimensional vector with values in the interval (0,1)!",call.=FALSE)
  }
  if(dist %in% c("Skew-normal","Skew-Student-t") & is.null(delta))
    stop("For 'Skew-normal' and 'Skew-Student-t' distributions an skewness parameter must be specified!",call.=FALSE)

  myseries <- matrix(rnorm((n+ps)*k),n+ps,k)
  for(i in 1:regim){
    if(ncol(parms[[i]]$location)!=k | nrow(parms[[i]]$location)!=(Intercept+ars$p[i]*k+ars$q[i]*r+ars$d[i]))
      stop(paste0("Dimension of location matrix in regime ",i," must be ",(Intercept+ars$p[i]*k+ars$q[i]*r+ars$d[i]),"X",k,"!"),call.=FALSE)
    parms[[i]]$scale2  <- try(chol(parms[[i]]$scale),silent=TRUE)
    if(!is.matrix(parms[[i]]$scale2)) stop(paste0("Scale matrix in regime ",i," is not positive definite!"),call.=FALSE)
  }
  for(i in 1:n){
    current <- ps + i
    if(regim > 1){
      if(!is.null(setar)) t.series.c <- myseries[current-delay,setar]
      else t.series.c <- t.series[current-delay]
      regimeni <- cut(t.series.c,breaks=c(-Inf,thresholds,Inf),labels=1:regim)
    }
    else regimeni <- 1
    if(Intercept) X <- 1 else X <-  vector()
    for(j in 1:ars$p[regimeni]) X <- c(X,myseries[current-j,])
    if(ars$q[regimeni] > 0) for(j in 1:ars$q[regimeni]) X <- c(X,ex.series[current-j,])
    if(ars$d[regimeni] > 0) for(j in 1:ars$d[regimeni]) if(!is.null(setar)) X <- c(X,myseries[current-j,setar]) else X <- c(X,t.series[current-j])
    Theta <- parms[[regimeni]]$location
    mu <- colSums(matrix(X,nrow(Theta),ncol(Theta))*Theta)
    u <- 1
    if(dist %in% c("Student-t","Skew-Student-t"))  u <- 1/rgamma(1,shape=extra/2,rate=extra/2)
    if(dist=="Slash")  u <- 1/rbeta(1,shape1=extra/2,shape2=1)
    if(dist=="Contaminated normal")  if(runif(1)<=extra[1]) u <- 1/extra[2]
    if(dist=="Laplace")  u <- rexp(1,rate=1/8)
    if(dist=="Hyperbolic") u <- rgig(n=1,lambda=1,chi=1,psi=extra^2)
    if(dist %in% c("Skew-normal","Skew-Student-t")){
      offset <- sqrt(u)*matrix(delta*qnorm(0.5 + runif(k)*0.5),k,1)
    }
    else offset <- matrix(0,k,1)
    myseries[current,] <- crossprod(parms[[regimeni]]$scale2,matrix(sqrt(u)*rnorm(k),k,1)) + matrix(mu,k,1) + offset
  }
  for(i in 1:regim) parms[[i]]$scale2 <- NULL
  datos <- data.frame(myseries)
  colnames(datos) <- paste("Y",1:k,sep="")
  if(regim > 1 & is.null(setar)) datos <- data.frame(datos,t.series) else datos <- data.frame(datos)
  if(max(ars$q)>0){
    colnames(ex.series) <- paste("X",1:r,sep="")
    datos <- data.frame(datos,ex.series)
  }
  return(invisible(datos))
}

#' @method plot mtar
#' @export
plot.mtar <- function(x,...){
  temp <- residuals(x)
  dev.new()
  par(mfrow=c(1,3))
  qqnorm(temp$full,col="blue",main="Quantile-type residuals")
  abline(0,1,lty=3)
  acf(temp$full,main="Quantile-type residuals")
  pacf(temp$full,main="Quantile-type residuals")
}

#' @method residuals mtar
#' @export
residuals.mtar <- function(object,...,type=c("quantile","mahalanobis"),plot.it=FALSE,identify){
  type <- match.arg(type)
  n.sim <- object$n.sim
  dist <- object$dist
  if(dist %in% c("Skew-normal","Skew-Student-t"))
    stop("Models with Skew-normal or Skew-t as noise distribution are not supported!",call.=FALSE)
  k <- ncol(object$data[[1]]$y)
  resi <- matrix(NA,1,k)
  idsi <- vector()
  if(object$regim > 1) lims <- (object$ps+1):(length(object$threshold.series))
  if(object$dist %in% c("Student-t","Hyperbolic","Slash","Contaminated normal"))
    extra <- rowMeans(matrix(object$chains$extra,ncol=n.sim))
  ids <- 1:nrow(object$data[[1]]$X)
  if(object$regim > 1){
    h <- as.integer(round(mean(object$chains$h)))
    thresholds <- rowMeans(matrix(object$chains$thresholds,object$regim-1,ncol=n.sim))
    Z <- object$threshold.series[lims-h]
    regs <- cut(Z,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
  }else regs <- matrix(1,nrow(object$data[[1]]$y),1)
  for(i in 1:object$regim){
    location <- matrix(0,nrow(object$chains[[i]]$location),k)
    scale <- matrix(0,k,k)
    for(j in 1:k){
      location[,j] <- rowMeans(matrix(object$chains[[i]]$location[,seq(j,n.sim,by=k)],nrow(location),n.sim))
      scale[,j] <- rowMeans(matrix(object$chains[[i]]$scale[,seq(j,n.sim,by=k)],k,n.sim))
    }
    places <- regs==i
    y <- object$data[[i]]$y[places,]
    X <- object$data[[i]]$X[places,]
    resu <- y-X%*%location
    ei <- eigen(scale,symmetric=TRUE)
    resu <- resu%*%(ei$vectors%*%diag(1/sqrt(as.vector(ei$values)))%*%t(ei$vectors))
    resi <- rbind(resi,resu)
    idsi <- c(idsi,ids[places])
  }
  resi <- matrix(na.omit(resi),nrow(object$data[[1]]$X),k)
  TT <- 50000
  sim <- matrix(rnorm(TT*k),TT,k)
  if(dist=="Gaussian") u <- 1
  if(dist=="Student-t") u <- 1/rgamma(TT,shape=extra/2,rate=extra/2)
  if(dist=="Slash") u <- 1/rbeta(TT,shape1=extra/2,shape2=1)
  if(dist=="Contaminated normal") u <- 1/(1 - (1-extra[2])*rbinom(TT,1,extra[1]))
  if(dist=="Hyperbolic") u <- rgig(n=TT,lambda=1,chi=1,psi=extra^2)
  if(dist=="Laplace") u <- rexp(TT,rate=1/8)
  sim2 <- sort(rowSums(sim^2)*u)
  resi2 <- matrix(rowSums(resi^2),nrow(resi),1)
  resij <- resi2
  if(type=="quantile"){
    resi2 <- apply(resi2,1,function(x) mean(sim2<=x))
    resi3 <- ifelse(resi2>=0.5,1-resi2,resi2)
    resij <- qnorm(ifelse(.Machine$double.xmin>=resi3,.Machine$double.xmin,resi3))*ifelse(resi2>0.5,-1,1)
  }
  idsx <- sort(idsi,index=TRUE)$ix
  resij <- matrix(resij[idsx],length(resij),1)
  if(type=="quantile"){
    sim <- sim*matrix(sqrt(u),TT,k)
    for(i in 1:k){
      temp <- apply(matrix(resi[,i],nrow(resi),1),1,function(x) mean(sim[,i]<=x))
      temp2 <- ifelse(temp>=0.5,1-temp,temp)
      resi[,i] <- qnorm(ifelse(.Machine$double.xmin>=temp2,.Machine$double.xmin,temp2))*ifelse(temp>0.5,-1,1)
    }
  }
  resi <- matrix(resi[idsx,],nrow(resi),k)
  row.names(resij) <- row.names(resi) <- row.names(object$data[[1]]$y)
  colnames(resij) <- " "
  colnames(resi) <- colnames(object$data[[1]]$y)
  if(plot.it==TRUE){
    nano_ <- list(...)
    nano_$y <- resij
    nano_$type <- "p"
    if(is.null(nano_$pch)) nano_$pch <- 20
    if(is.null(nano_$ylim)) nano_$ylim <- 1.1*range(resij)
    if(is.null(nano_$pch)) nano_$pch <- 20
    if(is.null(nano_$col)) nano_$col <- "black"
    if(is.null(nano_$main)) nano_$main <- ""
    if(is.null(nano_$xlab)) nano_$xlab <- "Theoretical quantiles"
    if(is.null(nano_$ylab)) nano_$ylab <- "Sample quantiles"
    if(is.null(nano_$labels)) labels <- row.names(object$data[[1]]$y)
    else{
      labels <- nano_$labels
      nano_$labels <- NULL
    }
    dev.new()
    outm <- do.call("qqnorm",nano_)
    abline(0,1,lty=3)
    if(!missingArg(identify)) identify(outm$x,outm$y,n=max(1,floor(abs(identify))),labels=row.names(resij),labels=labels)
  }
  return(invisible(list(full=resij,by.component=resi)))
}
