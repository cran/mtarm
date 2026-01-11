#'
#' @title Forecasting for multivariate TAR models
#' @description Computes forecasts from a fitted multivariate Threshold Autoregressive (TAR) model.
#'
#' @param object An object of class \code{mtar} obtained from a call to \code{mtar()}.
#' @param ... Additional arguments that may affect the prediction method.
#' @param newdata An optional \code{data.frame} containing future values of the threshold
#' series (if included in the fitted model), the exogenous series (if included in the fitted
#' model), and, when \code{out.of.sample = TRUE}, the realized values of the output series.
#' @param n.ahead A positive integer specifying the number of steps ahead to forecast.
#' @param row.names An optional variable in \code{newdata} specifying labels for the time
#  points corresponding to the rows of \code{newdata}.
#' @param credible An optional numeric value in \eqn{(0,1)} specifying the level of the
#' required credible intervals. By default, \code{credible} is set to \code{0.95}.
#' @param out.of.sample An optional logical indicator. If \code{TRUE} then the log-score,
#' Absolute Error (AE), Absolute Percentage Error (APE) and Squared Error (SE) are computed
#' as measures of predictive accuracy. In this case, \code{newdata} must include the observed
#' values of the output series.
#'
#' @return A list containing the forecast summaries and, when requested, measures of predictive accuracy.
#'
#' \tabular{ll}{
#' \code{ypred}   \tab a matrix with the results of the forecasting,\cr
#' \tab \cr
#' \code{summary} \tab a matrix with the mean and credible intervals of the forecasting,\cr
#' }
#' @references Nieto, F.H. (2005) Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data.
#'             Communications in Statistics - Theory and Methods, 34, 905-930.
#' @references Romero, L.V. and Calderon, S.A. (2021) Bayesian estimation of a multivariate TAR model when the noise
#'             process follows a Student-t distribution. Communications in Statistics - Theory and Methods, 50, 2508-2530.
#' @references Calderon, S.A. and Nieto, F.H. (2017) Bayesian analysis of multivariate threshold autoregressive models
#'             with missing data. Communications in Statistics - Theory and Methods, 46, 296-318.
#' @references Karlsson, S. (2013) Chapter 15-Forecasting with Bayesian Vector Autoregression. In Elliott, G. and
#'             Timmermann, A. Handbook of Economic Forecasting, Volume 2, 791–89, Elsevier.
#' @references Vanegas, L.H. and Calderón, S.A. and Rondón, L.M. (2025) Bayesian estimation of a multivariate tar model when the
#'             noise process distribution belongs to the class of gaussian variance mixtures. International Journal
#'             of Forecasting.
#' @examples
#' \donttest{
#' ###### Example 1: Returns of the closing prices of three financial indexes
#' data(returns)
#' fit1 <- mtar(~ COLCAP + BOVESPA | SP500, data=returns, row.names=Date,
#'              subset={Date<="2016-03-14"}, dist="Student-t",
#'              ars=ars(nregim=3,p=c(1,1,2)), n.burnin=2000, n.sim=3000,
#'              n.thin=2, ssvs=TRUE)
#' out1 <- predict(fit1, newdata=subset(returns,Date>"2016-03-14"), n.ahead=10)
#' out1$summary
#'
#' ###### Example 2: Rainfall and two river flows in Colombia
#' data(riverflows)
#' fit2 <- mtar(~ Bedon + LaPlata | Rainfall, data=riverflows, row.names=Date,
#'              subset={Date<="2009-04-04"}, dist="Laplace",
#'              ars=ars(nregim=3,p=5), n.burnin=2000, n.sim=3000, n.thin=2)
#' out2 <- predict(fit2, newdata=subset(riverflows,Date>"2009-04-04"), n.ahead=10)
#' out2$summary
#'
#' ###### Example 3: Temperature, precipitation, and two river flows in Iceland
#' data(iceland.rf)
#' fit3 <- mtar(~ Jokulsa + Vatnsdalsa | Temperature | Precipitation,
#'              data=iceland.rf, subset={Date<="1974-12-21"}, row.names=Date,
#'              ars=ars(nregim=2,p=15,q=4,d=2), n.burnin=2000, n.sim=3000,
#'              n.thin=2)
#' out3 <- predict(fit3, newdata=subset(iceland.rf,Date>"1974-12-21"), n.ahead=10)
#' out3$summary
#' }
#' @method predict mtar
#' @export
#'
predict.mtar <- function(object,...,newdata,n.ahead=1,row.names,credible=0.95,out.of.sample=FALSE){
  # Number of regimes in the fitted model
  regim <- object$regim
  # Check if newdata is required (TAR models or models with exogenous variables)
  if(((regim>1 & is.null(object$setar)) | max(object$ars$q) > 0) & missing(newdata))
     stop("The argument 'newdata' is required!",call.=FALSE)
  # Check validity of forecast horizon
  if(n.ahead<=0 | n.ahead!=floor(n.ahead))
     stop("The value of the argument 'n.ahead' must be a positive integer!",call.=FALSE)
  # Check validity of credible level
  if(credible <=0 | credible >=1)
     stop("The value of the argument 'credible' must be within the interval (0,1)!",call.=FALSE)
  # newdata is required if out-of-sample predictive accuracy is requested
  if(out.of.sample & missing(newdata))
     stop("The argument 'new data' is required for out-of-sample forecast accuracy measures!",call.=FALSE)
  # Helper function to generate innovations from the specified distribution
  gendist <- function(dist,Sigma,extra,delta){
    Sigma <- as.matrix(Sigma)
    nor <- crossprod(chol(Sigma),matrix(rnorm(k),k,1))
    u <- switch(dist,
                "Gaussian"=,"Skew-normal"={1},
                "Student-t"=,"Skew-Student-t"={1/rgamma(1,shape=extra/2,rate=extra/2)},
                "Slash"={1/rbeta(1,shape1=extra/2,shape2=1)},
                "Contaminated normal"={1/(1 - (1-extra[2])*rbinom(1,1,extra[1]))},
                "Hyperbolic"={rgig(n=1,lambda=1,chi=1,psi=extra^2)},
                "Laplace"={rexp(1,rate=1/8)})
    if(dist %in% c("Skew-normal","Skew-Student-t"))
       out_ <- sqrt(u)*matrix(delta*qnorm(0.5 + runif(k)*0.5) + nor,k,1)
    else out_ <- nor*sqrt(u)
    return(t(out_))
  }
  # Build model frame for newdata if provided
  if(!missing(newdata)){
     mmf <- match.call(expand.dots = FALSE)
     m <- match(c("newdata","row.names"), names(mmf), 0)
     mmf <- mmf[c(1,m)]
     mmf$drop.unused.levels <- TRUE
     names(mmf)[names(mmf)=="newdata"] <- "data"
     mmf[[1]] <- as.name("model.frame")
     mmf <- eval(mmf, parent.frame())
     if(!missingArg(row.names)) row.names <- as.vector(as.character(model.extract(mmf,row.names)))
  }
  # Build threshold variable for multi-regime models
  if(regim>1 & is.null(object$setar)){
     mz <- model.part(Formula(object$formula), data = mmf, rhs = 2, terms = TRUE)
     Z <- model.matrix(mz, data = mmf)
     if(attr(terms(mz),"intercept")){
        Znames <- colnames(Z)
        Z <- as.matrix(Z[,-1])
     }
     if(nrow(Z)<n.ahead)
        stop(paste0("The data.frame 'newdata' must have at least ",n.ahead," rows!"),call.=FALSE)
     Z <- rbind(matrix(object$threshold.series[-c(1:object$ps),],nrow(object$threshold.series)-object$ps,1),matrix(Z[1:n.ahead,],n.ahead,1))
  }
  # Build exogenous regressors if present
  if(max(object$ars$q) > 0){
     mx2 <- model.part(Formula(object$formula), data = mmf, rhs = 3, terms = TRUE)
     X2 <- model.matrix(mx2, data = mmf)
     if(attr(terms(mx2),"intercept")){
        X2names <- colnames(X2)
        X2 <- as.matrix(X2[,-1])
        colnames(X2) <- X2names[-1]
     }
     if(nrow(X2)<n.ahead)
        stop(paste0("The data.frame 'newdata' must have at least ",n.ahead," rows!"),call.=FALSE)
     X2 <- rbind(matrix(object$exogenous.series[-c(1:object$ps),],nrow(object$exogenous.series)-object$ps,ncol(object$exogenous.series)),matrix(X2[1:n.ahead,],n.ahead,ncol(object$exogenous.series)))
  }
  # Extract observed output series and its dimensions
  y <- object$data[[1]]$y
  name <- colnames(y)
  n <- nrow(y)
  k <- ncol(y)
  # Storage matrix for simulated paths
  ysim <- rbind(matrix(y,n,k*object$n.sim),matrix(0,n.ahead,k*object$n.sim))
  # Define likelihood for out-of-sample predictive evaluation
  if(out.of.sample){
     Lik <- function(dist,y,M,Sigma,delta,log,nu){
       Sigma <- as.matrix(Sigma)
       resu <- matrix(y-M,1,k)
       if(dist %in% c("Skew-normal","Skew-Student-t")){
          A <- chol2inv(chol(Sigma + as.vector(delta^2)*diag(k)))
          out <- resu%*%tcrossprod(A,resu)
          muv <- resu%*%(A*matrix(delta,k,k,byrow=TRUE))
          Sigmav <- diag(k) - matrix(delta,k,k)*A*matrix(delta,k,k,byrow=TRUE)
          if(dist=="Skew-normal"){
             out <- out + k*log(2*pi)
             if(k > 1) out <- out - 2*log(max(.Machine$double.xmin,pmvnorm(lower=-as.numeric(muv),upper=rep(Inf,k),sigma=Sigmav)))
             else out <- out - 2*log(max(.Machine$double.xmin,pnorm(-muv/sqrt(Sigmav),lower.tail=FALSE)))
          }
          if(dist=="Skew-Student-t"){
             out <- (nu+k)*log(1 + out/nu) - 2*lgamma((nu+k)/2) + 2*lgamma(nu/2) + k*log(nu*pi)
             muv <- muv*matrix(sqrt((nu + k)/(nu + out)),nrow(muv),k)
             if(k > 1) out <- out - 2*log(max(.Machine$double.xmin,pmvt(lower=-as.numeric(muv),upper=rep(Inf,k),sigma=Sigmav,df=round(nu,0)+k)))
             else out <- out - 2*log(max(.Machine$double.xmin,pt(-muv/sqrt(Sigmav),df=nu,lower.tail=FALSE)))
          }
       }else{
        out <- resu%*%tcrossprod(chol2inv(chol(Sigma)),resu)
        out <- switch(dist,
                      "Gaussian"={out},
                      "Laplace"={-2*log(besselK(sqrt(out)/2,(2-k)/2)) - ((2-k)/2)*log(out)	+ (k+2)*log(2)},
                      "Student-t"={(nu+k)*log(1 + out/nu) - 2*lgamma((nu+k)/2) + 2*lgamma(nu/2) + k*log(nu/2)},
                      "Contaminated normal"={-2*log(nu[1]*exp(-out*nu[2]/2)*nu[2]^(k/2) + (1-nu[1])*exp(-out/2))},
                      "Slash"={ifelse(out==0,-2*log(nu/(nu+k)),-2*lgamma((k+nu)/2) + (k+nu)*log(out/2) -2*log((nu/2)*pgamma(1,shape=(k+nu)/2,rate=out/2)))},
                      "Hyperbolic"={-2*log(besselK(nu*sqrt(1+out),(2-k)/2)) -((2-k)/2)*log(1+out) - k*log(nu) +2*log(besselK(nu,1))})
        out <- out + k*log(2*pi)
       }
       out <- -(out + log(det(Sigma)))/2
       if(log) out <- out + sum(y)
       return(exp(out))
     }
    # Build true future observations for evaluation
    mx <- model.part(Formula(object$formula), data = mmf, rhs = 1, terms = TRUE)
    D <- model.matrix(mx, data = mmf)
    if(attr(terms(mx),"intercept")){
       Dnames <- colnames(D)
       D <- matrix(D[,-1],ncol=length(Dnames)-1)
       colnames(D) <- Dnames[-1]
    }
    ytrue <- D
    if(nrow(ytrue)<n.ahead)
       stop(paste0("The data.frame 'newdata' must have at least ",n.ahead," rows!"),call.=FALSE)
    else ytrue <- matrix(D[1:n.ahead,],n.ahead,k)
    if(object$log) ytrue <- log(ytrue)
    prek <- matrix(0,n.ahead,object$n.sim)
    ytrue <- rbind(matrix(y,n,k),matrix(ytrue,nrow(ytrue),k))
  }
  # Build deterministic regressors (time trend and seasonality)
  t.ids <- matrix(seq(n+1,n+n.ahead),n.ahead,1)
  switch(object$trend,
         "none"={Xd <- matrix(1,n.ahead,1)},
         "linear"={Xd <- matrix(cbind(1,t.ids),n,2)},
         "quadratic"={Xd <- matrix(cbind(1,t.ids,t.ids^2),n,3)})
  if(!is.null(object$nseason)){
     Itil <- matrix(diag(object$nseason)[,-1],object$nseason,object$nseason-1)
     Xd <- cbind(Xd,t(matrix(apply(t.ids,1,function(x) Itil[x%%object$nseason + 1,]),object$nseason-1,n.ahead)))
  }
  if(!object$Intercept) Xd <- Xd[,-1]
  # Main simulation loop over posterior draws
  for(i in 1:object$n.sim){
      if(regim > 1){
         h <- object$chains$h[i]
         thresholds <- object$chains$thresholds[,i]
      }
      for(j in (n+1):(n+n.ahead)){
          if(regim > 1){
             if(!is.null(object$setar)) iz <- ysim[j-h,(i-1)*k + object$setar] else iz <- Z[j-h]
             regs <- cut(iz,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
          }else regs <- 1
          X <- Xd[j-n,]
          for(l in 1:object$ars$p[regs]) X <- c(X,ysim[j-l,((i-1)*k+1):(i*k)])
          if(object$ars$q[regs] > 0) for(l in 1:object$ars$q[regs]) X <- c(X,X2[j-l,])
          if(object$ars$d[regs] > 0) for(l in 1:object$ars$d[regs]) X <- c(X,Z[j-l,])
          if(object$ssvs) zetai <- object$chains[[regs]]$zeta[,i] else zetai <- rep(1,length(X))
          M <- matrix(X[zetai==1],1,sum(zetai))%*%object$chains[[regs]]$location[zetai==1,((i-1)*k+1):(i*k)]
          if(object$dist %in% c("Gaussian","Laplace")) extrai <- NULL else extrai <- object$chains$extra[,i]
          if(object$dist %in% c("Skew-normal","Skew-Student-t")) deltai <- object$chains$delta[,i] else deltai <- rep(0,k)
          ysim[j,((i-1)*k+1):(i*k)] <- M + gendist(object$dist,object$chains[[regs]]$scale[,((i-1)*k+1):(i*k)],extrai,deltai)
          if(out.of.sample){
             prek[j-n,i] <- Lik(object$dist,ytrue[j,],M,object$chains[[regs]]$scale[,((i-1)*k+1):(i*k)],deltai,object$log,extrai)
          }
      }
  }
  # Remove in-sample observations and reshape simulations
  ysim <- matrix(ysim[-c(1:n),],n.ahead,k*object$n.sim)
  colnames(ysim) <- rep(colnames(y),object$n.sim)
  if(object$log) ysim <- exp(ysim)
  # Compute point forecasts and highest density intervals
  out_ <- vector()
  predi <- function(x){
     x <- matrix(x,1,length(x))
     ks <- seq(credible,1,length=object$n.sim*(1-credible))
     lis <- t(apply(x,1,quantile,probs=ks-credible))
     lss <- t(apply(x,1,quantile,probs=ks))
     dif <- apply(abs(lss-lis),1,which.min)
     out_ <- c(ifelse(object$log,median(x),mean(x)),lis[dif],lss[dif])
     names(out_) <- paste0(name[i],c(".Mean",".HDI_Low",".HDI_high"))
     return(out_)
  }
  # Apply summary function to each series
  for(i in 1:k) out_ <- cbind(out_,t(apply(matrix(ysim[,seq(i,k*object$n.sim,k)],n.ahead,object$n.sim),1,predi)))
  if(missingArg(row.names)) row.names <- 1:n.ahead
  rownames(ysim) <- rownames(out_) <- row.names[1:n.ahead]
  # Build output object
  out__ <- list(summary=out_,output.series=object$data[[1]]$y,rownames=!is.null(object$call$row.names))
  class(out__) <- "predict_mtar"
  # Add out-of-sample accuracy measures if requested
  if(out.of.sample){
     out__$log.score <- matrix(-log(rowMeans(prek)),n.ahead,1)
     if(object$log) ytrue <- exp(ytrue)
     ytrue <- matrix(ytrue[-c(1:n),],n.ahead,k)
     out__$AE <- matrix(abs(out_[,seq(1,3*k,3)] - ytrue),nrow(ytrue),k)
     out__$APE <- 100*matrix(abs((out_[,seq(1,3*k,3)] - ytrue)/ytrue),nrow(ytrue),k)
     out__$SE <- matrix((out_[,seq(1,3*k,3)] - ytrue)^2,nrow(ytrue),k)
     rownames(out__$AE) <- rownames(out__$log.score) <- rownames(out__$APE) <- rownames(out__$SE) <- rownames(out_)
     colnames(out__$log.score) <- "log.score"
     colnames(out__$AE) <- paste0(colnames(y),sep=".","AE")
     colnames(out__$APE) <- paste0(colnames(y),sep=".","APE")
     colnames(out__$SE) <- paste0(colnames(y),sep=".","SE")
  }
  # Return prediction object
  return(out__)
}

#' @method plot predict_mtar
#' @export
plot.predict_mtar <- function(x,...,last,historical=list(),forecasts=list(),forecasts.PI=list(),main=NULL){
    # Number of components of output series
    k <- ncol(x$output.series)
    n <- nrow(x$output.series)
    # Extract prediction summary (mean and credible intervals)
    out_ <- x$summary
    # Forecast horizon
    n.ahead <- nrow(out_)
    # Original output series
    y <- matrix(x$output.series,n,k)
    # Number of historical observations to display
    if(missingArg(last)) last <- n
    else{
       if(floor(last)!=last | last<=0 | last>n) stop(paste0("Argument 'last' must be a positive integer lower than or equal to ",n),call.=FALSE)
    }
    # Combine last n observations with forecasted values
    y2 <- rbind(matrix(y[(n-last+1):n,],ncol=k),matrix(out_[,seq(1,3*k,3)],ncol=k))
    # Define x-axis limits for historical data and forecasts
    forecasts$xlim <- historical$xlim <- c(n-last+1,n+n.ahead)
    # X positions for forecast points
    historical$x <- (n-last+1):n
    forecasts$x <- (n+1):(n+n.ahead)
    forecasts$xaxt <- historical$xaxt <- "n"
    if(is.null(forecasts$main)) main <- colnames(x$output.series)
    else{
      if(length(forecasts$main)<k) stop(paste0("Argument 'main' has an incorrect length.It must be of length ",k," but an object of length ",length(fitted$main)," was provided!"),call.=FALSE)
      else main <- forecasts$main
    }
    if(is.null(forecasts$ylab)) ylab <- rep("",k)
    else{
      if(length(forecasts$ylab)<k) stop(paste0("Argument 'ylab' has an incorrect length. It must be of length ",k," but an object of length ",length(fitted$ylab)," was provided!"),call.=FALSE)
      else ylab <- forecasts$ylab
    }
    # Default axis labels
    forecasts$xlab <- forecasts.PI$ylab <- forecasts.PI$xlab <- historical$xlab <- historical$ylab <- ""
    # Default graphical parameters for historical series
    if(is.null(historical$type)) historical$type <- "b"
    if(is.null(historical$pch)) historical$pch <- 20
    if(is.null(historical$lty)) historical$lty <- 1
    if(is.null(historical$col)) historical$col <- "black"
    # Default graphical parameters for forecasts
    if(is.null(forecasts$type)) forecasts$type <- "b"
    if(is.null(forecasts$pch)) forecasts$pch <- 20
    if(is.null(forecasts$lty)) forecasts$lty <- 1
    if(is.null(forecasts$col)) forecasts$col <- "blue"
    # Default graphical parameters for prediction intervals
    if(is.null(forecasts.PI$density)) forecasts.PI$density <- NA
    if(is.null(forecasts.PI$col)) forecasts.PI$col <- "light gray"
    # Loop over each component of the output series
    for(i in 1:k){
        dev.new()
        # Set common y-axis limits based on observed data and forecasts
        historical$ylim <- forecasts.PI$ylim <- forecasts$ylim <- range(y2[,i],out_[,1:3+3*(i-1)])
        # Historical data for the i-th component of the output series
        historical$y <- y[(n-last+1):n,i]
        # Plot historical series
        do.call("plot",historical)
        # Coordinates for the prediction interval polygon
        xs <- c((n+1):(n+n.ahead),(n+n.ahead):(n+1))
        ys <- c(out_[,2+3*(i-1)],out_[n.ahead:1,3+3*(i-1)])
        # Add prediction interval polygon
        par(new=TRUE)
        forecasts.PI$x <- xs
        forecasts.PI$y <- ys
        do.call("polygon",forecasts.PI)
        # Set plot title
        forecasts$main <- main[i]
        forecasts$ylab <- ylab[i]
        # Add forecasted values line
        forecasts$y <- out_[,1+3*(i-1)]
        par(new=TRUE)
        do.call("plot",forecasts)
    }
}
#'
#'
#' @method out.of.sample mtar
#' @export
out.of.sample.mtar <- function(...,newdata,n.ahead=1,by.component=FALSE){
  # Collect all fitted mtar models passed through ...
  another <- list(...)
  # Store the original function call
  call. <- match.call()
  # Number of components of the output series
  k <- ncol(another[[1]]$data[[1]]$y)
  # Initialize output matrix:
  # rows = number of models, columns depend on whether results are by component
  out <- matrix(0,length(another),ifelse(by.component,3*k+1,4))
  # Vector to store row names (model identifiers)
  outnames <- vector()
  # Loop over each model
  for(ii in 1:length(another)){
      # Obtain out-of-sample predictions and accuracy measures
      temp <- predict(another[[ii]],newdata=newdata,n.ahead=n.ahead,out.of.sample=TRUE)
      # Average log predictive score
      out[ii,1] <- mean(temp$log.score)
      if(by.component){
         # Mean Absolute Error by component
         out[ii,2:(k+1)] <- apply(temp$AE,2,mean)
         # Mean Absolute Percentage Error by component
         out[ii,(k+2):(2*k+1)] <- apply(temp$APE,2,mean)
         # Mean Squared Error by component
         out[ii,(2*k+2):(3*k+1)] <- apply(temp$SE,2,mean)
      }else{
         # Overall Mean Absolute Error
         out[ii,2] <- mean(apply(temp$AE,1,mean))
         # Overall Mean Absolute Percentage Error
         out[ii,3] <- mean(apply(temp$APE,1,mean))
         # Overall Mean Squared Error
         out[ii,4] <- mean(apply(temp$SE,1,mean))
      }
      # Use the model call as the row name
      outnames[ii] <- as.character(call.[ii+1])
  }
  # Assign row names to the output matrix
  rownames(out) <- outnames
  # Assign column names depending on the output format
  if(!by.component) colnames(out) <- c("log.score","MAE","MAPE","MSE")
  else{
    name <- c("log.score")
    me <- c("MAE","MAPE","MSE")
    for(i in 1:length(me)) name <- c(name,paste0(colnames(another[[1]]$data[[1]]$y),sep=".",me[i]))
    colnames(out) <- name
  }
  # Return the out-of-sample evaluation results
  return(out)
}
#'
#' @method out.of.sample listmtar
#' @export
out.of.sample.listmtar <- function(x,newdata,n.ahead=1,by.component=FALSE,...){
  # Number of components of the output series
  k <- ncol(x[[1]]$data[[1]]$y)
  # Initialize output matrix:
  # rows = number of models, columns depend on whether results are by component
  out <- matrix(0,length(x),ifelse(by.component,3*k+1,4))
  # Loop over each object of class mtar in the list
  for(i in 1:length(x)){
      # Obtain out-of-sample predictions and accuracy measures
      temp <- predict(x[[i]],newdata=newdata,n.ahead=n.ahead,out.of.sample=TRUE)
      # Average log predictive score
      out[i,1] <- mean(temp$log.score)
      if(by.component){
         # Mean Absolute Error by component (averaged over horizons)
         out[i,2:(k+1)] <- mean(apply(temp$AE,1,mean))
         # Mean Absolute Percentage Error by component (averaged over horizons)
         out[i,(k+2):(2*k+1)] <- mean(apply(temp$APE,1,mean))
         # Mean Squared Error by component (averaged over horizons)
         out[i,(2*k+2):(3*k+1)] <- mean(apply(temp$SE,1,mean))
      }else{
         # Overall Mean Absolute Error
         out[i,2] <- mean(apply(temp$AE,1,mean))
         # Overall Mean Absolute Percentage Error
         out[i,3] <- mean(apply(temp$APE,1,mean))
         # Overall Mean Squared Error
         out[i,4] <- mean(apply(temp$SE,1,mean))
      }
  }
  # Assign row names using the names of the model list
  rownames(out) <- names(x)
  # Assign column names depending on the output format
  if(!by.component) colnames(out) <- c("log.score","MAE","MAPE","MSE")
  else{
    name <- c("log.score")
    me <- c("MAE","MAPE","MSE")
    for(i in 1:length(me)) name <- c(name,paste0(colnames(x[[1]]$data[[1]]$y),sep=".",me[i]))
    colnames(out) <- name
  }
  # Return the out-of-sample evaluation results
  return(out)
}
#'
#'
