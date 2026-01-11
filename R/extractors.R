#'
#' @title coef method for objects of class \code{mtar}
#' @param object an object of class \code{mtar} obtained from a call to the \code{mtar()} function.
#' @param FUN a function to be applied to the MCMC chains associated with each model parameter.
#' By default, \code{FUN} is set to \code{mean}.
#' @param ... additional arguments passed to \code{FUN}.
#' @return A list containing the summary statistics obtained by applying \code{FUN} to the
#' MCMC chains of each model parameter.
#' @method coef mtar
#' @export
coef.mtar <- function(object,...,FUN=mean){
  # Number of endogenous variables
  k <- ncol(object$data[[1]]$y)
  # Number of MCMC simulations
  n.sim <- object$n.sim
  # Output list to store results
  out_ <- list()
  # Loop over regimes
  for(i in 1:object$regim){
      # Initialize regime-specific container
      out_[[i]] <- list()
      # Initialize matrices for location and scale summaries
      out <- outs <- vector()
      # Loop over variables
      for(j in 1:k){
          # Extract MCMC draws of location parameters for variable j
          temp <- object$chains[[i]]$location[,seq(j,n.sim*k,k)]
          # Extract MCMC draws of scale parameters
          temps <- matrix(matrix(object$chains[[i]]$scale,k,n.sim*k)[,seq(j,n.sim*k,k)],nrow=k)
          # Apply summarizing function to location draws
          out <- cbind(out,apply(temp,1,FUN,...))
          # Apply summarizing function to scale draws
          outs <- cbind(outs,apply(temps,1,FUN,...))
      }
      # If SSVS procedure is enabled
      if(object$ssvs){
         # Base indices for OS lags
         quien <- object$deterministic + seq(k,object$ars$p[i]*k,k)
         namezeta <- paste0("OS.lag(",1:object$ars$p[i],")")
         # Add ES lags if present
         if(object$ars$q[i]>0){
            quien <- c(quien,max(quien)+seq(object$r,object$ars$q[i]*object$r,object$r))
            namezeta <- c(namezeta,paste0("ES.lag(",1:object$ars$q[i],")"))
         }
         # Add TS lags if present
         if(object$ars$d[i]>0){
            quien <- c(quien,max(quien)+seq(1,object$ars$d[i],1))
            namezeta <- c(namezeta,paste0("TS.lag(",1:object$ars$d[i],")"))
         }
         # Summaries for selected zeta parameters
         temp <- matrix(apply(matrix(object$chains[[i]]$zeta[quien,],length(quien),n.sim),1,FUN,...),1,length(quien))
         # Assign names to zeta output
         colnames(temp) <- namezeta
         rownames(temp) <- "zeta"
         # Store zeta summaries
         out_[[i]]$zeta <- temp
      }
      # Store summarized location and scale parameters
      out_[[i]]$location <- out
      out_[[i]]$scale <- outs
      # Add column / row names based on variable names
      colnames(out_[[i]]$location) <- colnames(out_[[i]]$scale) <- rownames(out_[[i]]$scale) <- colnames(object$data[[1]]$y)
  }
  if(object$regim > 1){
     # Summaries of delay parameter
     out_$delay <- apply(matrix(object$chains$h,1,n.sim),1,FUN,...)
     # Summaries of threshold parameters
     out_$thresholds <- matrix(apply(object$chains$thresholds,1,FUN,...),object$regim-1,1)
     rownames(out_$thresholds) <- paste0("Threshold",1:(object$regim-1))
     colnames(out_$thresholds) <- ""
  }
  # Skewness parameters for skewed distributions
  if(object$dist %in% c("Skew-normal","Skew-Student-t")){
     out_$skewness <- matrix(apply(object$chains$delta,1,FUN,...),k,1)
     rownames(out_$skewness) <- paste0("delta",1:k)
     colnames(out_$skewness) <- ""
  }
  # Extra tail-thickness parameters for heavy-tailed distributions
  if(object$dist %in% c("Slash","Contaminated normal","Student-t","Hyperbolic","Skew-Student-t")){
     out_$extra <- matrix(apply(object$chains$extra,1,FUN,...),nrow(object$chains$extra),1)
     rownames(out_$extra) <- paste0("nu",1:nrow(object$chains$extra))
     colnames(out_$extra) <- ""
  }
  # Return all summarized outputs
  return(out_)
}

#' @title vcov method for objects of class \code{mtar}
#' @description Computes estimates of the variance–covariance matrices for the scale
#' parameters of a fitted multivariate TAR model.
#' @param object an object of class \code{mtar}, typically obtained from a call to
#' \code{mtar()}.
#' @param FUN a function to be applied to the MCMC chains of the scale parameters in order
#' to obtain point estimates. By default, \code{FUN} is set to \code{mean()}.
#' @param ... additional arguments passed to \code{FUN}.
#' @return A list containing the variance–covariance estimates obtained by applying
#'         \code{FUN} to the MCMC chains associated with the scale parameters.
#' @method vcov mtar
#' @export
vcov.mtar <- function(object,...,FUN=mean){
  # Number of endogenous variables
  k <- ncol(object$data[[1]]$y)
  # Number of MCMC simulations
  n.sim <- object$n.sim
  # List to store variance–covariance matrices by regime
  out_ <- list()
  # Loop over regimes
  for(i in 1:object$regim){
      # Placeholder for the covariance estimates
      outs <- vector()
      # Loop over variables to extract the relevant MCMC draws
      for(j in 1:k){
          temps <- matrix(matrix(object$chains[[i]]$scale,k,n.sim*k)[,seq(j,n.sim*k,k)],nrow=k)
          # Apply summary function (mean by default) to each row
          outs <- cbind(outs,apply(temps,1,FUN,...))
      }
      # Apply summary function (mean by default) to each row
      out_[[i]] <- outs
      # Assign row/column names corresponding to variable names
      colnames(out_[[i]]) <- rownames(out_[[i]]) <- colnames(object$data[[1]]$y)
  }
  return(out_)
}


#' @method summary mtar
#' @export
summary.mtar <- function(object, credible=0.95,...){
  # Number of endogenous variables
  k <- ncol(object$data[[1]]$y)
  # Number of MCMC simulations
  n.sim <- object$n.sim
  # Output list
  out_ <- list()
  # Internal helper function to compute summary statistics for MCMC draws
  resumen <- function(x){
        x <- matrix(x,ifelse(is.null(nrow(x)),1,nrow(x)),ifelse(is.null(ncol(x)),length(x),ncol(x)))
        # Initialize output: mean, 2(1-PD), HDI low, HDI high
        y <- matrix(0,nrow(x),4)
        # Posterior means
        y[,1] <- rowMeans(x)
        # Compute 2(1-PD): twice one minus posterior probability of direction
        y[,2] <- apply(x,1,function(x) min(mean(sign(median(x))*x > 0),1-1/200000))
        y[,2] <- 2*(1 - y[,2])
        # Build high-density interval (HDI) by scanning quantile pairs
        ks <- seq(credible,1,length=n.sim*(1-credible))
        lis <- t(apply(x,1,quantile,probs=ks-credible))
        lss <- t(apply(x,1,quantile,probs=ks))
        # Pick smallest-width interval
        dif <- apply(abs(lss-lis),1,which.min)
        y[,3] <- lis[cbind(1:nrow(x),dif)]
        y[,4] <- lss[cbind(1:nrow(x),dif)]
        # Naming columns
        colnames(y) <- c("   Mean"," 2(1-PD) ","HDI_low","HDI_high")
        return(y)
  }
  # Threshold summary if more than one regime
  if(object$regim > 1){
     thresholds <- matrix(resumen(matrix(object$chains$thresholds,nrow=object$regim-1))[,c(1,3,4)],ncol=3)
     out_$h <- mean(object$chains$h)
     out_$thresholds <- thresholds
     out_$threshold.series <- object$ts
  }
  # Store names of endogenous series
  out_$output.series <- ifelse(length(colnames(object$data[[1]]$y))==1,colnames(object$data[[1]]$y),paste(colnames(object$data[[1]]$y),collapse="    |    "))
  # Exogenous series (if present)
  if(min(object$ars$q)>0) out_$exogenous.series <- ifelse(length(colnames(object$exogenous.series))==1,colnames(object$exogenous.series),paste(colnames(object$exogenous.series),collapse="    |    "))
  out_$dist <- object$dist
  out_$ssvs <- object$ssvs
  out_$location <- list()
  out_$scale <- list()
  out_$ars <- object$ars
  # Construct description of deterministic components
  a <- ifelse(object$Intercept,"Intercept","")
  b <- ifelse(object$trend=="none","",paste0(ifelse(a=="","a ","+ a "),object$trend," time trend"))
  c <- ifelse(is.null(object$nseason),"",paste0("+ ",object$nseason," seasonal periods"))
  out_$deterministics <- paste(a,b,c)
  # Sample size
  out_$sample.size <- nrow(object$data[[1]]$X)
  # Initialize SSVS storage (if used)
  if(object$ssvs) out_$zeta <- list()
  # Loop over regimes
  for(i in 1:object$regim){
      out <- outs <- vector()
      # Loop over endogenous variables
      for(j in 1:k){
          temp <- object$chains[[i]]$location[,seq(j,n.sim*k,k)]
          temps <- matrix(matrix(object$chains[[i]]$scale,k,n.sim*k)[,seq(j,n.sim*k,k)],nrow=k)
          # Add NA columns to separate blocks after the first variable
          if(j > 1){
             out <- cbind(out,matrix(NA,nrow(out),1),resumen(temp))
             outs <- cbind(outs,resumen(temps)[,c(1,3,4)])
          }else{
               out <- resumen(temp)
               outs <- resumen(temps)[,c(1,3,4)]
          }
      }
      # Reshape scale matrix
      outs <- matrix(outs,k,3*k)
      rownames(out) <- object$name[[i]]
      # Drop coefficients excluded by SSVS
      if(object$ssvs){
         temp <- apply(object$chains[[i]]$zeta,1,mean)
         out <- out[temp>0.5,]
      }
      # Reorder scale matrix with separators
      outs <- matrix(cbind(outs,NA),k,3*k+1)
      outs <- outs[,c(seq(1,3*k,3),3*k+1,seq(2,3*k,3),3*k+1,seq(3,3*k,3))]
      outs <- matrix(outs,k,length(outs)/k)
      rownames(outs) <- colnames(object$data[[1]]$y)
      colnames(outs) <- c(rownames(outs),"",rownames(outs),"",rownames(outs))
      # SSVS summaries for lag selection
      if(object$ssvs){
         quien <- object$deterministic + seq(k,object$ars$p[i]*k,k)
         namezeta <- paste0("OS.lag(",1:object$ars$p[i],")")
         # Exogenous lags
         if(object$ars$q[i]>0){
            quien <- c(quien,max(quien)+seq(object$r,object$ars$q[i]*object$r,object$r))
            namezeta <- c(namezeta,paste0("ES.lag(",1:object$ars$q[i],")"))
         }
         # Threshold lags
         if(object$ars$d[i]>0){
            quien <- c(quien,max(quien)+seq(1,object$ars$d[i],1))
            namezeta <- c(namezeta,paste0("TS.lag(",1:object$ars$d[i],")"))
         }
         # Compute mean posterior inclusion probabilities
         temp <- matrix(apply(matrix(object$chains[[i]]$zeta[quien,],length(quien),n.sim),1,mean),1,length(quien))
         colnames(temp) <- namezeta
         rownames(temp) <- "SSVS"
      }
      # Store results for regime i
      out_$location[[i]] <- out
      out_$scale[[i]] <- outs
      if(object$ssvs) out_$zeta[[i]] <- temp
  }
  # Skewness parameter summaries
  if(object$dist %in% c("Skew-normal","Skew-Student-t")){
     out <- resumen(object$chains$delta)
     rownames(out) <- paste0(paste0("delta",1:k),paste0(rep("",max(nchar(object$name[[1]]))-6),collapse=" "))
     out_$delta <- out
  }
  # Extra parameters for heavy-tailed distributions
  if(object$dist %in% c("Slash","Contaminated normal","Student-t","Hyperbolic","Skew-Student-t")){
     out <- resumen(object$chains$extra)
     out[,2] <- NA
     # Label rows appropriately depending on distribution
     if(object$dist %in% c("Slash","Student-t","Hyperbolic","Skew-Student-t")) rownames(out)[nrow(out)] <- paste0("nu",paste0(rep("",max(nchar(object$name[[1]]))-1),collapse=" "))
     else rownames(out)[nrow(out):(nrow(out)-1)] <- paste0(c("nu2","nu1"),paste0(rep("",max(nchar(object$name[[1]]))-2),collapse=" "))
     out_$extra <- out
  }
  # Attach row names if provided by user
  if(!is.null(object$row.names)) out_$row.names <- object$row.names
  # Assign summary class
  class(out_) <- "summary_mtar"
  return(out_)
}
#'
#'
#' @method print summary_mtar
#' @export
print.summary_mtar <- function(x, digits=max(3, getOption("digits") - 2),...){
  # Print sample size and optional row names
  cat("\n\nSample size          :",paste0(x$sample.size," time points",ifelse(is.null(x$row.names),"",x$row.names)))
  # Print sample size and optional row names
  cat(ifelse(x$ssvs,"\nOutput Series (OS)   :","\nOutput Series        :"),x$output.series)
  # Print threshold series when more than one regime exists
  if(x$ars$nregim > 1) cat(ifelse(max(x$ars$d)>0,"\nThreshold Series (TS):","\nThreshold Series     :"),x$threshold.series)
  # Print exogenous series if present
  if(max(x$ars$q)>0) cat("\nExogenous Series (ES):",x$exogenous.series)
  # Print error distribution and number of regimes
  cat("\nError Distribution   :",x$dist)
  cat("\nNumber of regimes    :",x$ars$nregim)
  # Print deterministic components
  cat("\nDeterministics       :",x$deterministics)
  # Print autoregressive lag orders for each regime
  if(min(x$ars$p)==max(x$ars$p)) cat("\nAutoregressive orders:",paste0(x$ars$p[1]," in each regime"))
  else cat("\nAutoregressive orders:",paste(x$ars$p,collapse=", "))
  # Print maximum lags for exogenous variables if present
  if(max(x$ars$q)>0){
     if(min(x$ars$q)==max(x$ars$q)) cat("\nMaximum lags for ES  :",paste0(x$ars$q[1]," in each regime"))
     else cat("\nMaximum lags for ES  :",paste(x$ars$q,collapse=", "))
  }
  # Print maximum lags for threshold variable if present
  if(max(x$ars$d)>0){
     if(min(x$ars$d)==max(x$ars$d)) cat("\nMaximum lags for TS  :",paste0(x$ars$d[1]," in each regime"))
     else cat("\nMaximum lags for TS  :",paste(x$ars$d,collapse=", "))
  }
  cat("\n\n")
  # If multiple regimes, print the threshold intervals (Mean, HDI_low, HDI_high)
  if(x$ars$nregim > 1){
     # Round threshold summary statistics
     thresholds <- round(x$thresholds,digits=digits)
     # Construct interval strings for printing
     thresholds1 <- paste0(c("(-Inf",paste0("(",thresholds[,1])),",",c(paste0(thresholds[,1],"]"),"Inf)"))
     thresholds2 <- paste0(c("(-Inf",paste0("(",thresholds[,2])),",",c(paste0(thresholds[,2],"]"),"Inf)"))
     thresholds3 <- paste0(c("(-Inf",paste0("(",thresholds[,3])),",",c(paste0(thresholds[,3],"]"),"Inf)"))
     # Build a data frame of the printed intervals
     d <- data.frame(cbind(thresholds1,thresholds2,thresholds3))
     rownames(d) <- paste("Regime",1:nrow(d))
     colnames(d) <- rep(" ",3)
     # Print threshold summary table
     cat("\nThresholds (Mean, HDI_low, HDI_high)\n")
     print(d)
  }
  # Loop through regimes and print location and scale summaries
  for(i in 1:x$ars$nregim){
      cat("\n\nRegime",i,":\n")
      # Print SSVS inclusion probabilities for regime i if present
      if(x$ssvs) print(round(x$zeta[[i]],digits=2),digits=2)
      # Print autoregressive coefficients summary
      cat("\nAutoregressive coefficients\n")
      print(round(x$location[[i]],digits=digits),na.print="   |   ")
      # Print scale parameter summary
      cat("\nScale parameter (Mean, HDI_low, HDI_high)\n")
      print(round(x$scale[[i]],digits=digits),na.print="   .")
  }
  # Print skewness parameter if present
  if(x$dist %in% c("Skew-normal","Skew-Student-t")){
     cat("\n\nSkewness parameter","\n")
     print(round(x$delta,digits=digits),na.print="   .   ",)
  }
  # Print extra parameters if present
  if(x$dist %in% c("Slash","Contaminated normal","Student-t","Hyperbolic","Skew-Student-t")){
     cat("\n\nExtra parameter","\n")
     print(round(x$extra,digits=digits),na.print="   .   ")
  }
  cat("\n\n")
}
#'
#'
#' @method plot mtar
#' @export
plot.mtar <- function(x,...){
  plot(resid(x),...)
}
#'
#' @method fitted mtar
#' @export
fitted.mtar <- function(object,...){
  # Number of posterior simulation draws
  n.sim <- object$n.sim
  # Error distribution used in the model
  dist <- object$dist
  # Number of endogenous variables
  k <- ncol(object$data[[1]]$y)
  # Sample size
  n <- nrow(object$data[[1]]$X)
  Mi <- idsi <- vector()
  # Original row index vector
  ids <- 1:n
  a <- coef(object)
  # Determine regime membership for each observation (if more than one regime)
  if(object$regim > 1){
     lims <- (object$ps+1):(length(object$threshold.series))
     # Estimated delay
     h <- a$delay
     # Posterior mean of thresholds
     thresholds <- a$thresholds
     # Delayed threshold variable
     Z <- object$threshold.series[lims-h]
     # Assign each observation to a regime using threshold cut points
     regs <- cut(Z,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
  }else regs <- matrix(1,n,1) # Single regime: all observations belong to regime 1
  # Loop over regimes to compute fitted values
  for(i in 1:object$regim){
      # Logical index: which rows belong to regime i
      places <- regs==i
      y <- object$data[[i]]$y[places,]
      X <- object$data[[i]]$y[places,]
      location <- a[[i]]$location
      # Linear predictor X*beta for regime i
      if(object$ssvs){
         zeta <- rowMeans(object$chains[[i]]$zeta)
         X <- object$data[[i]]$X[places,zeta>0.5]
         location <- location[zeta>0.5,]
      }
      M <- X%*%location
      # Store original row indices for sorting later
      idsi <- c(idsi,ids[places])
      # Append fitted values
      Mi <- rbind(Mi,M)
  }
  # Reorder fitted values to the original data ordering
  idsx <- sort(idsi,index=TRUE)$ix
  Mi <- matrix(Mi[idsx,],n,k)
  # Add appropriate row and column names
  rownames(Mi) <- rownames(object$data[[1]]$y)
  colnames(Mi) <- colnames(object$data[[1]]$y)
  # Build output list with observed and fitted values
  out <- list(observed=object$data[[1]]$y,fitted=Mi)
  if(object$log){
     out$observed <- exp(out$observed)
     out$fitted <- exp(out$fitted)
  }
  out$rownames <- !is.null(object$call$row.names)
  # Assign S3 class
  class(out) <- "fittedmtar"
  return(out)
}
#'
#' @method plot fittedmtar
#' @export
plot.fittedmtar <- function(x,...,last=NULL,observed=list(),fitted=list()){
      k <- ncol(x$fitted)
      n <- nrow(x$fitted)
      if(missingArg(last)) last <- n
      else{
        if(floor(last)!=last | last<=0 | last>n) stop(paste0("Argument 'last' must be a positive integer lower than or equal to ",n),call.=FALSE)
      }
      if(x$rownames) observed$x <- as.Date(rownames(x$observed))[(n-last+1):n] else observed$x <- (n-last+1):n
      observed$main <- observed$ylab <- observed$xlab <- ""
      if(is.null(observed$col)) observed$col <- "black"
      if(is.null(observed$type)) observed$type <- "b"
      if(is.null(observed$pch)) observed$pch <- 20
      if(is.null(observed$lty)) observed$lty <- 3
      fitted$x <- observed$x
      if(is.null(fitted$col)) fitted$col <- "blue"
      if(is.null(fitted$type)) fitted$type <- "l"
      if(is.null(fitted$pch)) fitted$pch <- 20
      if(is.null(fitted$lty)) fitted$lty <- 3
      if(is.null(fitted$main)) fitted.main <- colnames(x$observed)
      else{
         if(length(fitted$main)<k) stop(paste0("Argument 'main' has an incorrect length.It must be of length ",k," but an object of length ",length(fitted$main)," was provided!"),call.=FALSE)
         else fitted.main <- fitted$main
      }
      if(is.null(fitted$ylab)) fitted.ylab <- rep("",k)
      else{
         if(length(fitted$ylab)<k) stop(paste0("Argument 'ylab' has an incorrect length. It must be of length ",k," but an object of length ",length(fitted$ylab)," was provided!"),call.=FALSE)
         else fitted.ylab <- fitted$ylab
      }
      if(is.null(fitted$xlab)) fitted.xlab <- rep("",k)
      else{
        if(length(fitted$xlab)<k) stop(paste0("Argument 'xlab' has an incorrect length. It must be of length ",k," but an object of length ",length(fitted$xlab)," was provided!"),call.=FALSE)
        else fitted.xlab <- fitted$xlab
      }
      leg <- list()
      leg$x <- "topleft"
      leg$legend <- c("Observed","Fitted")
      leg$lty <- leg$pch <- c(NA,NA)
      if(observed$type %in% c("p","b")) leg$pch[1] <- observed$pch
      if(observed$type %in% c("l","b")) leg$lty[1] <- observed$lty
      if(fitted$type %in% c("p","b")) leg$pch[2] <- fitted$pch
      if(fitted$type %in% c("l","b")) leg$lty[2] <- fitted$lty
      leg$col <- c(observed$col,fitted$col)
      leg$bty <- "n"
      for(j in 1:k){
          dev.new()
          observed$y <- x$observed[(n-last+1):n,j]
          observed$ylim <- range(cbind(x$observed[(n-last+1):n,j],x$fitted[(n-last+1):n,j]))
          observed$ylim[2] <- observed$ylim[2]*1.1
          fitted$y <- x$fitted[(n-last+1):n,j]
          fitted$ylim <- observed$ylim
          fitted$main <- fitted.main[j]
          fitted$ylab <- fitted.ylab[j]
          fitted$xlab <- fitted.xlab[j]
          do.call("plot",observed)
          par(new=TRUE)
          do.call("plot",fitted)
          do.call("legend",leg)
      }
}
#' @method print listmtar
#' @export
print.listmtar <- function(x,...){
  x2 <- x[[1]]
  # Build a string with the names of the output series
  # If there is only one series, keep its name, otherwise, concatenate all names
  output.series <- ifelse(length(colnames(x2$data[[1]]$y))==1,colnames(x2$data[[1]]$y),paste(colnames(x2$data[[1]]$y),collapse="    |    "))
  # Build a string with the names of the exogenous series, if present
  if(min(x2$ars$q)>0) exogenous.series <- ifelse(length(colnames(x2$exogenous.series))==1,colnames(x2$exogenous.series),paste(colnames(x2$exogenous.series),collapse="    |    "))
  # Print the sample size and optional row names information
  cat("\n\nSample size          :",paste0(nrow(x2$data[[1]]$X)," time points",ifelse(is.null(x2$row.names),"",x2$row.names)))
  # Print the output series names
  cat("\nOutput Series        :",output.series)
  # If there is more than one regime, print the threshold series information
  if(x2$ars$nregim > 1) cat(ifelse(max(x2$ars$d)>0,"\nThreshold Series (TS):","\nThreshold Series     :"),x2$ts)
  # If exogenous variables are included, print their names
  if(max(x2$ars$q)>0) cat("\nExogenous Series (ES):",exogenous.series)
  # Print the noise process distribution
  cat("\nError Distribution   :",paste0(x2$dist,collapse=", "))
  # Print the range of the number of regimes across models in the list
  cat("\nNumber of regimes    :",paste0(c(x[[1]]$ars$nregim,x[[length(x)]]$ars$nregim),collapse=" to "))
  # Identify the deterministic components in the model: intercept, trend, and seasonality
  a <- ifelse(x2$Intercept,"Intercept","")
  b <- ifelse(x2$trend=="none","",paste0(ifelse(a=="","a ","+ a "),x2$trend," time trend"))
  c <- ifelse(is.null(x2$nseason),"",paste0("+ ",x2$nseason," seasonal periods"))
  # Print the deterministic terms included in the model
  cat("\nDeterministics       :",paste(a,b,c))
  # Print the range of autoregressive orders p across models
  cat("\nAutoregressive orders:",paste0(c(x[[1]]$ars$p[1],x[[length(x)]]$ars$p[1]),collapse=" to "))
  # If exogenous variables are included, print the maximum lag orders q
  if(max(x2$ars$q)>0){
     cat("\nMaximum lags for ES  :",paste0(c(x[[1]]$ars$q[1],x[[length(x)]]$ars$q[1]),collapse=" to "))
  }
  # If threshold effects are included, print the maximum lag orders d
  if(max(x2$ars$d)>0){
     cat("\nMaximum lags for TS  :",paste0(c(x[[1]]$ars$d[1],x[[length(x)]]$ars$d[1]),collapse=" to "))
  }
}
#' @method print mtar
#' @export
print.mtar <- function(x,...){
    # Build a string with the names of the output series
    # If there is only one component, keep its name, otherwise, concatenate all names
    output.series <- ifelse(length(colnames(x$data[[1]]$y))==1,colnames(x$data[[1]]$y),paste(colnames(x$data[[1]]$y),collapse="    |    "))
    # Build a string with the names of the exogenous series, if present
    if(min(x$ars$q)>0) exogenous.series <- ifelse(length(colnames(x$exogenous.series))==1,colnames(x$exogenous.series),paste(colnames(x$exogenous.series),collapse="    |    "))
    # Print the sample size and optional row names information
    cat("\n\nSample size          :",paste0(nrow(x$data[[1]]$X)," time points",ifelse(is.null(x$row.names),"",x$row.names)))
    # Print the output series names
    cat("\nOutput Series        :",output.series)
    # If there is more than one regime, print the threshold series information
    if(x$ars$nregim > 1) cat(ifelse(max(x$ars$d)>0,"\nThreshold Series (TS):","\nThreshold Series     :"),x$ts)
    # If exogenous variables are included in the model, print their names
    if(max(x$ars$q)>0) cat("\nExogenous Series (ES):",exogenous.series)
    # Print the noise process distribution
    cat("\nError Distribution   :",x$dist)
    # Print the number of regimes in the model
    cat("\nNumber of regimes    :",x$regim)
    # Identify the deterministic components in the model: intercept, trend, and seasonality
    a <- ifelse(x$Intercept,"Intercept","")
    b <- ifelse(x$trend=="none","",paste0(ifelse(a=="","a ","+ a "),x$trend," time trend"))
    c <- ifelse(is.null(x$nseason),"",paste0("+ ",x$nseason," seasonal periods"))
    # Print the deterministic terms included in the model
    cat("\nDeterministics       :",paste(a,b,c))
    # Print the autoregressive orders p by regime
    if(min(x$ars$p)==max(x$ars$p)) cat("\nAutoregressive orders:",paste0(x$ars$p[1]," in each regime"))
    else cat("\nAutoregressive orders:",paste(x$ars$p,collapse=", "))
    # If exogenous variables are included, print the maximum lag orders q by regime
    if(max(x$ars$q)>0){
       if(min(x$ars$q)==max(x$ars$q)) cat("\nMaximum lags for ES  :",paste0(x$ars$q[1]," in each regime"))
       else cat("\nMaximum lags for ES  :",paste(x$ars$q,collapse=", "))
    }
    # If threshold effects are included, print the maximum lag orders d by regime
    if(max(x$ars$d)>0){
       if(min(x$ars$d)==max(x$ars$d)) cat("\nMaximum lags for TS  :",paste0(x$ars$d[1]," in each regime"))
       else cat("\nMaximum lags for TS  :",paste(x$ars$d,collapse=", "))
    }
}
#'
#' @method residuals mtar
#' @export
residuals.mtar <- function(object,...){
  # Number of MCMC simulations
  n.sim <- object$n.sim
  # Noise process distribution
  dist <- object$dist
  # Number of response variables
  k <- ncol(object$data[[1]]$y)
  # Vectors to store residuals and their original indices
  idsi <- resi <- vector()
  # Number of Monte Carlo draws for residual transformation
  TT <- 100000
  # Original observation indices
  ids <- 1:nrow(object$data[[1]]$X)
  # Regime classification
  a <- coef(object)
  if(object$regim > 1){
     lims <- (object$ps+1):(length(object$threshold.series))
     h <- as.integer(a$delay)
     thresholds <- a$thresholds
     Z <- object$threshold.series[lims-h]
     regs <- cut(Z,breaks=c(-Inf,sort(thresholds),Inf),labels=FALSE)
  }else regs <- matrix(1,nrow(object$data[[1]]$y),1)
  # Storage for skew-distribution residuals by component
  resus <- vector()
  # Loop over regimes
  for(i in 1:object$regim){
      location <- a[[i]]$location
      scale <- a[[i]]$scale
      # Select observations belonging to the current regime
      places <- regs==i
      y <- object$data[[i]]$y[places,]
      X <- object$data[[i]]$X[places,]
      if(object$ssvs){
         zeta <- rowMeans(object$chains[[i]]$zeta)
         X <- object$data[[i]]$X[places,zeta>0.5]
         location <- location[zeta>0.5,]
      }
      # Raw residuals
      resu <- y - X%*%location
      # Special treatment for skewed distributions
      if(dist %in% c("Skew-normal","Skew-Student-t")){
         resus0 <- vector()
         # Latent scale for Skew-Student-t distribution
         if(dist=="Skew-Student-t") u <- 1/rgamma(TT,shape=a$extra/2,rate=a$extra/2) else u <- rep(1,TT)
         # Posterior mean skewness parameters
         resu2 <- matrix(qnorm(0.5 + 0.5*runif(TT*k)),TT,k)*matrix(a$skewness,TT,k,byrow=TRUE) + crossprod(matrix(rnorm(TT*k),k,TT),chol(scale))
         # Simulated reference residuals
         resu2 <- resu2*matrix(sqrt(u),length(u),k)
         # Probability integral transform by component
         for(ii in 1:k){
             temp <- apply(matrix(resu[,ii],sum(places),1),1,function(x) mean(resu2[,ii]<=x))
             temp2 <- ifelse(temp>=0.5,1-temp,temp)
             temp <- qnorm(ifelse(.Machine$double.xmin>=temp2,.Machine$double.xmin,temp2))*ifelse(temp>0.5,-1,1)
             resus0 <- cbind(resus0,temp)
         }
         # Store componentwise residuals
         resus <- rbind(resus,resus0)
         # Radial residuals
         resu2 <- rowSums(resu2^2)
         resu <-  matrix(apply(matrix(rowSums(resu^2),sum(places),1),1,function(x) mean(resu2<=x)),sum(places),1)
      }else{
           # Standardization for symmetric distributions
           if(nrow(scale)==1) resu <- resu/as.numeric(sqrt(scale))
           else{
              ei <- eigen(scale,symmetric=TRUE)
              resu <- resu%*%(ei$vectors%*%diag(1/sqrt(as.vector(ei$values)))%*%t(ei$vectors))
           }
      }
      # Collect residuals and corresponding indices
      resi <- rbind(resi,resu)
      idsi <- c(idsi,ids[places])
  }
  # Radial residuals for symmetric distributions
  if(!(dist %in% c("Skew-normal","Skew-Student-t"))){
     resi <- matrix(resi,nrow(object$data[[1]]$X),k)
     # Simulated reference sample
     sim <- matrix(rnorm(TT*k),TT,k)
     # Latent scale depending on the distribution
     u <- switch(dist,"Gaussian"={rep(1,TT)},
                 "Student-t"={1/rgamma(TT,shape=a$extra/2,rate=a$extra/2)},
                 "Slash"={1/rbeta(TT,shape1=a$extra/2,shape2=1)},
                 "Contaminated normal"={1/(1 - (1-a$extra[2])*rbinom(TT,1,a$extra[1]))},
                 "Hyperbolic"={rgig(n=TT,lambda=1,chi=1,psi=a$extra^2)},
                 "Laplace"={rexp(TT,rate=1/8)})
     sim2 <- sort(rowSums(sim^2)*u)
     resi2 <- matrix(rowSums(resi^2),nrow(resi),1)
     resi2 <- apply(resi2,1,function(x) mean(sim2<=x))
  }else resi2 <- matrix(resi,nrow(object$data[[1]]$X),1)
  # Symmetrized probability integral transform
  resi3 <- ifelse(resi2>=0.5,1-resi2,resi2)
  resij <- qnorm(ifelse(.Machine$double.xmin>=resi3,.Machine$double.xmin,resi3))*ifelse(resi2>0.5,-1,1)
  # Restore original observation order
  idsx <- sort(idsi,index=TRUE)$ix
  resij <- matrix(resij[idsx],length(resij),1)
  # Componentwise residuals for symmetric distributions
  if(!(dist %in% c("Skew-normal","Skew-Student-t"))){
     sim <- sim*matrix(sqrt(u),TT,k)
     for(i in 1:k){
         temp <- apply(matrix(resi[,i],nrow(resi),1),1,function(x) mean(sim[,i]<=x))
         temp2 <- ifelse(temp>=0.5,1-temp,temp)
         resi[,i] <- qnorm(ifelse(.Machine$double.xmin>=temp2,.Machine$double.xmin,temp2))*ifelse(temp>0.5,-1,1)
     }
  }else resi <- matrix(resus,nrow(object$data[[1]]$y),k)
  # Reorder componentwise residuals
  resi <- matrix(resi[idsx,],nrow(resi),k)
  # Set row and column names
  row.names(resij) <- row.names(resi) <- row.names(object$data[[1]]$y)
  colnames(resij) <- " "
  colnames(resi) <- colnames(object$data[[1]]$y)
  # Output list with full and componentwise residuals
  out_ <- list(full=resij,by.component=resi)
  class(out_) <- "residuals_mtar"
  return(out_)
}
#'
#' @method plot residuals_mtar
#' @export
plot.residuals_mtar <- function(x,...){
  # Collect additional graphical arguments
  nano_ <- list(...)
  # Set default plotting character if not provided
  if(is.null(nano_$pch)) nano_$pch <- 20
  # Set default y-axis limits based on the range of full residuals
  if(is.null(nano_$ylim)) nano_$ylim <- range(x$full)
  # Set default color for points if not provided
  if(is.null(nano_$col)) nano_$col <- "black"
  # Open a new graphics device
  dev.new()
  # Set plotting layout to one row and two columns
  par(mfrow=c(1,2))
  # Normal Q–Q plot of the full (quantile-type) residuals
  qqnorm(x$full,pch=nano_$pch,col=nano_$col,main="",ylim=nano_$ylim)
  # Add a reference line with slope 1 and intercept 0
  abline(0,1,lty=3)
  # Extract full residuals for histogram
  x <- x$full
  # Histogram of quantile-type residuals with density scaling
  hist(x,freq=FALSE,xlab="Quantile-type residual",ylab="Density",main="",xlim=nano_$ylim)
  # Overlay the standard normal density curve
  curve(dnorm(x), col="blue", add=TRUE)
}

#' @method print residuals_mtar
#' @export
print.residuals_mtar <- function(x, digits=max(3, getOption("digits") - 2), ...){
  # Compute selected quantiles for the full residuals and each component residual
  x2 <- apply(cbind(x$full,x$by.component),2,quantile,probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
  # Set column names: one for full residuals and the rest for each component
  colnames(x2) <- c("Full",colnames(x$by.component))
  # Prefix row names to indicate quantile levels
  rownames(x2) <- paste0("Quantile ",rownames(x2))
  # Print the quantile table with the specified number of digits
  print(x2,digits=digits,na.print="")
}

#' @method model.matrix mtar
#' @export
#'
model.matrix.mtar <- function(object,...){
  out_ <- list()
  for(j in 1:object$regim) out_[[j]] <- object$data[[j]]$X
  return(out_)
}
