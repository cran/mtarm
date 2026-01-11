#' @title Bayesian Estimation of Multivariate TAR Models
#'
#' @description
#' This function is a wrapper that applies \code{mtar()} over a grid of model specifications
#' defined by all combinations of the noise distribution (\code{dist}), the number of regimes
#' (from \code{nregim.min} to \code{nregim.max}), the autoregressive order within each regime
#' (from \code{p.min} to \code{p.max}), the maximum lag of the exogenous series within each regime
#' (from \code{q.min} to \code{q.max}), and the maximum lag of the threshold series within each
#' regime (from \code{d.min} to \code{d.max}).
#' In all calls to \code{mtar()}, the same set of time points is used for model fitting. This is
#' achieved by appropriately adjusting the \code{subset} argument of \code{mtar()} for each model
#' specification, thereby ensuring comparability across models.
#'
#' @param formula A three-part expression of class \code{Formula} describing the TAR model to be fitted.
#' The first part specifies the variables in the multivariate output series, the second part
#' defines the threshold series, and the third part specifies the variables in the multivariate
#' exogenous series.
#'
#' @param nregim.min An optional integer specifying the minimum number of regimes. By default,
#' \code{nregim.min} is set to \code{1}.
#'
#' @param nregim.max An integer specifying the maximum number of regimes.
#'
#' @param p.min	An optional integer specifying the minimum autoregressive order within each regime.
#' By default, \code{p.min} is set to \code{1}.
#'
#' @param p.max	An integer specifying the maximum autoregressive order within each regime.
#'
#' @param q.min	An optional integer specifying the minimum value of the maximum lag of the exogenous
#' series within each regime. By default, \code{q.min} is set to \code{0}.
#'
#' @param q.max	An optional integer specifying the maximum value of the maximum lag of the exogenous
#' series within each regime. By default, \code{q.max} is set to \code{0}.
#'
#' @param d.min	An optional integer specifying the minimum value of the maximum lag of the threshold
#' series within each regime. By default, \code{d.min} is set to \code{0}.
#'
#' @param d.max	An optional integer specifying the maximum value of the maximum lag of the threshold
#' series within each regime. By default, \code{d.max} is set to \code{0}.
#'
#' @param Intercept An optional logical indicating whether an intercept should be included within each regime.
#'
#' @param trend	An optional character string specifying the degree of deterministic time trend to be
#' included in each regime. Available options are \code{"linear"}, \code{"quadratic"}, and
#' \code{"none"}. By default, \code{trend} is set to \code{"none"}.
#'
#' @param nseason An optional integer, greater than or equal to 2, specifying the number of seasonal periods.
#' When provided, \code{nseason - 1} seasonal dummy variables are added to the regressors within each regime.
#' By default, \code{nseason} is set to \code{NULL}, thereby indicating that the TAR model has no seasonal effects.
#'
#' @param data A data frame containing the variables in the model. If not found in \code{data}, the
#' variables are taken from \code{environment(formula)}, typically the environment from
#' which \code{mtar_grid()} is called.
#'
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#'
#' @param dist A character vector specifying the multivariate distributions used to model the noise
#' process. Available options are \code{"Gaussian"}, \code{"Student-t"}, \code{"Slash"},
#' \code{"Hyperbolic"}, \code{"Laplace"}, \code{"Contaminated normal"},
#' \code{"Skew-normal"}, and \code{"Skew-Student-t"}. By default, \code{dist} is set to
#' \code{"Gaussian"}.
#'
#' @param n.sim	An optional positive integer specifying the number of simulation iterations after the
#' burn-in period. By default, \code{n.sim} is set to \code{500}.
#'
#' @param n.burnin An optional positive integer specifying the number of burn-in iterations. By default,
#' \code{n.burnin} is set to \code{100}.
#'
#' @param n.thin An optional positive integer specifying the thinning interval. By default,
#' \code{n.thin} is set to \code{1}.
#'
#' @param row.names An optional variable in \code{data} labelling the time points corresponding to each row of the data set.
#'
#' @param prior	An optional list specifying the hyperparameter values that define the prior
#' distribution. This list can be validated using the \code{priors()} function. By default,
#' \code{prior} is set to an empty list, thereby indicating that the hyperparameter values
#' should be set so that a non-informative prior distribution is obtained.
#'
#' @param ssvs An optional logical indicating whether the Stochastic Search Variable Selection (SSVS)
#' procedure should be applied to identify relevant lags of the output, exogenous, and threshold
#' series. By default, \code{ssvs} is set to \code{FALSE}.
#'
#' @param setar An optional positive integer indicating the component of the output series used as the
#' threshold variable. By default, \code{setar} is set to \code{NULL}, indicating that the
#' fitted model is not a SETAR model.
#'
#' @param plan_strategy An optional character string specifying the execution strategy for parallel
#' computation. Available options are \code{"sequential"} and \code{"multisession"}. By default,
#' \code{plan_strategy} is set to \code{"sequential"}.
#'
#' @param progress An optional logical indicating whether a progress bar should be displayed during
#' execution. By default, \code{progress} is set to \code{TRUE}.
#'
#' @return A list whose elements are objects of class \code{mtar}, each corresponding to a distinct
#' model specification considered in the grid.
#'
#' @seealso mtar
#' @export mtar_grid
mtar_grid <- function(formula, data, subset, Intercept=TRUE, trend=c("none","linear","quadratic"), nseason=NULL,
                      nregim.min=1, nregim.max=NULL, p.min=1, p.max=NULL, q.min=0, q.max=0, d.min=0,
                      d.max=0, row.names, dist="Gaussian", prior=list(), n.sim=500, n.burnin=100,
                      n.thin=1, ssvs=FALSE, setar=NULL, plan_strategy=c("sequential","multisession"), progress=TRUE){
  # Check logical arguments for validity and correct length
  if(!is.logical(Intercept) | length(Intercept)!= 1) stop("'Intercept' must be a single logical value",call.=FALSE)
  if(!is.logical(ssvs) | length(ssvs)!=1) stop("'ssvs' must be a single logical value",call.=FALSE)
  if(!is.logical(progress) | length(progress)!=1) stop("'progress' must be a single logical value",call.=FALSE)
  # Select parallelization strategy
  plan_strategy <- match.arg(plan_strategy)
  # Remove duplicated distribution names
  dist <- unique(dist)
  # Check supported error distributions
  if(!all(dist %in% c("Gaussian","Student-t","Hyperbolic","Laplace","Slash","Contaminated normal","Skew-Student-t","Skew-normal")))
     stop("Only 'Gaussian', 'Student-t', 'Hyperbolic', 'Laplace', 'Slash', 'Contaminated normal',
          'Skew-Student-t' and 'Skew-normal' distributions are supported",call.=FALSE)
  # Reset gamma0 and eta0 if incompatible distributions are jointly requested
  if((sum(dist %in% c("Student-t","Hyperbolic","Slash"))>1 | sum(dist %in% c("Skew-Student-t","Hyperbolic","Slash"))>1)
     & (!is.null(prior$gamma0) | !is.null(prior$eta0))){
    message("\nThe hyperparameters 'gamma0' and 'eta0' are set to their by-default values\n")
    prior$gamma0 <- prior$eta0 <- NULL
  }
  # Validate minimum and maximum number of regimes
  if(nregim.min!=floor(nregim.min) | nregim.min<=0)
     stop("The argument 'nregim.min' must be a positive integer value!",call.=FALSE)
  if(is.null(nregim.max)) stop("The argument 'nregim.max' is required!",call.=FALSE)
  if(nregim.max!=floor(nregim.max) | nregim.max<=0 | nregim.max<nregim.min)
     stop("The argument 'nregim.max' must be a positive integer greater than or equal to nregim.min!",call.=FALSE)
  # Validate minimum and maximum autoregressive orders
  if(p.min!=floor(p.min) | p.min<=0)
     stop("The argument 'p.min' must be a positive integer value!",call.=FALSE)
  if(is.null(p.max)) stop("The argument 'p.max' is required!",call.=FALSE)
  if(p.max!=floor(p.max) | p.max<=0 | p.max<p.min)
     stop("The argument 'p.max' must be a positive integer greater than or equal to p.min!",call.=FALSE)
  # Validate minimum and maximum q
  if(q.min!=floor(q.min) | q.min<0)
     stop("The argument 'q.min' must be a non-negative integer value!",call.=FALSE)
  if(q.max!=floor(q.max) | q.max<0 | q.max<q.min)
     stop("The argument 'q.max' must be a non-negative integer greater than or equal to q.min!",call.=FALSE)
  # Validate minimum and maximum d
  if(d.min!=floor(d.min) | d.min<0)
     stop("The argument 'd.min' must be a non-negative integer value!",call.=FALSE)
  if(d.max!=floor(d.max) | d.max<0 | d.max<d.min)
     stop("The argument 'd.max' must be a non-negative integer greater than or edual to d.min!",call.=FALSE)
  # Build the grid of models to be estimated
  grid <- expand.grid(dist=dist,l=nregim.min:nregim.max,p=p.min:p.max,q=q.min:q.max,d=d.min:d.max,stringsAsFactors=FALSE)
  grid <- grid[order(grid$dist, grid$l, grid$p, grid$q, grid$d), ]
  # Prepare a call to mtar() to extract data structure and dimensions
  mycall <- match.call()
  mycall[[1]] <- as.name("mtar")
  mycall$ars <- ars(nregim=1,p=1,q=q.max,d=d.max)
  mycall$stop <- TRUE
  mycall$dist <- "Gaussian"
  temp <- eval(mycall)
  # Extract data and dimension of the response
  mydata <- temp$data
  k <- temp$k
  # Clean temporary arguments from the call
  mycall$subset <- mycall$stop <- NULL
  if(!is.null(temp$mynames)) mycall$row.names <- as.name("mynames")
  mycall$data <- data.frame(mydata,mynames=temp$mynames)
  mycall$progress <- FALSE
  # Determine maximum lag needed to align estimation samples
  hmax <- ifelse(is.null(prior$hmax),3,prior$hmax)
  ps <- max(max(grid$p,grid$q,grid$d),ifelse(max(grid$l)>1,hmax,0))
  # Compute required trimming of initial observations for each grid point
  grid$subset <- ps-apply(matrix(cbind(grid$p,grid$q,grid$d,ifelse(grid$l>1,hmax,0)),nrow(grid),4),1,max)
  grid$subset <- ifelse(grid$subset>0,grid$subset,NA)
  mycall$ssvs <- ifelse(missingArg(ssvs),FALSE,ssvs)
  # Configure progress bar handling
  if(progress) handlers("txtprogressbar") else handlers("void")
  # Set parallel execution plan
  future::plan(plan_strategy)
  # Run model estimation over the grid with progress reporting
  with_progress({
    pbg <- progressor(steps=nrow(grid))
    out <- future_lapply(1:nrow(grid),
    function(i){
      # Extract current grid configuration
      nowpar <- grid[i,]
      # Set prior distribution for current model
      mycall$prior <- priors(prior,regim=nowpar$l,k=k,dist=nowpar$dist,setar=mycall$setar,ssvs=mycall$ssvs)
      mycall$dist <- nowpar$dist
      # Specify AR structure
      mycall$ars <- ars(nregim=nowpar$l,p=nowpar$p,q=nowpar$q,d=nowpar$d)
      # Apply subset trimming if needed
      if(is.na(nowpar$subset)) mycall$subset <- NULL
      else mycall$subset <- -c(1:nowpar$subset)
      # Estimate the model
      fitm <- eval(mycall)
      # Update progress bar if specified
      if(progress) pbg()
      return(fitm)
    },future.seed=TRUE)
  })
  # Assign names and class to output list
  if(q.max>0){
     if(d.max>0) names(out) <- paste0(grid$dist,".",grid$l,".",grid$p,".",grid$q,".",grid$d)
     else names(out) <- paste0(grid$dist,".",grid$l,".",grid$p,".",grid$q)
  }else{
     if(d.max>0) names(out) <- paste0(grid$dist,".",grid$l,".",grid$p,".",grid$d)
     else names(out) <- paste0(grid$dist,".",grid$l,".",grid$p)
  }
  class(out) <- "listmtar"
  # Return fitted models
  return(out)
}
#'
#' @title Geweke's convergence diagnostic for \code{mtar} objects
#' @description This function computes Geweke's convergence diagnostic for Markov chain Monte Carlo
#' (MCMC) output obtained from Bayesian estimation of multivariate TAR models. It is a
#' wrapper around \code{geweke.diag()} that applies the diagnostic to the posterior chains
#' returned by a call to \code{mtar()}.
#'
#' @param x An object of class \code{mtar} returned by the function \code{mtar()}.
#' @param frac1 A numeric value in \eqn{(0,1)} specifying the fraction of the initial part
#'              of each chain to be used in the diagnostic.
#' @param frac2 A numeric value in \eqn{(0,1)} specifying the fraction of the final part
#'              of each chain to be used in the diagnostic.
#'
#' @return
#' A list containing the Geweke z-scores for the parameters of the \code{mtar} model.
#' @seealso \code{\link[coda]{geweke.diag}}
#' @export
#' geweke.diagTAR
geweke.diagTAR <- function(x,frac1=0.1,frac2=0.5){
  # Check that the input object is of class 'mtar'
  if(!inherits(x,"mtar")) stop("Only objects of class 'mtar' are supported!",call.=FALSE)
  # Convert the mtar object into an mcmc object
  x2 <- as.mcmc(x)
  # Remove the delay parameter from the MCMC output
  x2$delay <- NULL
  # If the model has more than one regime, compute Geweke diagnostics for thresholds
  if(x$ars$nregim > 1){
     x2$thresholds <- geweke.diag(x2$thresholds,frac1=frac1,frac2=frac2)$z
     temp <- matrix(x2$thresholds,length(x2$thresholds),1)
     rownames(temp) <- names(x2$thresholds); colnames(temp) <- ""
     x2$thresholds <- temp
  }
  # Initialize a data frame to collect Geweke statistics for location parameters
  temp <- data.frame(ns=NA)
  # Initialize a matrix to collect Geweke statistics for scale parameters
  temp2 <- matrix(NA,ncol(x2$scale[[1]]),1)
  rownames(temp2) <- colnames(x2$scale[[1]])
  # Define name mappings to ensure proper ordering of parameters
  cc <- c(colnames(x$data[[1]]$y),":Time",":Season")
  cc2 <- c(paste0(1:ncol(x$data[[1]]$y),colnames(x$data[[1]]$y)),":.1Time",":.2Season")
  # Loop over regimes to compute Geweke diagnostics for location and scale parameters
  for(j in 1:x$ars$nregim){
      x2$location[[j]] <- geweke.diag(x2$location[[j]],frac1=frac1,frac2=frac2)$z
      x2$location[[j]] <- data.frame(ns=names(x2$location[[j]]),ps=x2$location[[j]])
      temp <- merge(temp,x2$location[[j]],by.x="ns",by.y="ns",all.x=TRUE,all.y=TRUE)
      x2$scale[[j]] <- geweke.diag(x2$scale[[j]],frac1=frac1,frac2=frac2)$z
      temp2 <- cbind(temp2,matrix(x2$scale[[j]],length(x2$scale[[j]]),1))
  }
  # Remove rows with missing parameter names
  temp <- temp[!is.na(temp[,1]),]
  # Temporarily modify names to allow correct sorting
  for(jj in 1:(length(cc))) temp[,1] <- gsub(cc[jj],cc2[jj],temp[,1])
  temp <- temp[sort(temp[,1],index=TRUE)$ix,]
  for(jj in 1:(length(cc))) temp[,1] <- gsub(cc2[jj],cc[jj],temp[,1])
  # Store the Geweke statistics for location parameters as a matrix
  x2$location <- as.matrix(temp[,-1])
  rownames(x2$location) <- temp[,1]
  # Store the Geweke statistics for scale parameters as a matrix
  x2$scale <- matrix(temp2[,-1],nrow(temp2),x$ars$nregim)
  rownames(x2$scale) <- rownames(temp2)
  # Assign regime-specific column names
  colnames(x2$location) <- colnames(x2$scale) <- paste0("Regime ",1:x$ars$nregim)
  # Compute Geweke diagnostics for extra parameters if they are present
  if(!(x$dist %in% c("Gaussian","Laplace","Skew-normal"))){
     x2$extra <- geweke.diag(x2$extra,frac1=frac1,frac2=frac2)$z
     temp <- matrix(x2$extra,length(x2$extra),1)
     rownames(temp) <- names(x2$extra); colnames(temp) <- ""
     x2$extra <- temp
  }
  # Compute Geweke diagnostics for skewness parameters if they are present
  if(x$dist %in% c("Skew-normal","Skew-Student-t")){
     x2$skewness <- geweke.diag(x2$skewness,frac1=frac1,frac2=frac2)$z
     temp <- matrix(x2$skewness,length(x2$skewness),1)
     rownames(temp) <- names(x2$skewness); colnames(temp) <- ""
     x2$skewness <- temp
  }
  # Store the fractions used in the Geweke diagnostic
  x2$frac <- c(frac1,frac2)
  # Assign the output class
  class(x2) <- "gdmtar"
  return(x2)
}
#'
#'
#' @method print gdmtar
#' @export
print.gdmtar <- function(x, digits=max(3, getOption("digits") - 2), ...){
  # Print the fractions used for the first and second Geweke windows
  cat("\nFraction in 1st window = ",x$frac[1])
  cat("\nFraction in 2nd window = ",x$frac[2],"\n\n")
  # Print Geweke diagnostics for thresholds if they are present
  if(!is.null(x$thresholds)){
     cat("Thresholds:\n")
     print(x$thresholds,digits=digits,na.print="")
  }
  # Print Geweke diagnostics for location parameters
  cat("\n\nAutoregressive coefficients:\n")
  print(x$location,digits=digits,na.print="")
  # Print Geweke diagnostics for the scale parameters
  cat("\n\nScale parameter:\n")
  print(x$scale,digits=digits,na.print="")
  # Print Geweke diagnostics for skewness parameters if they are present
  if(!is.null(x$skewness)){
     cat("\n\nSkewness parameter:\n")
     print(x$skewness,digits=digits,na.print="")
  }
  # Print Geweke diagnostics for extra parameters if they are present
  if(!is.null(x$extra)){
     cat("\n\nExtra parameter:\n")
     print(x$extra,digits=digits,na.print="")
  }
}
#'
#' @title Geweke-Brooks plot for objects of class \code{mtar}
#' @description	This function is a wrapper around \code{geweke.plot()} that applies the
#' Geweke-Brooks convergence diagnostic to the MCMC chains obtained from a
#' fitted \code{mtar} model.
#' @param x An object of class \code{mtar} returned by a call to \code{mtar()}.
#' @param frac1	fraction to use from beginning of chain
#' @param frac2	fraction to use from end of chain
#' @param nbins Number of segments
#' @param pvalue p-value used to plot confidence limits for the null hypothesis
#' @param auto.layout	If \code{TRUE} then, set up own layout for plots, otherwise use existing one
#' @param ask If \code{TRUE} then prompt user before displaying each page of plots. Default is
#' \code{dev.interactive()}.
#' @param ... Additional graphical parameters passed to the plotting routines.
#' @seealso \code{\link[coda]{geweke.plot}}
#' @export geweke.plotTAR
geweke.plotTAR <- function(x, frac1=0.1, frac2=0.5, nbins=20, pvalue=0.05, auto.layout=TRUE, ask, ...){
  # Check that the input object is of class 'mtar'
  if(!inherits(x,"mtar")) stop("Only objects of class 'mtar' are supported!",call.=FALSE)
  # Convert the mtar object to an mcmc object for diagnostic analysis
  x2 <- as.mcmc(x)
  # Plot Geweke diagnostics for threshold parameters if they are present
  if(!is.null(x2$thresholds)){
     dev.new()
     cat("\nThresholds\n")
     geweke.plot(x2$thresholds,frac1=frac1,frac2=frac2,nbins=nbins,pvalue=pvalue,
                 auto.layout=auto.layout,ask=ask,...)
  }
  # Loop over regimes to plot diagnostics for each regime separately
  for(j in 1:x$regim){
      # Plot Geweke diagnostics for location parameters
      dev.new()
      cat(paste("\nAutoregressive coefficients Regime",j,"\n"))
      geweke.plot(x2$location[[j]],frac1=frac1,frac2=frac2,nbins=nbins,pvalue=pvalue,
                  auto.layout=auto.layout,ask=ask,...)
      # Plot Geweke diagnostics for scale parameters
      dev.new()
      cat(paste("\nScale parameter Regime",j,"\n"))
      geweke.plot(x2$scale[[j]],frac1=frac1,frac2=frac2,nbins=nbins,pvalue=pvalue,
                  auto.layout=auto.layout,ask=ask,...)
  }
  # Plot Geweke diagnostics for skewness parameters if they are present
  if(!is.null(x2$skewness)){
     dev.new()
     cat("\nSkewness parameter\n")
     geweke.plot(x2$skewness,frac1=frac1,frac2=frac2,nbins=nbins,pvalue=pvalue,
                 auto.layout=auto.layout,ask=ask,...)
  }
  # Plot Geweke diagnostics for extra parameters if they are present
  if(!is.null(x2$extra)){
     dev.new()
     cat("\nExtra parameter\n")
     geweke.plot(x2$extra,frac1=frac1,frac2=frac2,nbins=nbins,pvalue=pvalue,
                 auto.layout=auto.layout,ask=ask,...)
  }
}
