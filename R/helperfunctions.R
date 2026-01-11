#'
#' @title Auxiliary function for setting hyperparameter values
#' @description This function constructs and validates the list of hyperparameter values
#' used to define the prior distributions of the model parameters. Hyperparameters not
#' explicitly provided by the user are assigned their default values, which define
#' non-informative prior distributions.
#'
#' @param prior A list specifying user-defined values for the hyperparameters. Any
#' hyperparameters not included in this list are set to their default values.
#' @param regim A positive integer indicating the number of regimes in the model.
#' @param k     A positive integer indicating the dimension of the multivariate output series.
#' @param dist  A character string specifying the distribution chosen to model the noise process.
#' @param setar A positive integer indicating the component of the output series that acts as the
#' threshold variable in a SETAR specification. If \code{NULL} is specified then the
#' model is not a SETAR.
#' @param ssvs  A logical indicating whether the Stochastic Search Variable Selection (SSVS)
#' procedure should be applied.
#' @return A list containing the hyperparameter values defining the prior distributions of
#' all model parameters.
#' @export priors
#'
priors <- function(prior,regim,k,dist,setar,ssvs){
  # Set default prior mean for location parameter
  if(is.null(prior$mu0)) prior$mu0 <- 0
  if(!is.numeric(prior$mu0)) stop("The hyperparameter 'mu0' must be a numeric value",call.=FALSE)
  # Set default prior variance for location parameter
  if(is.null(prior$delta0)) prior$delta0 <- 1e9
  if(!is.numeric(prior$delta0) | prior$delta0<=0) stop("The hyperparameter 'delta0' must be a positive value",call.=FALSE)
  # Set default prior scale for scale parameter
  if(is.null(prior$omega0)) prior$omega0 <- 1e-9
  if(!is.numeric(prior$omega0) | prior$omega0<=0) stop("The hyperparameter 'omega0' must be a positive value",call.=FALSE)
  # Degrees of freedom parameter for inverse-Wishart prior
  if(is.null(prior$tau0)) prior$tau0 <- k
  if(prior$tau0 <= k-1) stop(paste("The value of the hyperparameter 'tau0' must be higher than",k-1),call.=FALSE)
  # Prior for skewness parameter if they are present
  if(dist %in% c("Skew-normal","Skew-Student-t")){
     if(is.null(prior$lambda0)) prior$lambda0 <- 1e9
  }
  # Priors related to regime switching and threshold parameters
  if(regim > 1){
     if(is.null(setar)){
        # Minimum delay for TAR models
        if(is.null(prior$hmin)) prior$hmin <- 0
     }
     else{
        # Minimum delay for SETAR models must be positive
        if(is.null(prior$hmin)) prior$hmin <- 1
        if(prior$hmin==0) stop("For SETAR models the delay parameter 'h' must be positive!",call.=FALSE)
     }
     # Maximum delay parameter
     if(is.null(prior$hmax)) prior$hmax <- 3
     # Priors for threshold distribution
     if(is.null(prior$alpha0)) prior$alpha0 <- 0
     if(is.null(prior$alpha1)) prior$alpha1 <- 1
     # Validate delay parameters
     if(floor(prior$hmin)!=prior$hmin | floor(prior$hmax)!=prior$hmax | prior$hmin < 0 | prior$hmax < 0)
        stop("The hyperparameters 'hmin' and 'hmax' must be positive integer values",call.=FALSE)
     if(prior$hmin > prior$hmax)
        stop("The hyperparameter 'hmin' must be lower than or equal to 'hmax'",call.=FALSE)
     # Validate threshold bounds
     if(prior$alpha0 < 0 | prior$alpha1 < 0 | prior$alpha0 > 1 | prior$alpha1 > 1)
        stop("The hyperparameters 'alpha0' and 'alpha1' must be within the interval (0,1)",call.=FALSE)
     if(prior$alpha0 > prior$alpha1)
        stop("The hyperparameter 'alpha0' must be lower than 'alpha1'",call.=FALSE)
  }
  # Prior inclusion probability for SSVS
  if(ssvs){
     if(is.null(prior$rho0)) prior$rho0 <- 0.5
     if(prior$rho0 < 0 | prior$rho0 > 1) stop("The hyperparameter 'rho0' must be within the interval (0,1)",call.=FALSE)
  }
  # Distribution-specific hyperparameters
  switch(dist,
         "Hyperbolic"={if(is.null(prior$gamma0)) prior$gamma0 <- 0.1
                      if(is.null(prior$eta0)) prior$eta0 <- 4},
         "Skew-Student-t"=,
         "Student-t"={if(is.null(prior$gamma0)) prior$gamma0 <- 2
                      if(is.null(prior$eta0)) prior$eta0 <- 100},
         "Slash"={if(is.null(prior$gamma0)) prior$gamma0 <- 1e-9
                  if(is.null(prior$eta0)) prior$eta0 <- 1e-9},
         "Contaminated normal"={if(is.null(prior$gamma01)) prior$gamma01 <- 1e-9
                                if(is.null(prior$eta01)) prior$eta01 <- 1e-9
                                if(is.null(prior$gamma02)) prior$gamma02 <- 1e-9
                                if(is.null(prior$eta02)) prior$eta02 <- 1e-9})
  # Validate hyperparameters for extra parameters
  if(dist %in% c("Student-t","Skew-Student-t","Hyperbolic","Slash")){
     if(!is.numeric(prior$gamma0) | prior$gamma0 <= 0) stop("The hyperparameter 'gamma0' must be a positive value",call.=FALSE)
     if(!is.numeric(prior$eta0) | prior$eta0 <= 0) stop("The hyperparameter 'eta0' must be a positive value",call.=FALSE)
  }
  # Validate hyperparameters for extra parameters of contaminated normal distribution
  if(dist=="Contaminated normal"){
     if(!is.numeric(prior$gamma01) | prior$gamma01 <= 0) stop("The hyperparameter 'gamma01' must be a positive value",call.=FALSE)
     if(!is.numeric(prior$eta01) | prior$eta01 <= 0) stop("The hyperparameter 'eta01' must be a positive value",call.=FALSE)
     if(!is.numeric(prior$gamma02) | prior$gamma02 <= 0) stop("The hyperparameter 'gamma02' must be a positive value",call.=FALSE)
     if(!is.numeric(prior$eta02) | prior$eta02 <= 0) stop("The hyperparameter 'eta02' must be a positive value",call.=FALSE)
  }
  # Return the validated and completed prior list
  return(prior)
}
#'
#' @title Auxiliary function to specify the number of regimes and lag orders
#' @description This auxiliary function defines the regime structure of a multivariate TAR model by
#' specifying the number of regimes and the corresponding lag orders for the endogenous,
#' exogenous, and threshold series in each regime.
#' @param nregim A positive integer indicating the total number of regimes.
#' @param p A list of positive integers specifying the autoregressive order of the output series
#' within each regime.
#' @param q A list of non-negative integers specifying the maximum lag of the exogenous series
#' within each regime.
#' @param d A list of non-negative integers specifying the maximum lag of the threshold series
#' within each regime.
#' @return A list containing the number of regimes and the regime-specific lag-order specifications.
#' @export ars
#'
ars <- function(nregim=1,p=1,q=0,d=0){
  # Check that the number of regimes is a positive integer
  if(nregim!=floor(nregim) | nregim<=0)
     stop("The argument 'nregim' must be a positive integer value!",call.=FALSE)
  # Truncate p, q, and d to at most nregim elements
  p <- p[1:min(nregim,length(p))]
  q <- q[1:min(nregim,length(q))]
  d <- d[1:min(nregim,length(d))]
  # Validate autoregressive order p (must be positive integers)
  if(any(p!=floor(p)) | any(p<=0))
     stop("The argument 'p' must be a vector of positive integer values",call.=FALSE)
  # Validate moving average (q) and delay (d) orders (must be non-negative integers)
  if(any(q!=floor(q)) | any(q<0) | any(d!=floor(d)) | any(d<0))
     stop("The arguments 'q' and 'd' must be vectors of non-negative integer values",call.=FALSE)
  # Replicate single values of p, q, and d across regimes if needed
  if(length(p)==1 & nregim > 1) p <- rep(p[1],nregim)
  if(length(q)==1 & nregim > 1) q <- rep(q[1],nregim)
  if(length(d)==1 & nregim > 1) d <- rep(d[1],nregim)
  # Fill missing regime values with defaults
  if(length(p) < nregim) p[(length(p)+1):nregim] <- rep(1,nregim-length(p))
  if(length(q) < nregim) q[(length(q)+1):nregim] <- rep(0,nregim-length(q))
  if(length(d) < nregim) q[(length(d)+1):nregim] <- rep(0,nregim-length(d))
  if(nregim==1) d <- rep(0,nregim)
  # Return the AR structure as a list
  return(list(nregim=nregim,p=p,q=q,d=d))
}
#'


