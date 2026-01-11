#' @importFrom grDevices dev.new dev.interactive
#' @importFrom graphics abline hist lines par identify curve legend
#' @importFrom methods missingArg
#' @importFrom stats pnorm pt density median model.extract model.matrix pgamma qgamma qnorm quantile rbeta rbinom rexp rgamma rnorm runif sd terms na.omit qqnorm dnorm predict resid fitted coef
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Formula Formula model.part
#' @importFrom mvtnorm pmvnorm pmvt
#' @importFrom GIGrvg rgig
#' @importFrom coda mcmc as.mcmc geweke.diag geweke.plot HPDinterval
#' @importFrom future.apply future_lapply
#' @importFrom progressr handlers progressor with_progress
#' @importFrom future plan
#'
#' @title Temperature, precipitation, and two river flows in Iceland
#'
#' @description This data set contains two daily river-flow time series, measured in cubic
#' meters per second, for rivers in Iceland from January 1, 1972, to December 12, 1974.
#' Additionally, daily measurements of precipitation (in millimeters) and temperature
#' (in degrees Celsius) were recorded at the Hveravellir meteorological station. The
#' precipitation values correspond to the accumulated precipitation up to 9:00 A.M.
#' relative to the same time on the previous day.
#'
#' @docType data
#'
#' @usage data(iceland.rf)
#'
#' @format A data frame with 1,096 rows and 5 variables:
#' \describe{
#'   \item{Vatnsdalsa}{Numeric vector representing the daily flow of the Vatnsdalsá river.}
#'   \item{Jokulsa}{Numeric vector representing the daily flow of the Jökulsá Eystri river.}
#'   \item{Precipitation}{Numeric vector of daily precipitation amounts (millimeters).}
#'   \item{Temperature}{Numeric vector of daily temperature measurements (degrees Celsius).}
#'   \item{Date}{Vector indicating the calendar date of each observation.}
#' }
#' @keywords datasets
#' @references Tong, Howell (1990) Non‑linear Time Series: A Dynamical System Approach.
#'             Oxford University Press. Oxford, UK.
#' @references Ruey S., Tsay (1998) Testing and Modeling Multivariate Threshold Models.
#'             Journal of the American Statistical Association, 93, 1188-1202.
#' @examples
#' data(iceland.rf)
#' dev.new()
#' plot(ts(as.matrix(iceland.rf[,-5])), main="Iceland")
#'
"iceland.rf"
#'
#' @title Returns of the closing prices of three financial indexes
#'
#' @description This dataset contains daily returns computed from the closing prices of
#' the COLCAP, BOVESPA, and S&P 500 stock market indexes over the period from
#' February 10, 2010, to March 31, 2016, comprising 1,505 observations.
#' The COLCAP index reflects the price dynamics of the 20 most liquid stocks traded on
#' the Colombian stock market. The BOVESPA index represents the Brazilian stock market,
#' the largest and most important exchange in Latin America and among the largest
#' worldwide. The Standard & Poor's 500 (S&P 500) index tracks the performance of
#' 500 large-cap companies listed in the United States.
#'
#' @docType data
#'
#' @usage data(returns)
#'
#' @format A data frame with 1,505 rows and 4 variables:
#' \describe{
#'   \item{Date}{A vector indicating the date of each observation.}
#'   \item{COLCAP}{A numeric vector containing the returns of the COLCAP index.}
#'   \item{SP500}{A numeric vector containing the returns of the S&P 500 index.}
#'   \item{BOVESPA}{A numeric vector containing the returns of the BOVESPA index.}
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
#' @description This dataset contains daily observations of rainfall (in millimeters)
#' and river flows (in \eqn{m^3}/s) for two rivers in southern Colombia. Rainfall was
#' recorded at a meteorological station located at an altitude of 2400 meters above sea
#' level. River flows were measured at two hydrological stations: El Trébol station,
#' which monitors the Bedón River at an altitude of 1720 meters, and Villalosada station,
#' which monitors the La Plata River at an altitude of 1300 meters.
#' The stations are located near the equator, a geographic feature that helps reduce
#' seasonal distortions and facilitates the analysis of the dynamic relationship between
#' rainfall and river flows. The sample period spans from January 1, 2006, to April 14, 2009.
#'
#' @docType data
#'
#' @usage data(riverflows)
#'
#' @format A data frame with 1200 rows and 4 variables:
#' \describe{
#'   \item{Date}{A vector indicating the date of each observation.}
#'   \item{Bedon}{A numeric vector giving the daily flow of the Bedón River.}
#'   \item{LaPlata}{A numeric vector giving the daily flow of the La Plata River.}
#'   \item{Rainfall}{A numeric vector indicating daily rainfall amounts.}
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
#' @title U.S. Stock Returns
#'
#' @description The dataset comprises observations of both continuously compounded and
#' simple returns derived from the S&P 500 index, along with the first difference of
#' the Chicago Board Options Exchange Market Volatility Index (VIX). The sample period
#' spans from January 5, 1990, to March 30, 2012.
#'
#' @docType data
#'
#' @usage data(US.returns)
#'
#' @format A data frame with 5606 rows and 4 variables:
#' \describe{
#'   \item{Date}{A vector indicating the date of each observation.}
#'   \item{CCR}{A numeric vector giving the continuously compounded returns.}
#'   \item{SR}{A numeric vector giving the simple returns.}
#'   \item{dVIX}{A numeric vector giving (the first difference of) the Chicago
#'               Board Options Exchange Market Volatility Index (VIX).}
#' }
#' @keywords datasets
#' @references Massacci, D. (2014) A two-regime threshold model with conditional skewed
#'             student-t distributions for stock returns. Economic Modelling, 43, 9-20.
#' @examples
#' data(US.returns)
#' dev.new()
#' plot(ts(as.matrix(US.returns[,-1])), main="Returns and (the first difference of) VIX")
#'
"US.returns"
#'

