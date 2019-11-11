#' @title Fixed-b sup Wald test for a single change in the mean of a long-memory time series.
#' @description This function performs a sup-Wald test on a change-in-mean, which is standardized by a non-parametric kernel-based long-run variance estimator.
#' Therefore, the test is robust under long-memory.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15} using the Bartlett kernel.
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param d integer that specifies the long-memory parameter.
#' @param bandw integer that determines the bandwidth parameter for the long-run variance estimator. It can take values in the range \code{bandw=[0.05,0.1,0.2]}. Default is
#' \code{bandw=0.1}, which is suggested by Iacone, Leybourne and Taylor (2014).
#' @param tau integer that defines the search area, which is \code{[tau,1-tau]}. Default is \code{tau=0.15} as suggested by Andrews (1993).
#' @return Returns a numeric vector containing the test statistic and the corresponding critical values of the test.
#' @seealso \code{\link{CUSUMfixed}}, \code{\link{snsupwald}}
#' @author Kai Wenger
#' @examples
#' # set model parameters
#' T        <- 500
#' d        <- 0.2
#' 
#' set.seed(410)
#' 
#' # generate a fractionally integrated (long-memory) time series
#' tseries  <- fracdiff::fracdiff.sim(n=T, d=d)$series
#'
#' # generate a fractionally integrated (long-memory) time series
#' #  with a change in mean in the middle of the series
#' changep  <- c(rep(0,T/2), rep(1,T/2))
#' tseries2 <- tseries+changep
#' 
#' # estimate the long-memory parameter of both series via local 
#' # Whittle approach. The bandwidth to estimate d is chosen 
#' # as T^0.65, which is usual in literature
#' d_est    <- LongMemoryTS::local.W(tseries, m=floor(1+T^0.65))$d
#' d_est2   <- LongMemoryTS::local.W(tseries2, m=floor(1+T^0.65))$d
#'
#' # perform the test on both time series
#' fixbsupw(tseries, d=d_est)
#' fixbsupw(tseries2, d=d_est2)
#' # For the series with no change in mean the test does not reject 
#' # the null hypothesis of a constant mean across time at any 
#' # reasonable significance level.
#' # For the series with a change in mean the test rejects the 
#' # null hypothesis at a 1% significance level.
#' @references
#' Iacone, F. and Leybourne, S. J. and Taylor, R. A. M. (2014): A fixed-b Test for a Break in Level at an unknown Time under Fractional Integration. Journal of Time Series Analysis, 35, pp. 40-54.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

fixbsupw <- function(x,d,bandw=0.1,tau=0.15)
{
  if(!(bandw %in% c(0.05,0.1,0.2))) stop("Use one of the following bandwidths: 0.05,0.1,0.2")
  if (any(is.na(x)))
    stop("x contains missing values")
  if(tau<=0 | tau>=1)
    stop("It must hold that 0<tau<1")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x) == 
      FALSE) 
    stop("x must be a univariate numeric vector")
  if(tau!=0.15)
    warning("Critical values are just implemented for tau=0.15")
  
  T     <- length(x)
  m     <- bandw*T
  cc    <- t(c(0,1))
  out   <- c()

  for(k in ((T*tau):(T*(1-tau))))
  {
    X     <- cbind(rep(1,T),c(rep(0,round(k)),rep(1,round(T-k))))
    reg   <- stats::lm(x~X-1)
    u_hat <- unname(reg$residuals)
    beta  <- cc%*%reg$coefficients
    sigm  <- fb_longrun(u_hat,m)
    out   <- c(out,(beta^2)/(sigm*(cc%*%solve(t(X)%*%X))%*%t(cc)))
  }

  crit_values  <- CV_shift(d=d,procedure="supwaldfixedb",param=bandw)
  testfixbsupw <- max(out)
  result       <- c(crit_values,testfixbsupw)
  names(result)<- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
