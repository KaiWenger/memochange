#' @title Fixed-b sup Wald test for a single change in the mean of a long-memory time series.
#' @description This function performs a sup-Wald test on a change-in-mean, which is standardized by a non-parametric kernel-based long-run variance estimator.
#' Therefore, the test is robust under long-memory.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15} using the Bartlett kernel.
#' @param tseries the univariate numeric time series to be investigated.
#' @param d the long-memory parameter.
#' @param bandw the bandwidth parameter for the long-run variance estimator. It can take values in the range \code{bandw=[0.05,0.1,0.2]}. Default is
#' \code{bandw=0.1}, which is suggested by Iacone, Leybourne and Taylor (2014).
#' @param tau defines the search area, which is \code{[tau,1-tau]}. Default is \code{tau=0.15} as suggested by Andrews (1993).
#' @return Returns the test statistic and the corresponding critical values of the test.
#' @author Kai Wenger
#' @examples
#' library(fracdiff)
#' library(LongMemoryTS)
#'
#' n        <- 500
#' d        <- 0.2
#' tseries  <- fracdiff.sim(n=n,d=d)$series
#' d_est    <- local.W(tseries, m=floor(1+n^0.65))$d
#'
#' changep  <- c(rep(0,n/2),rep(1,n/2))
#' tseries2 <- tseries+changep
#' d_est2   <- local.W(tseries2, m=floor(1+n^0.65))$d
#'
#' fixbsupw(tseries,d=d_est)
#' fixbsupw(tseries2,d=d_est2)
#' @references
#' Iacone, F. and Leybourne, S. J. and Taylor, R. A. M. (2014): A fixed-b Test for a Break in Level at an unknown Time under Fractional Integration. Journal of Time Series Analysis, 35, pp. 40-54.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

fixbsupw <- function(tseries,d,bandw=0.1,tau=0.15)
{
  if(!(bandw %in% c(0.05,0.1,0.2))) stop("Use one of the following bandwidths: 0.05,0.1,0.2")

  n     <- length(tseries)
  m     <- bandw*n
  cc    <- t(c(0,1))
  out   <- c()

  for(k in ((n*tau):(n*(1-tau))))
  {
    X     <- cbind(rep(1,n),c(rep(0,k),rep(1,(n-k))))
    reg   <- lm(tseries~X-1)
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
