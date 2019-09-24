#' @title Simple test for change-in-mean under long memory
#' @description This function performs a CUSUM test on a change-in-mean that is robust under long memory. It is based on the fractionally differenced series where
#' the long-memory parameter is estimated by a consistent estimator.
#' The function returns the test statistic as well as the p-value of the test.
#' @param tseries the univariate numeric time series to be investigated.
#' @param d the long-memory parameter.
#' @return Returns the test statistic and the p-value of the test.
#' @author Kai Wenger
#' @examples
#' library(fracdiff)
#' library(strucchange)
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
#' CUSUM_simple(tseries,d_est)
#' CUSUM_simple(tseries2,d_est2)
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): A simple test on structural change in long-memory time series. Economics Letters, 136, pp. 90-94.
#' @export

CUSUM_simple <- function(tseries,d)
{
  diff.series   <- diffseries(tseries,d=d)
  CUSUM_diff    <- gefp(diff.series ~ 1, fit = lm, vcov = kernHAC)
  result        <- sctest(CUSUM_diff)
  result        <- c(result$statistic, result$p.value)
  names(result) <- c("Teststatistic","p-value")
  return(round(result,3))
}
