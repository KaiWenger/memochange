#' @title Simple test for change-in-mean under long memory
#' @description This function performs a CUSUM test on a change-in-mean that is robust under long memory. It is based on the fractionally differenced series where
#' the long-memory parameter is estimated by a consistent estimator.
#' The function returns the test statistic as well as the p-value of the test.
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param d integer that specifies the long-memory parameter.
#' @return Returns a numeric vector containing the test statistic and the p-value of the test.
#' @seealso \code{\link{CUSUMLM}}, \code{\link{CUSUMfixed}}
#' @author Kai Wenger
#' @examples
#' # set model parameters
#' T        <- 500
#' d        <- 0.2
#' 
#' set.seed(410)
#' 
#' # generate a fractionally integrated (long-memory) time series without a change in mean
#' tseries  <- fracdiff::fracdiff.sim(n=T, d=d)$series
#'
#' # generate a fractionally integrated (long-memory) time series 
#' # with a change in mean in the middle of the series
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
#' CUSUM_simple(tseries, d_est)
#' CUSUM_simple(tseries2, d_est2)
#' # For the series with no change in mean the test does not 
#' # reject the null hypothesis of a constant mean across time 
#' # at any reasonable significance level.
#' # For the series with a change in mean the test rejects the 
#' # null hypothesis at a 5% significance level.
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): A simple test on structural change in long-memory time series. Economics Letters, 136, pp. 90-94.
#' @export

CUSUM_simple    <- function(x,d)
{
  if (any(is.na(x)))
    stop("x contains missing values")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x) == 
      FALSE) 
    stop("x must be a univariate numeric vector")
  
  diff.series   <- fracdiff::diffseries(x,d=d)
  CUSUM_diff    <- strucchange::gefp(diff.series ~ 1, fit = stats::lm, vcov = sandwich::kernHAC)
  result        <- strucchange::sctest(CUSUM_diff)
  result        <- c(result$statistic, result$p.value)
  names(result) <- c("Teststatistic","p-value")
  return(round(result,3))
}
