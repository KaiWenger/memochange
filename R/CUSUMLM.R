#' @title CUSUM long memory test for a single change in the mean of a long-memory time series.
#' @description This function performs a modified CUSUM test for a change-in-mean that is robust under long memory. It replaces the standardization
#' as well as the long-run variance estimator compared to the standard CUSUM test.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15}.
#' @param tseries the univariate numeric time series to be investigated.
#' @param d the long-memory parameter.
#' @param delta the bandwidth that is used to estimate the constant \code{G} that approximates the short run dynamics of the time series at the origin.
#' The same bandwidth should be used that is applied to estimate \code{d} before. See Wenger, Leschinski, Sibbertsen (2018) for details.
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
#' CUSUMLM(tseries,delta=0.65,d=d_est)
#' CUSUMLM(tseries2,delta=0.65,d=d_est2)
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Wang, L. (2008): Change-in-mean problem for long memory time series models with applications. Journal of Statistical Computation and Simulation, 78:7, pp. 653-668.
#'
#' Horvath, L. and Kokoszka, P. (1997): The effect of long-range dependence on change-point estimators. Journal of Statistical Planung and Inference, 64, pp. 57-81.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

CUSUMLM <- function(tseries,d,delta,tau=0.15)
{
  n            <- length(tseries)
  G            <- G.hat(as.matrix(tseries), d, m=floor(1+n^(delta)))
  if(d!=0)
  {
    C2         <- gamma(1-2*d)*G*2*sin(pi*d)/(d*(1+2*d))
    enumerator <- c(n^(-1/2-d)/sqrt(C2))*cumsum(tseries-mean(tseries))
  }else
  {
    C2         <- 2*pi*G
    enumerator <- c(n^(-1/2-d)/sqrt(C2))*cumsum(tseries-mean(tseries))
  }


  crit_values  <- CV_shift(d=d,procedure="cusumlm",param=0)
  testCUSUMLM  <- max(abs(enumerator[(n*tau):(n*(1-tau))]))
  result       <- c(crit_values,testCUSUMLM)
  names(result)<- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
