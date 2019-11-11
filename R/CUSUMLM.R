#' @title CUSUM long memory test for a single change in the mean of a long-memory time series.
#' @description This function performs a modified CUSUM test for a change-in-mean that is robust under long memory. It replaces the standardization
#' as well as the long-run variance estimator compared to the standard CUSUM test.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15}.
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param d integer that specifies the long-memory parameter.
#' @param delta integer that determines the bandwidth that is used to estimate the constant \code{G} that approximates the short run dynamics of the time series at the origin.
#' The same bandwidth should be used that is applied to estimate \code{d} before. See Wenger, Leschinski, Sibbertsen (2018) for details.
#' @param tau integer that defines the search area, which is \code{[tau,1-tau]}. Default is \code{tau=0.15} as suggested by Andrews (1993).
#' @return Returns a numeric vector containing the test statistic and the corresponding critical values of the test.
#' @seealso \code{\link{CUSUMfixed}}, \code{\link{CUSUM_simple}}
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
#' CUSUMLM(tseries, delta=0.65, d=d_est)
#' CUSUMLM(tseries2, delta=0.65, d=d_est2)
#' # For the series with no change in mean the test does not 
#' # reject the null hypothesis of a constant mean across time 
#' # at any reasonable significance level.
#' # For the series with a change in mean the test rejects the 
#' # null hypothesis at a 1% significance level.
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Wang, L. (2008): Change-in-mean problem for long memory time series models with applications. Journal of Statistical Computation and Simulation, 78:7, pp. 653-668.
#'
#' Horvath, L. and Kokoszka, P. (1997): The effect of long-range dependence on change-point estimators. Journal of Statistical Planung and Inference, 64, pp. 57-81.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

CUSUMLM <- function(x,d,delta,tau=0.15)
{
  if (any(is.na(x)))
    stop("x contains missing values")
  if(tau<=0 | tau>=1)
    stop("It must hold that 0<tau<1")
  if(delta<0.6 | delta>0.8)
    stop("It must hold that 0.6<delta<0.8")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x) == 
      FALSE) 
    stop("x must be a univariate numeric vector")
  if(tau!=0.15)
    warning("Critical values are just implemented for tau=0.15")
  
  T            <- length(x)
  G            <- LongMemoryTS::G.hat(as.matrix(x), d, m=floor(1+T^(delta)))
  if(d!=0)
  {
    C2         <- gamma(1-2*d)*G*2*sin(pi*d)/(d*(1+2*d))
    enumerator <- c(T^(-1/2-d)/sqrt(C2))*cumsum(x-mean(x))
  }else
  {
    C2         <- 2*pi*G
    enumerator <- c(T^(-1/2-d)/sqrt(C2))*cumsum(x-mean(x))
  }


  crit_values  <- CV_shift(d=d,procedure="cusumlm",param=0)
  testCUSUMLM  <- max(abs(enumerator[(T*tau):(T*(1-tau))]))
  result       <- c(crit_values,testCUSUMLM)
  names(result)<- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
