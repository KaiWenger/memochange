#' @title Wilcoxon long memory test for a single change in the mean of a long-memory time series.
#' @description This function performs a Wilcoxon type test for a change-in-mean that is robust under long memory. It applies a consistent estimator of the
#' long-run variance under long memory and uses a different normalization compared to a standard Wilcoxon test.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15}.
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param d integer that specifies the long-memory parameter.
#' @param tau integer that defines the search area, which is \code{[tau,1-tau]}. Default is \code{tau=0.15} as suggested by Andrews (1993).
#' @return Returns a numeric vector containing the test statistic and the corresponding critical values of the test.
#' @seealso \code{\link{snwilcoxon}}
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
#' # estimate the long-memory parameter of both series via local Whittle approach.
#' # The bandwidth to estimate d is chosen as T^0.65, which is usual in literature
#' d_est    <- LongMemoryTS::local.W(tseries, m=floor(1+T^0.65))$d
#' d_est2   <- LongMemoryTS::local.W(tseries2, m=floor(1+T^0.65))$d
#'
#' # perform the test on both time series
#' wilcoxonLM(tseries, d=d_est)
#' wilcoxonLM(tseries2, d=d_est2)
#' # For the series with no change in mean the test does not reject the null hypothesis
#' # of a constant mean across time at any reasonable significance level.
#' # For the series with a change in mean the test rejects the null hypothesis 
#' # at a 5% significance level.
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Dehling, H. and Rooch, A. and Taqqu, M. S. (2012): Non-Parametric Change-Point Tests for Long-Range Dependent Data. Scandinavian Journal of Statistics, 40, pp. 153-173.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

wilcoxonLM <- function(x,d,tau=0.15)
{
  if (any(is.na(x)))
    stop("x contains missing values")
  if(tau<=0 | tau>=1)
    stop("It must hold that 0<tau<1")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x) == 
      FALSE) 
    stop("x must be a univariate numeric vector")
  if(tau!=0.15)
    warning("Critical values are just implemented for tau=0.15")
  
  T   <- length(x)
  out <- vector("numeric",(T-1))
  for(k in 1:(T-1))
  {
    X1     <- x[1:k]
    X2     <- x[(k+1):T]
    X1_m   <- matrix(X1,length(X1),length(X2))
    X2_m   <- matrix(X2,length(X1),length(X2),byrow=TRUE)
    out[k] <- abs(sum((X1_m<=X2_m)-0.5))
  }

  crit_values   <- CV_shift(d=d,procedure="wilcoxonLM",param=0)
  testwilco     <- (1/T^(1.5+d))*max(out[(T*tau):(T*(1-tau))])
  result        <- c(crit_values,testwilco)
  names(result) <- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
