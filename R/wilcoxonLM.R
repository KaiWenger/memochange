#' @title Wilcoxon long memory test for a single change in the mean of a long-memory time series.
#' @description This function performs a Wilcoxon type test for a change-in-mean that is robust under long memory. It applies a consistent estimator of the
#' long-run variance under long memory and uses a different normalization compared to a standard Wilcoxon test.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15}.
#' @param tseries the univariate numeric time series to be investigated.
#' @param d the long-memory parameter.
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
#' wilcoxonLM(tseries,d=d_est)
#' wilcoxonLM(tseries2,d=d_est2)
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Dehling, H. and Rooch, A. and Taqqu, M. S. (2012): Non-Parametric Change-Point Tests for Long-Range Dependent Data. Scandinavian Journal of Statistics, 40, pp. 153-173.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

wilcoxonLM <- function(tseries,d,tau=0.15)
{
  n   <- length(tseries)
  out <- vector("numeric",(n-1))
  for(k in 1:(n-1))
  {
    X1     <- tseries[1:k]
    X2     <- tseries[(k+1):n]
    X1_m   <- matrix(X1,length(X1),length(X2))
    X2_m   <- matrix(X2,length(X1),length(X2),byrow=TRUE)
    out[k] <- abs(sum((X1_m<=X2_m)-0.5))
  }

  crit_values   <- CV_shift(d=d,procedure="wilcoxonLM",param=0)
  testwilco     <- (1/n^(1.5+d))*max(out[(n*tau):(n*(1-tau))])
  result        <- c(crit_values,testwilco)
  names(result) <- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
