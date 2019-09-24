#' @title Self-normalized sup Wald test for a single change in the mean of a long-memory time series.
#' @description This function performs a sup Wald test for a change-in-mean that is robust under long memory. In contrast to a standard sup Wald test
#' it applies a self-normalization approach to estimate the long-run variance.
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
#' snwilcoxon(tseries,d=d_est)
#' snwilcoxon(tseries2,d=d_est2)
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Shao, X. (2011): A simple test of changes in mean in the possible presence of long-range dependence. Journal of Time Series Analysis, 32, pp. 598-606.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

snsupwald <- function(tseries,d,tau=0.15)
{

  n          <- length(tseries)
  grid       <- 1:n

  x_bar_1k   <- (cumsum(tseries)/grid)
  x_bar_k1n  <- (cumsum(tseries[n:1])/grid[1:n])[n:1]
  enumerator <- sqrt(n)*abs(x_bar_1k[-n]-x_bar_k1n[-1])

  S_1k                   <- t(matrix(rep(cumsum(tseries),n),n,n))
  S_1k                   <- S_1k-t(apply(matrix(rep(x_bar_1k,n),n,n),1,function(x){x*(1:n)}))
  S_1k[upper.tri(S_1k)]  <- 0
  sumS2_1k               <- rowSums(S_1k^2)
  S_k1n                  <- matrix(rep(tseries,n),n,n,byrow = TRUE)
  S_k1n[lower.tri(S_k1n)]<- 0
  S_k1n                  <- (apply(S_k1n,1,function(x){cumsum(x)}))
  helpf                  <- matrix(rep(x_bar_k1n[n:1],n),n,n,byrow=TRUE)#*(1:n)
  helpf                  <- helpf[,n:1]
  helpf2                 <- matrix(rep_len(c(1:n,1),n*n),n,n)
  helpf                  <- helpf2*helpf
  helpf[upper.tri(helpf)]<- 0
  sumS2_k1n              <- colSums((S_k1n-helpf)^2)
  denominator            <- n^(-1)*sqrt(sumS2_1k[-n]+sumS2_k1n[-1])

  crit_values   <- CV_shift(d=d,procedure="snsupwald",param=0)
  testsnsupwald <- max((enumerator/denominator)[(n*tau):(n*(1-tau))])
  result        <- c(crit_values,testsnsupwald)
  names(result) <- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
