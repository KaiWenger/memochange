#' @title Self-normalized Wilcoxon test for a single change in the mean of a long-memory time series.
#' @description This function performs a Wilcoxon test for a change-in-mean that is robust under long memory. In contrast to a standard Wilcoxon test
#' it applies a self-normalization approach to estimate the long-run variance.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15}. Furthermore, it is assumed that we have a 1st-order Hermite process. For details see Betken (2016).
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
#' Betken, A. (2016): Testing for change-points in long-range dependent time series by means of a self-normalized wilcoxon test. Journal of Time Series Analysis, 37, pp. 785-908.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

snwilcoxon <- function(tseries,d,tau=0.15)
{
  n          <- length(tseries)
  R          <- vector("numeric",n)
  grid       <- 1:n

  X1<-X2     <- tseries
  X1_m       <- matrix(X1,length(X1),length(X2))
  X2_m       <- matrix(X2,length(X1),length(X2),byrow=TRUE)
  R          <- colSums((X1_m<=X2_m))

  x_bar_1k   <- (cumsum(R)/grid)
  x_bar_k1n  <- (cumsum(R[n:1])/grid[1:n])[n:1]

  S_1k                   <- t(matrix(rep(cumsum(R),n),n,n))
  S_1k                   <- S_1k-t(apply(matrix(rep(x_bar_1k,n),n,n),1,function(x){x*(1:n)}))
  S_1k[upper.tri(S_1k)]  <- 0
  sumS2_1k               <- rowSums(S_1k^2)
  S_k1n                  <- matrix(rep(R,n),n,n,byrow = TRUE)
  S_k1n[lower.tri(S_k1n)]<- 0
  S_k1n                  <- (apply(S_k1n,1,function(x){cumsum(x)}))
  helpf                  <- matrix(rep(x_bar_k1n[n:1],n),n,n,byrow=TRUE)#*(1:n)
  helpf                  <- helpf[,n:1]
  helpf2                 <- matrix(rep_len(c(1:n,1),n*n),n,n)
  helpf                  <- helpf2*helpf
  helpf[upper.tri(helpf)]<- 0
  sumS2_k1n              <- colSums((S_k1n-helpf)^2)
  denominator            <- n^(-1/2)*sqrt(sumS2_1k[-n]+sumS2_k1n[-1])
  enumerator             <- ((1:n)*((cumsum(R)/(1:n))-mean(R)))[-n]

  testsnwilcoxon <- max(abs(enumerator/denominator))
  crit_values    <- CV_shift(d=d,procedure="snwilcox",param=0)
  result         <- c(crit_values,testsnwilcoxon)
  names(result)  <- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
