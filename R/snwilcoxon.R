#' @title Self-normalized Wilcoxon test for a single change in the mean of a long-memory time series.
#' @description This function performs a Wilcoxon test for a change-in-mean that is robust under long memory. In contrast to a standard Wilcoxon test
#' it applies a self-normalization approach to estimate the long-run variance.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15}. Furthermore, it is assumed that we have a 1st-order Hermite process. For details see Betken (2016).
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param d integer that specifies the long-memory parameter.
#' @param tau integer that defines the search area, which is \code{[tau,1-tau]}. Default is \code{tau=0.15} as suggested by Andrews (1993).
#' @return Returns a numeric vector containing the test statistic and the corresponding critical values of the test.
#' @seealso \code{\link{wilcoxonLM}}, \code{\link{snsupwald}}
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
#' snwilcoxon(tseries, d=d_est)
#' snwilcoxon(tseries2, d=d_est2)
#' # For the series with no change in mean the test does not reject the null hypothesis 
#' # of a constant mean across time at any reasonable significance level.
#' # For the series with a change in mean the test rejects the null hypothesis 
#' # at a 1% significance level.
#' @references
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Betken, A. (2016): Testing for change-points in long-range dependent time series by means of a self-normalized wilcoxon test. Journal of Time Series Analysis, 37, pp. 785-908.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

snwilcoxon <- function(x,d,tau=0.15)
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
  
  T          <- length(x)
  R          <- vector("numeric",T)
  grid       <- 1:T

  X1<-X2     <- x
  X1_m       <- matrix(X1,length(X1),length(X2))
  X2_m       <- matrix(X2,length(X1),length(X2),byrow=TRUE)
  R          <- colSums((X1_m<=X2_m))

  x_bar_1k   <- (cumsum(R)/grid)
  x_bar_k1n  <- (cumsum(R[T:1])/grid[1:T])[T:1]

  S_1k                   <- t(matrix(rep(cumsum(R),T),T,T))
  S_1k                   <- S_1k-t(apply(matrix(rep(x_bar_1k,T),T,T),1,function(x){x*(1:T)}))
  S_1k[upper.tri(S_1k)]  <- 0
  sumS2_1k               <- rowSums(S_1k^2)
  S_k1n                  <- matrix(rep(R,T),T,T,byrow = TRUE)
  S_k1n[lower.tri(S_k1n)]<- 0
  S_k1n                  <- (apply(S_k1n,1,function(x){cumsum(x)}))
  helpf                  <- matrix(rep(x_bar_k1n[T:1],T),T,T,byrow=TRUE)
  helpf                  <- helpf[,T:1]
  helpf2                 <- matrix(rep_len(c(1:T,1),T*T),T,T)
  helpf                  <- helpf2*helpf
  helpf[upper.tri(helpf)]<- 0
  sumS2_k1n              <- colSums((S_k1n-helpf)^2)
  denominator            <- T^(-1/2)*sqrt(sumS2_1k[-T]+sumS2_k1n[-1])
  enumerator             <- ((1:T)*((cumsum(R)/(1:T))-mean(R)))[-T]

  testsnwilcoxon <- max(abs(enumerator/denominator))
  crit_values    <- CV_shift(d=d,procedure="snwilcox",param=0)
  result         <- c(crit_values,testsnwilcoxon)
  names(result)  <- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
