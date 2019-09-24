#' @title Self-normalized CUSUM tests for structural change under long memory.
#' @description This function performs a family of CUSUM tests for a change-in-mean that are robust under long memory. They apply non-parametric kernel-based
#' fixed-b and fixed-m long-run variance estimators in the denominator of the test statistics.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15} using the Bartlett kernel for the fixed-b tests or averaging the first m periodogram
#' ordinates (which corresponds to the Daniell kernel) for the fixed-m tests.
#' @param tseries the univariate numeric time series to be investigated.
#' @param d the long-memory parameter.
#' @param procedure specifies whether the CUSUM fixed-b or fixed-m type A or type B tests are used. It can be chosen between
#' \code{"CUSUMfixedb_typeA"}, \code{"CUSUMfixedb_typeB"}, \code{"CUSUMfixedm_typeA"}, and \code{"CUSUMfixedm_typeB"} (see Wenger, Leschinski (2019) for details).
#' @param bandw the bandwidth used for estimation of the long-run variance. For the fixed-b tests \code{b=[0.05,0.1,0.2,0.3,...,0.9,1]}, for the
#' fixed-m tests \code{m=[1,2,3,4,10,25,50,100,150,200]}. Recommended bandwidth by Wenger, Leschinski (2019) are \code{b=0.1} and \code{m=10}.
#' @param tau defines the search area, which is \code{[tau,1-tau]}. Default is \code{tau=0.15} as suggested by Andrews (1993).
#' @return Returns the test statistic and the corresponding critical values of the test.
#' @author Kai Wenger
#' @examples
#' library(fracdiff)
#' library(longmemo)
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
#' CUSUMfixed(tseries,d=d_est,procedure="CUSUMfixedb_typeA",bandw=0.1)
#' CUSUMfixed(tseries,d=d_est,procedure="CUSUMfixedb_typeB",bandw=0.1)
#' CUSUMfixed(tseries,d=d_est,procedure="CUSUMfixedm_typeA",bandw=10)
#' CUSUMfixed(tseries,d=d_est,procedure="CUSUMfixedm_typeB",bandw=10)
#'
#' CUSUMfixed(tseries2,d=d_est2,procedure="CUSUMfixedb_typeA",bandw=0.1)
#' CUSUMfixed(tseries2,d=d_est2,procedure="CUSUMfixedb_typeB",bandw=0.1)
#' CUSUMfixed(tseries2,d=d_est2,procedure="CUSUMfixedm_typeA",bandw=10)
#' CUSUMfixed(tseries2,d=d_est2,procedure="CUSUMfixedm_typeB",bandw=10)
#' @references
#' Wenger, K. and Leschinski, C. (2019): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Hualde, J. and Iacone, F. (2017): Fixed bandwidth asymptotics for the studentized mean of fractionally integrated processes. Economics Letters, 150, pp. 39-43.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

CUSUMfixed     <- function(tseries,d,procedure,bandw,tau=0.15)
{
  if(procedure   != "CUSUMfixedb_typeA" & procedure != "CUSUMfixedb_typeB"
     & procedure != "CUSUMfixedm_typeA" & procedure != "CUSUMfixedm_typeB")  stop("You misspecified which procedure to use")
  if((procedure == "CUSUMfixedb_typeA" | procedure == "CUSUMfixedb_typeB")
     & !(bandw %in% c(0.05,seq(0.1,1,0.1)))) stop("Use one of the following bandwidths: 0.05,0.1,0.2,0.3,...,0.9,1")
  if((procedure == "CUSUMfixedm_typeA" | procedure == "CUSUMfixedm_typeB")
     & !(bandw %in% c(1,2,3,4,10,25,50,100,150,200))) stop("Use one of the following bandwidths: 1,2,3,4,10,25,50,100,150,200")

  n            <- length(tseries)
  m            <- bandw*n
  u_hat        <- tseries-mean(tseries)
  cumulated    <- cumsum(u_hat)
  sigma        <- c()

  if(procedure=="CUSUMfixedb_typeA") sigma <- fb_longrun(u_hat,m)

  if(procedure=="CUSUMfixedb_typeB"){
    for(k in ((n*tau):(n*(1-tau))))
    {X      <- cbind(rep(1,n),c(rep(0,k),rep(1,(n-k))))
    reg    <- lm(tseries~X-1)
    u_hat2 <- unname(reg$residuals)
    sigm   <- fb_longrun(u_hat2,m)
    sigma  <- c(sigma,sigm)}}

  if(procedure=="CUSUMfixedm_typeA"){
    m            <- bandw
    sigma        <- (2*pi/m)*sum(per(u_hat)[2:(m+1)])}

  if(procedure=="CUSUMfixedm_typeB"){
    m            <- bandw
    for(k in ((n*tau):(n*(1-tau))))
    {First      <- c(rep(mean(tseries[1:k]),k),rep(0,(n-k)))
    Second     <- c(rep(0,k),rep(mean(tseries[(k+1):n]),(n-k)))
    X_tilde    <- tseries-First-Second
    sigma      <- c(sigma,(2*pi/m)*sum(per(X_tilde)[2:(m+1)]))}}

  crit_values             <- CV_shift(d=d,procedure=procedure,param=bandw)
  testCUSUMfixedb_typeB   <- max(abs(cumulated[((n*tau):(n*(1-tau)))]/(sqrt(n*sigma))))
  result                  <- c(crit_values,testCUSUMfixedb_typeB)
  names(result)           <- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
