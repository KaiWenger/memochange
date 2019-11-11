#' @title Self-normalized CUSUM tests for structural change under long memory.
#' @description This function performs a family of CUSUM tests for a change-in-mean that are robust under long memory. They apply non-parametric kernel-based
#' fixed-b and fixed-m long-run variance estimators in the denominator of the test statistics.
#' The function returns the test statistic as well as critical values.
#' @details
#' Note that the critical values are generated for \code{tau=0.15} using the Bartlett kernel for the fixed-b tests or averaging the first m periodogram
#' ordinates (which corresponds to the Daniell kernel) for the fixed-m tests.
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param d integer that specifies the long-memory parameter.
#' @param procedure string that specifies whether the CUSUM fixed-b or fixed-m type A or type B tests are used. It can be chosen between
#' \code{"CUSUMfixedb_typeA"}, \code{"CUSUMfixedb_typeB"}, \code{"CUSUMfixedm_typeA"}, and \code{"CUSUMfixedm_typeB"} (see Wenger, Leschinski (2019) for details).
#' @param bandw integer that determines the bandwidth used for estimation of the long-run variance. For the fixed-b tests \code{b=[0.05,0.1,0.2,0.3,...,0.9,1]}, for the
#' fixed-m tests \code{m=[1,2,3,4,10,25,50,100,150,200]}. Recommended bandwidth by Wenger, Leschinski (2019) are \code{b=0.1} and \code{m=10}.
#' @param tau integer that defines the search area, which is \code{[tau,1-tau]}. Default is \code{tau=0.15} as suggested by Andrews (1993).
#' @return Returns a numeric vector containing the test statistic and the corresponding critical values of the test.
#' @seealso \code{\link{CUSUMLM}}, \code{\link{CUSUM_simple}}, \code{\link{fixbsupw}}
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
#' # perform the different types of the test on both time series
#' CUSUMfixed(tseries, d=d_est, procedure="CUSUMfixedb_typeA", bandw=0.1)
#' CUSUMfixed(tseries, d=d_est, procedure="CUSUMfixedb_typeB", bandw=0.1)
#' CUSUMfixed(tseries, d=d_est, procedure="CUSUMfixedm_typeA", bandw=10)
#' CUSUMfixed(tseries, d=d_est, procedure="CUSUMfixedm_typeB", bandw=10)
#'
#' CUSUMfixed(tseries2, d=d_est2, procedure="CUSUMfixedb_typeA", bandw=0.1)
#' CUSUMfixed(tseries2, d=d_est2, procedure="CUSUMfixedb_typeB", bandw=0.1)
#' CUSUMfixed(tseries2, d=d_est2, procedure="CUSUMfixedm_typeA", bandw=10)
#' CUSUMfixed(tseries2, d=d_est2, procedure="CUSUMfixedm_typeB", bandw=10)
#' # For the series with no change in mean all tests do not reject 
#' # the null hypothesis of a constant mean across time at 
#' # any reasonable significance level.
#' # For the series with a change in mean all tests reject the 
#' # null hypothesis at a 1% significance level.
#' @references
#' Wenger, K. and Leschinski, C. (2019): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Hualde, J. and Iacone, F. (2017): Fixed bandwidth asymptotics for the studentized mean of fractionally integrated processes. Economics Letters, 150, pp. 39-43.
#'
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#' @export

CUSUMfixed     <- function(x,d,procedure,bandw,tau=0.15)
{
  if(procedure   != "CUSUMfixedb_typeA" & procedure != "CUSUMfixedb_typeB"
     & procedure != "CUSUMfixedm_typeA" & procedure != "CUSUMfixedm_typeB")  stop("You misspecified which procedure to use")
  if((procedure == "CUSUMfixedb_typeA" | procedure == "CUSUMfixedb_typeB")
     & !(bandw %in% c(0.05,seq(0.1,1,0.1)))) stop("Use one of the following bandwidths: 0.05,0.1,0.2,0.3,...,0.9,1")
  if((procedure == "CUSUMfixedm_typeA" | procedure == "CUSUMfixedm_typeB")
     & !(bandw %in% c(1,2,3,4,10,25,50,100,150,200))) stop("Use one of the following bandwidths: 1,2,3,4,10,25,50,100,150,200")
  if (any(is.na(x)))
    stop("x contains missing values")
  if(tau<=0 | tau>=1)
    stop("It must hold that 0<tau<1")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x) == 
      FALSE) 
    stop("x must be a univariate numeric vector")
  if(tau!=0.15)
    warning("Critical values are just implemented for tau=0.15")
  
  T            <- length(x)
  m            <- bandw*T
  u_hat        <- x-mean(x)
  cumulated    <- cumsum(u_hat)
  sigma        <- c()

  if(procedure=="CUSUMfixedb_typeA") sigma <- fb_longrun(u_hat,m)

  if(procedure=="CUSUMfixedb_typeB"){
    for(k in ((T*tau):(T*(1-tau))))
    {X      <- cbind(rep(1,T),c(rep(0,round(k)),rep(1,round(T-k))))
    reg     <- stats::lm(x~X-1)
    u_hat2  <- unname(reg$residuals)
    sigm    <- fb_longrun(u_hat2,m)
    sigma   <- c(sigma,sigm)}}

  if(procedure=="CUSUMfixedm_typeA"){
    m            <- bandw
    sigma        <- (2*pi/m)*sum(longmemo::per(u_hat)[2:(m+1)])}

  if(procedure=="CUSUMfixedm_typeB"){
    m            <- bandw
    for(k in ((T*tau):(T*(1-tau))))
    {First      <- c(rep(mean(x[1:k]),round(k)),rep(0,round(T-k)))
    Second      <- c(rep(0,round(k)),rep(mean(x[(k+1):T]),round(T-k)))
    X_tilde     <- x-First-Second
    sigma       <- c(sigma,(2*pi/m)*sum(longmemo::per(X_tilde)[2:(m+1)]))}}

  crit_values             <- CV_shift(d=d,procedure=procedure,param=bandw)
  testCUSUMfixedb_typeB   <- max(abs(cumulated[((T*tau):(T*(1-tau)))]/(sqrt(T*sigma))))
  result                  <- c(crit_values,testCUSUMfixedb_typeB)
  names(result)           <- c("90%","95%","99%","Teststatistic")
  return(round(result,3))
}
