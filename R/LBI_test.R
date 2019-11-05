#' @title Locally best invariant test against a change in persistence
#' @description This function performs the locally best invariant test against a change in persistence as suggested by Busetti and Taylor (2004). Under the null hypothesis the time series is I(0) throughout and
#' under the alternative a change from either I(1) to I(0) or I(0) to I(1) has occured. 
#' @details
#' The critical values of the tests vary with the sample size. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Busetti and Taylor (2004). These are, however, only valid for \code{tau=0.2}. 
#' In case that another value is chosen for \code{tau}, it is recommended to set \code{simu=1} which means that critical values are simulated based on the given data using M replications.
#' For a time series of length \code{T=100} and \code{M=10,000} replications this takes approximately five minutes with increasing duration for higher T or M.  
#' It should be noted, however, that M smaller than 10,000 make the results unreliable.
#'
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param trend whether the time series exhibits a trend, \code{"none"} implies no trend and \code{"linear"} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature. Note that if \code{T*tau<=1+as.numeric(trend=="linear")} the test statistic cannot be calculated.
#' @param statistic which type of test statistic should be used, \code{"mean"} corresponds to Hansen's (1991) mean score, \code{"max"} to Andrews' (1993) maximum statistic, and \code{"exp"} to Andrews and Ploberger's (1994) mean-exponential statistic.
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values should be simulated. Default is \code{M=10000}.
#' @return Returns a matrix that consists of test statistic and critical values (corresponding to \code{alpha=0.1,0.05,0.01}) for testing against a change from I(1) to I(0), I(0) to I(1), and against a change in an unknown direction.
#' @seealso \code{\link{cusum_test}}, \code{\link{LKSN_test}}, \code{\link{MR_test}}, \code{\link{ratio_test}}.
#' @author Janis Becker
#' @examples
#' set.seed(410)
#' 
#' # generate dummy-data
#' series <- c(rnorm(100), cumsum(rnorm(100)))
#' 
#' # test for a break in persistence
#' LBI_test(series, trend="none", statistic="mean")
#' @references
#' Busetti, F. and Taylor, R. (2004): Tests of stationarity against a change in persistence. Journal of Econometrics, 123, pp. 33-66.
#' @export
LBI_test<-function(x,trend=c("none","linear"),tau=0.2,statistic=c("mean","max","exp"),simu=0,M=10000)
{
  statistic<-match.arg(statistic,c("mean","max","exp"))
  trend<-match.arg(trend,c("none","linear"))
  if(tau<=0 | tau>=0.5)
    stop("It must hold that 0<tau<0.5")
  if (any(is.na(x)))
    stop("x contains missing values")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric vector")
  T<-length(x)
  f<-as.numeric(trend=="linear") 
  if ((T*tau)<=(f+1))
    stop("increase T*tau to guarantee that the test statistic can be calculated")
  if(tau!=0.2 & simu==0)
    warning("Note that the critical values stated are not valid for a tau different from 0.2")
  stat<-LBI(x=x,trend=trend,tau=tau)
  if(statistic=="mean") t_stats<-c(mean(stat$tstat1),mean(stat$tstat2),max(mean(stat$tstat1),mean(stat$tstat2)))
  if(statistic=="max") t_stats<-c(max(stat$tstat1),max(stat$tstat2),max(stat$tstat1,stat$tstat2))
  if(statistic=="exp") t_stats<-c(log(mean(exp(.5*stat$tstat1))),log(mean(exp(.5*(stat$tstat2)))),max(log(mean(exp(.5*stat$tstat1))),log(mean(exp(.5*(stat$tstat2))))))
  if(simu==1){Crit<-CV(x=x,statistic=statistic,trend=trend,type="LBI",M=M,tau=tau)}
  else{
    if(trend=="none" & statistic=="mean") Crit<-getCV()$cv_lbi_test[1:3,c(1,4,7)]
    if(trend=="none" & statistic=="exp") Crit<-getCV()$cv_lbi_test[1:3,c(2,5,8)]
    if(trend=="none" & statistic=="max") Crit<-getCV()$cv_lbi_test[1:3,c(3,6,9)]
    if(trend=="linear" & statistic=="mean") Crit<-getCV()$cv_lbi_test[4:6,c(1,4,7)]
    if(trend=="linear" & statistic=="exp") Crit<-getCV()$cv_lbi_test[4:6,c(2,5,8)]
    if(trend=="linear" & statistic=="max") Crit<-getCV()$cv_lbi_test[4:6,c(3,6,9)]
    Crit<-t(Crit)
  }
  result<-cbind(Crit,t_stats)
  colnames(result)<-c("90%","95%","99%","Teststatistic")
  rownames(result)<-c("Against change from I(0) to I(1)","Against change from I(1) to I(0)","Against change in unknown direction")
  return(result)
}


#' function to calculate sequence of LBI test statistics by Busetti and Taylor (2004). For internal use only.
#' @keywords internal
LBI<-function(x,trend,tau)
{
  T<-length(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  tstat1<-rep(NA,length(Ttau))
  tstat2<-rep(NA,length(Ttau))
  if(trend=="none"){p<-0}else{p<-1}
  tr<-(1:T)^p
  resi<-stats::lm(x~tr)$residuals
  var<-mean(resi^2)
  q<-1
  for(i in Ttau)
  {
    tstat1[q]<-1/(var*(T-i)^2)*sum((cumsum(rev(resi[(i+1):T])))^2)
    tstat2[q]<-1/(var*(i)^2)*sum((rev(cumsum(rev(resi)))[1:i])^2)
    q<-q+1
  }
  return(list(tstat1=tstat1,tstat2=tstat2))
}
