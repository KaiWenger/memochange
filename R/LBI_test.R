#' @title Locally best invariant test against a change in persistence
#' @description This function performs the locally best invariant test against a change in persistence as suggested by Busetti and Taylor (2004). Under the null hypothesis the time series is I(0) throughout and
#' under the alternative a change from either I(1) to I(0) or I(0) to I(1) has occured. 
#' @details
#' The critical values of the tests vary with the sample size. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Busetti and Taylor (2004). These are based on \code{tau=0.2}. If \code{simu=1},
#' the critical values are simulated based on the given data using M replications. Caution, for large M this might take a while,
#' small M, however, make the results unreliable.
#'
#' @param x the univariate numeric time series to be investigated.
#' @param trend whether the time series exhibits a trend, \code{none} implies no trend and \code{linear} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature.
#' @param statistic which type of test statistic should be used, \code{mean} corresponds to Hansen's (1991) mean score, \code{max} to Andrews' (1993) maximum statistic, and \code{exp} to Andrews and Ploberger's (1994) mean-exponential statistic
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values should be simulated. Default is \code{M=10000}.
#' @return Returns a matrix that consists of test statistic and critical values (corresponding to \code{alpha=0.1,0.05,0.01}) for testing against a change from I(1) to I(0), I(0) to I(1), and against a change in an unknown direction.
#' @author Janis Becker
#' @examples
#' series<- c(rnorm(200),cumsum(rnorm(200)))
#' LBI_test(series,trend="none",statistic="mean")
#' @references
#' Busetti, F. and Taylor, R. (2004): Tests of stationarity against a change in persistence. Journal of Econometrics, 123, pp. 33-66.
#' @export
LBI_test<-function(x,trend=c("none","linear"),tau=0.2,statistic=c("mean","max","exp"),simu=0,M=10000)
{
  statistic<-statistic[1]
  trend<-trend[1]
  if ((statistic %in% c("mean","max","exp")) == FALSE)
    stop("statistic must be one of mean, max, exp. See details.")
  if ((trend %in% c("none","linear")) == FALSE)
    stop("trend must be one of none, linear. See details.")
  if (any(is.na(x)))
    stop("missing values not allowed in time series")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric time series")
  T<-length(x)
  f<-as.numeric(trend=="linear") 
  if ((T*tau)<=(f+1))
    stop("increase T*tau to guarantee that the test statistic can be calculated")
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
  resi<-lm(x~tr)$residuals
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
