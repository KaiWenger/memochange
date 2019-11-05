#' @title Simulates persistence-break process
#' @description This function simulates a fractional white noise process that exhibits a break in persistence using \code{FI.sim} from the \code{LongMemoryTS} package. In the first part of the series the noise is integrated with order \code{d_1} and in the second part with order \code{d_2}.
#' @param T length of the time series.
#' @param trend whether the time series exhibits a trend, \code{"none"} implies no trend and \code{"linear"} implies a linear trend.
#' @param tau break fraction, \code{T*tau} yields the break point. It needs to hold that \code{0<tau<1}.
#' @param tp trend parameter, \code{t*tp} yields the contribution of the trend component if \code{trend="linear"}.
#' @param d1 order of integration of the first part of the series.
#' @param d2 order of integration of the second part of the series.
#' @param mean mean of the series. Default is \code{mean=0}.
#' @param var variance of the innovations. Default is \code{var=1}
#' @return Returns a vector containing the simulated time series.
#' @author Janis Becker
#' @examples
#' set.seed(410)
#' 
#' # generate persistence-break time series
#' series <- pb_sim(500, 0.5, "none", d1=0.2, d2=0.8, mean=0, var=1)
#' 
#' # plot generated series
#' stats::ts.plot(series)
#' @export

pb_sim<-function(T,tau,trend=c("none","linear"),tp=0,d1,d2,mean=0,var=1){
  trend<-match.arg(trend,c("none","linear"))
  if(tau<=0 | tau>=1)
    stop("It must hold that 0<tau<1")
  if(trend=="linear" & tp==0)
    warning("The series is supposed to exhibit a trend but tp is set to be zero")
  T1<-round(T*tau)
  T2<-round(T*(1-tau))
  if((T1+T2)>T){T2<-T2-1}
  if((T1+T2)<T){T2<-T2+1}
  S1=LongMemoryTS::FI.sim(T1,q=1,rho=0,d=d1,var=var,burnin=0)+mean
  S2=LongMemoryTS::FI.sim(T2,q=1,rho=0,d=d2,var=var,burnin=0)+S1[length(S1)]
  series=c(S1,S2)
  if(trend=="linear") series<-series + 1:T*tp
  return(series)
}
