#' @title Simulates persistence-break process
#' @description This function simulates a process that exhibits a break in persistence.
#' @param T length of the time series.
#' @param trend whether the time series exhibits a trend, \code{none} implies no trend and \code{linear} implies a linear trend.
#' @param tau break fraction, \code{T*tau} yields the break point. It needs to hold that \code{0<tau<1}.
#' @param tp trend parameter, \code{t*tp} yields the contribution of the trend component if \code{trend="linear"}.
#' @param d1 order of integration of the first part of the series.
#' @param d2 order of integration of the second part of the series.
#' @param mean mean of the series.
#' @param var variance of the innovations.
#' @return Returns a vector containing the simulated time series.
#' @author Janis Becker
#' @examples
#' library(LongMemoryTS)
#' series<- pb_sim(500,0.5,"none",d1=0.2,d2=0.8,mean=0,var=1)
#' ts.plot(series)
#' @export

pb_sim<-function(T,tau,trend=c("none","linear"),tp=1,d1,d2,mean,var){
  T1<-round(T*tau)
  T2<-round(T*(1-tau))
  if((T1+T2)>T){T2<-T2-1}
  if((T1+T2)<T){T2<-T2+1}
  series<-c(FI.sim(T1,q=1,rho=0,d=d1,var=var),FI.sim(T2,q=1,rho=0,d=d2,var=var))+mean
  if(trend=="linear") series<-series + 1:T*tp
  return(series)
}
