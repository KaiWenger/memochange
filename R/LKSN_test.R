#' @title DF-type test against a change in persistence
#' @description This function performs the DF-type test against a change in persistence as suggested by Leybourne, Kim, Smith, and Newbold (2003). Under the null hypothesis the time series is I(1) throughout and
#' under the alternative a change from either I(1) to I(0) or I(0) to I(1) has occured. 
#' @details
#' The critical values of the tests vary with the sample size. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Leybourne, Kim, Smith, and Newbold (2003). These are based on \code{tau=0.2}. If \code{simu=1},
#' the critical values are simulated based on the given data using M replications. Caution, for large M this might take a while,
#' small M, however, make the results unreliable.
#'
#' @param x the univariate numeric time series to be investigated.
#' @param trend whether the time series exhibits a trend, \code{none} implies no trend and \code{linear} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature.
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values should be simulated. Default is \code{M=10000}.
#' @return Returns a matrix that consists of test statistic and critical values (corresponding to \code{alpha=0.1,0.05,0.01}) for testing against a change from I(1) to I(0), I(0) to I(1), and against a change in an unknown direction.
#' @author Janis Becker
#' @examples
#' library(urca)
#' series<- c(rnorm(200),cumsum(rnorm(200)))
#' LKSN_test(series,trend="none")
#' @references
#' Leybourne, S., Kim, T., Smith, V., and Newbold, P. (2003): Tests for a change in persistence against the null of difference-stationarity. Econometrics Journal, 6, pp. 291-311.
#' @export
LKSN_test<-function(x,trend=c("none","linear"),tau=0.2,simu=0,M=10000)
{
  trend<-trend[1]
  if ((trend %in% c("none","linear")) == FALSE)
    stop("trend must be one of none, linear. See details.")
  if (any(is.na(x)))
    stop("missing values not allowed in time series")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric time series")
  if ((T*tau)<=10)
    stop("T*tau needs to be at least 11 to guarantee that the test statistic can be calculated")
  stat<-LKSN(x=x,trend=trend,tau=tau)
  t_stats<-c(min(stat$tstat1),min(stat$tstat2),min(stat$tstat1,stat$tstat2))
  if(simu==1){Crit<-CV(x=x,trend=trend,type="LKSN",M=M,tau=tau)}
  else{
    if(trend=="none") Crit<-getCV()$cv_LKSN_test[1:3,]
    if(trend=="linear") Crit<-getCV()$cv_LKSN_test[4:6,]
    T<-length(x)
    if(T<100) Crit<-Crit[,1:2]
    if(T>1000) Crit<-Crit[,9:10]
    if(T>99 & T<1001){
      if(min(abs(as.numeric(colnames(Crit))-T))==0){Crit<-Crit[,rank(abs(as.numeric(colnames(Crit))-T))<2]}
      else{
        Tdif<-as.numeric(colnames(Crit))-T
        if(Tdif[which.min(abs(Tdif))]<0){
          Crit<-Crit[,(which.min(abs(Tdif))):(which.min(abs(Tdif))+3)]
          Tdif<-Tdif[(which.min(abs(Tdif))):(which.min(abs(Tdif))+3)]
        }
        else{
          Crit<-Crit[,(which.min(abs(Tdif))-2):(which.min(abs(Tdif))+1)]
          Tdif<-Tdif[(which.min(abs(Tdif))-2):(which.min(abs(Tdif))+1)]
        }
        Tdif<-abs(Tdif)[c(1,3)]
        Crit[,1:2]<-(sum(Tdif)-Tdif[1])/(sum(Tdif))*Crit[,1:2]+(sum(Tdif)-Tdif[2])/(sum(Tdif))*Crit[,3:4]
        Crit<-Crit[,1:2]
      }
    }
  }
  result<-cbind(Crit,t_stats)
  colnames(result)<-c("90%","95%","Teststatistic")
  rownames(result)<-c("Against change from I(0) to I(1)","Against change from I(1) to I(0)","Against change in unknown direction")
  return(result)
}


#' function to calculate sequence of LKSN test statistics. For internal use only
#' @keywords internal
LKSN<-function(x,trend,tau)
{
  T<-length(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  T1<-rep(NA,length(Ttau))
  T2<-rep(NA,length(Ttau))
  q<-1
  if(trend=="linear"){
    for(i in Ttau)
    {
      T1[q]<-ur.ers(x[1:i],model="trend")@teststat
      T2[q]<-ur.ers(rev(x)[1:i],model="trend")@teststat
      q<-q+1
    }
  }
  else{
    for(i in Ttau)
    {
      T1[q]<-ur.ers(x[1:i])@teststat
      T2[q]<-ur.ers(rev(x)[1:i])@teststat
      q<-q+1
    }
  }
  return(list(tstat1=T1,tstat2=T2))
}

