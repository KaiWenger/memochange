#' @title Cusum-type test against a change in persistence
#' @description This function performs a cusum-type test for a change in persistence as suggested by Leybourne, Taylor, and Kim (2006) and extended by Sibbertsen and Kruse (2009).
#' Under the null hypothesis the time series is nonstationary throughout and under the alternative a change from nonstationary to stationary or vice versa has occured.
#' @details
#' Leybourne, Taylor, and Kim (2006) introduced a cusum-type test that is able to identify when time series
#' exhibit changes in persistence. Under the null
#' hypothesis, the series is throughout I(1), i.e. nonstationary. Under the alternative the series exhibits a
#' break either from I(0) to I(1) or vice versa. Sibbertsen and Kruse (2009) extended the test such that under the null hypothesis
#' the time series is I(d) throughout, with d>1/2 and under the alternative a change from I(d1) to I(d2), where d1<1/2 and 1/2<d2<3/2, or vice versa has occured.
#' While the test statistic remains the same, the critical values of the extended test change as they depend on the order of integration.
#' Furthermore, the procedure by SK integrates the series if d is estimated to be smaller than 1/2. This allows to overcome the problem of the approach by LKT which is that is has a degenerated limiting distribution when the series is stationary.
#' To determine the order of integration (the memory parameter d) the semiparametric estimator by Geweke and Porter-Hudak (1983) is used.
#'
#' The critical values of the tests vary with sample size and d. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Leybourne, Taylor, and Kim (2006) respectively the response curves by Sibbertsen and Kruse (2009).
#' These are, however, only valid for \code{tau=0.2} and \code{m=0}. 
#' In case that non-default values are chosen for \code{tau} or \code{m}, it is recommended to set \code{simu=1} which means that critical values are simulated based on the given data using M replications. 
#' Caution, for a time series of length \code{T=100} and \code{M=10,000} replications this takes approximately thirty minutes with increasing duration for higher T or M.  
#' It should be noted, however, that M smaller than 10,000 make the results unreliable.
#'
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param trend whether the time series exhibits a trend, \code{"none"} implies no trend and \code{"linear"} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature. Note that if \code{T*tau<=1+as.numeric(trend=="linear") + (m>3)*(m-3)} the test statistic cannot be calculated.
#' @param type which type of cusum test should be performed, \code{"LKT"} for the cusum test by Leybourne, Taylor, and Kim (2006) and \code{"SK"} for the extension by Sibbertsen and Kruse (2009). See details.
#' @param m Number of covariances used for the estimation of the long run variance. Default is \code{m=0}.
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation based on critical values for \code{tau=0.2}. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values should be simulated. Default is \code{M=10000}.
#' @return Returns a matrix that consists of test statistic and critical values for testing against a change from nonstationary to stationary, stationary to nonstationary, and against a change in an unknown direction.
#' @author Janis Becker
#' @examples
#' set.seed(410)
#'
#' # generate dummy-data
#' series <- c(rnorm(200), cumsum(rnorm(200)))
#' 
#' # test for a break in persistence
#' cusum_test(series, trend="none", type="SK")
#' @references
#' Leybourne, S., Kim, T., and Taylor, R. (2007): Cusum of squares-based tests for a change in persistence. Journal of Time Series Analysis, 28, pp. 408-433.
#'
#' Sibbertsen, P. and Kruse, R. (2009): Testing for a break in persistence under long-range dependencies. Journal of Time Series Analysis, 30, pp. 263-285.
#' @export

cusum_test<-function(x,trend=c("none","linear"),tau=0.2,type=c("LKT","SK"),m=0,simu=0,M=10000)
{
  type<-match.arg(type,c("LKT","SK"))
  trend<-match.arg(trend,c("none","linear"))
  if(tau<=0 | tau>=0.5)
    stop("It must hold that 0<tau<0.5")
  if (any(is.na(x)))
    stop("x contains missing values")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric vector")
  T<-length(x)
  f<-as.numeric(trend=="linear") + as.numeric(m>3)*(m-3)
  if ((T*tau)<=(f+1))
    stop("increase T*tau to guarantee that the test statistic can be calculated")
  if(tau!=0.2 & simu==0)
    warning("Note that the critical values stated are not valid for a tau different from 0.2")
  if(type=="SK"){
    d<-fracdiff::fdGPH(x,.8)$d
    if(d<.5){
      x<-cumsum(x)
      d<-d+1
    }
  }
  t_stats<-cusum(x=x,trend=trend,m=m,tau=tau)$tstat
  if(simu==1 & type=="SK") Crit<-CV(x=x,trend=trend,type="cusum_SK",m=m,M=M,d=d,tau=tau)
  if(simu==1 & type=="LKT") Crit<-CV(x=x,trend=trend,type="cusum",m=m,M=M,tau=tau)
  if(simu==0 & type=="LKT"){
    if(trend=="none") Crit<-getCV()$cv_cusum_test[1:5,]
    if(trend=="linear") Crit<-getCV()$cv_cusum_test[6:10,]
    T<-length(x)
    if(T<50) Crit<-Crit[1,]
    if(T>500) Crit<-Crit[5,]
    if(T>49 & T<501){
      if(min(abs(as.numeric(rownames(Crit))-T))==0){Crit<-Crit[which.min(rank(abs(as.numeric(rownames(Crit))-T))),]}
      else{
        Tdif<-as.numeric(rownames(Crit))-T
        if(Tdif[which.min(abs(Tdif))]<0){
          Crit<-Crit[which.min(abs(Tdif)):(which.min(abs(Tdif))+1),]
          Tdif<-Tdif[which.min(abs(Tdif)):(which.min(abs(Tdif))+1)]
        }
        else{
          Crit<-Crit[(which.min(abs(Tdif))-1):which.min(abs(Tdif)),]
          Tdif<-Tdif[(which.min(abs(Tdif))-1):which.min(abs(Tdif))]
        }
        Tdif<-abs(Tdif)
        Crit[1,]<-(sum(Tdif)-Tdif[1])/(sum(Tdif))*Crit[1,]+(sum(Tdif)-Tdif[2])/(sum(Tdif))*Crit[2,]
        Crit<-Crit[1,]
      }
    }
    Crit<-matrix(Crit,ncol=2,nrow=3)
    Crit[,2]<-rev(Crit[,2])
  }
  if(simu==0 & type=="SK"){
    dmat<-matrix(nrow=1,ncol=10)
    for(f in 0:9) dmat[1,f+1]<-d^f
    if(trend=="none") Crit<-cbind(getCV()$respcurve[1:3,]%*%t(dmat),rev(getCV()$respcurve[4:6,]%*%t(dmat)))
    if(trend=="linear") Crit<-cbind(getCV()$respcurve[7:9,]%*%t(dmat),rev(getCV()$respcurve[10:12,]%*%t(dmat)))
  }
  result<-cbind(t(Crit),t_stats)
  for(i in 1:2) result[i,1:3]=rev(result[i,1:3])
  colnames(result)<-c("90%","95%","99%","Teststatistic")
  rownames(result)<-c("Lower","Upper")
  return(result)
}


#' function to calculate sequence of cusum test statistics by Leybourne, Taylor, and Kim (2006). For internal use only
#' @keywords internal
cusum<-function(x,trend,m=0,tau)
{
  T<-length(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  q<-1
  if(trend=="none"){p<-0}
  else{p<-1}
  T1<-rep(NA,length(Ttau))
  T2<-rep(NA,length(Ttau))
  for(i in Ttau)
  {
    tr<-(1:i)^p
    resi_1<-stats::lm(x[1:i]~tr)$residuals
    resi_2<-stats::lm(rev(x)[1:i]~tr)$residuals
    if(m>0){
      index <- 1:m
      cov_1<- sapply(index, function(x) t(diff(resi_1[-c(1:x)])) %*%
                       diff(resi_1[-c((length(resi_1) - x + 1):length(resi_1))]))
      cov_2<- sapply(index, function(x) t(diff(resi_2[-c(1:x)])) %*%
                       diff(resi_2[-c((length(resi_2) - x + 1):length(resi_2))]))
      bartlett <- 1 - index/(m + 1)
      var_1 <- 1/i*sum(diff(resi_1)^2)+ 2/length(resi_1) * t(bartlett) %*% cov_1
      var_2 <- 1/i*sum(diff(resi_2)^2)+ 2/length(resi_2) * t(bartlett) %*% cov_2
    }
    else{var_1 <- 1/i*sum(diff(resi_1)^2);var_2 <- 1/i*sum(diff(resi_2)^2)}
    T1[q]<-(i^(-2)*sum(resi_1^2))/var_1
    T2[q]<-(i^(-2)*sum(resi_2^2))/var_2
    q<-q+1
  }
  return(list(tstat=min(T1)/min(T2),T1=T1,T2=T2))
}
