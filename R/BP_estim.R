#' @title Breakpoint estimators for a change in persistence
#' @description This function estimates the location where the investigated time series exhibits a break in persistence. It requires
#' knowledge of the direction of the break, i.e. an increase or decrease in persistence. This
#' needs to be determined beforehand using one of the various persistence change tests provided in this package.
#' @details
#' The estimators BT and LKSN are only consistent for changes from I(0) to I(1) or vice versa, the LKT estimator is consistent for changes from stationary to nonstationary memory or vice versa (cf. also Sibbertsen and Kruse (2009)), and the MR estimator is consistent for changes in d in general.
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param trend whether the time series exhibits a trend, \code{"none"} implies no trend and \code{"linear"} implies a linear trend.
#' @param tau the function searches in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature. Note that if \code{type="BT"} and \code{T*tau<=2+ as.numeric(trend=="linear")}, \code{type="LT"} and \code{T*tau<=2+ as.numeric(trend=="linear") + (m>3)*(m-3)}, \code{type="LKSN"} and \code{T*tau<=10}, or \code{type="MR"} and \code{T*tau<=2+(p>1)*p}the break point cannot be found.
#' @param type which type of break point estimator should be used, \code{"LKSN"} for the estimator by Leybourne, Kim, Smith, and Newbold (2003), \code{"BT"} for the estimator by Busetti and Taylor (2004), \code{"LKT"} for the estimator by Leybourne, Kim, and Taylor (2006),
#'  and \code{MR} for the estimator by Martins and Rodrigues (2014). See details.
#' @param direction direction of the change in persistence, \code{"01"} implies an increase in persistence over time and \code{"10"} a decrease. See details.
#' @param d_estim which estimator should be used to determine the order of integration in the two regimes, \code{"GPH"} corresponds to the estimator by Geweke and Porter-Hudak (1983) and \code{"ELW"} corresponds to the exact local Whittle estimator by Shimotsu and Phillips (2005).
#' @param d_bw bandwidth used for estimating the order of integration d. Default is \code{d_bw=0.7}. Note that the estimation of the memory parameter can only be performed for 0<d_bw<=1 and it is even recommended to choose 0.5<=d_bw<=0.8 as otherwise the estimators might be inconsistent.")
#' @param m Number of covariances used for the estimation of the long run variance when considering the LKT estimator. Default is \code{m=0}.
#' @param serial boolean, indicating whether to account for serial correlation of the errors when considering the MR estimator. Default is \code{serial=FALSE} implying no correction for serial correlation.
#' @return Returns a list that contains break point, estimate of the order of integration in the two regimes (the memory parameter d) and standard deviation of this estimate.
#' @author Janis Becker
#' @seealso \code{\link{cusum_test}}, \code{\link{LKSN_test}}, \code{\link{MR_test}}, \code{\link{ratio_test}}.
#' @examples
#' set.seed(410)
#' 
#' # generate dummy-data
#' series <- c(rnorm(200), cumsum(rnorm(200)))
#' 
#' # estimate the break point
#' BP_estim(series, trend="none", type="BT", direction="01", d_estim="ELW")
#' @references
#' Leybourne, S., Kim, T., Smith, V., and Newbold, P. (2003): Tests for a change in persistence against the null of difference-stationarity. Econometrics Journal, 6, pp. 291-311.
#'
#' Busetti, F. and Taylor, R. (2004): Tests of stationarity against a change in persistence. Journal of Econometrics, 123, pp. 33-66.
#'
#' Leybourne, S., Kim, T., and Taylor, R. (2007): Cusum of squares-based tests for a change in persistence. Journal of Time Series Analysis, 28, pp. 408-433.
#'
#' Martins, L.. and Rodrigues, P. (2014): Testing for persistence change in fractionally integrated models: An application to world inflation rates Cusum of squares-based tests for a change in persistence. Computational Statistics and Data Analysis, 76, pp. 502-522.
#' @export


BP_estim<-function(x,trend=c("none","linear"),tau=0.2,type=c("BT","LKT","LKSN","MR"),direction=c("01","10"),d_estim=c("ELW","GPH"),d_bw=0.7,m=0,serial=c(FALSE,TRUE))
{
  trend<-match.arg(trend,c("none","linear"))
  type<-match.arg(type,c("BT","LKT","LKSN","MR"))
  direction<-match.arg(direction,c("01","10"))
  d_estim<-match.arg(d_estim,c("ELW","GPH"))
  serial <- serial[1]
  if ((serial %in% c(FALSE, TRUE)) == FALSE) 
    stop("serial must be one of FALSE, TRUE. See details.")
  if(tau<=0 | tau>=0.5)
    stop("It must hold that 0<tau<0.5.")
  if(d_bw<=0 | d_bw>1)
    stop("It needs to hold that 0<d_bw<=1 as otherwise the estimation cannot be performed.")
  if(d_bw<0.5 | d_bw>0.8)
    warning("It is recommended to choose 0.5<=d_bw<=0.8 as otherwise the estimators might be inconsistent.")
  if (any(is.na(x)))
    stop("x contains missing values.")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric vector.")
  T<-length(x)
  f<-as.numeric(trend=="linear")
  if ((T*tau)<(f+2) & type=="BT")
    stop("T*tau needs to be at least 2 when type is 'BT'.")
  if ((T*tau)<(f+2+ as.numeric(m>3)*(m-3)) & type=="LKT")
    stop("increase T*tau to guarantee that the break point can be found.")
  if ((T*tau)<=10& type=="LKSN")
    stop("T*tau needs to be at least 11 to guarantee that the break point can be found when type is 'LKSN'")
  if ((T*tau)<=2 & type=="MR")
    stop("increase T*tau to guarantee that the break point can be found")
  if(trend=="none"){b<-0}else{b<-1}
  if(type=="BT"){
    Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
    bp<-rep(NA,length(Ttau))
    q<-1
    for(i in Ttau)
    {
      trend_1<-(1:i)^b
      trend_2<-((i+1):T)^b
      resi_1<-stats::lm(x[1:i]~trend_1)$residuals
      resi_2<-stats::lm(x[(i+1):T]~trend_2)$residuals
      bp[q]<-(T-i)^(-2)*sum(resi_2^2)/(i^(-2)*sum(resi_1^2))
      q<-q+1
    }
    if(direction=="01"){breakpoint<-which.max(bp)+floor(tau*T)+1}
    else{breakpoint<-which.min(bp)+floor(tau*T)+1}
  }
  if(type=="LKT"){
    if(direction=="01"){breakpoint<-which.min(cusum(x,trend=trend,m=m,tau=tau)$T1)+floor(tau*T)+1}
    else{breakpoint<-which.min(rev(cusum(x,trend=trend,tau=tau,m=m)$T2))+floor(tau*T)+1}
  }
  if(type=="LKSN"){
    bp<-LKSN(x,trend=trend,tau=tau)
    if(direction=="01"){breakpoint<-which.min(bp$tstat1)+floor(tau*T)+1}
    else{breakpoint<-which.min(rev(bp$tstat2))+floor(tau*T)+1}
  }
  if(type=="MR"){
    stat<-MR(x,trend=trend,serial=serial,tau=tau)
    if(direction=="01"){breakpoint<-stat$i[which.min(stat$tstat1)]+1}
    else{breakpoint<-length(x)-stat$i[which.min(stat$tstat2)]+1}
  }
  if(d_estim=="GPH"){
    d1<-fracdiff::fdGPH(x[1:(breakpoint-1)],d_bw)
    d2<-fracdiff::fdGPH(x[breakpoint:T],d_bw)
    return(list(Breakpoint=breakpoint,d_1=d1$d,sd_1=d1$sd.reg,d_2=d2$d,sd_2=d2$sd.reg))
  }
  else{
    d1<-LongMemoryTS::ELW(x[1:(breakpoint-1)],floor(1+length(x[1:(breakpoint-1)])^d_bw))
    d2<-LongMemoryTS::ELW(x[breakpoint:T],floor(1+length(x[breakpoint:T])^d_bw))
    return(list(Breakpoint=breakpoint,d_1=d1$d,sd_1=d1$s.e.,d_2=d2$d,sd_2=d2$s.e.))
  }
}
