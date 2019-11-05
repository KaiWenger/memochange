#' @title DF-type test against a change in persistence
#' @description This function performs the DF-type test against a change in persistence as suggested by Leybourne, Kim, Smith, and Newbold (2003). Under the null hypothesis the time series is I(1) throughout and
#' under the alternative a change from either I(1) to I(0) or I(0) to I(1) has occured. 
#' @details
#' The critical values of the tests vary with the sample size. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Leybourne, Kim, Smith, and Newbold (2003). These are, however, only valid for \code{tau=0.2} and \code{lmax=0}. 
#' In case that non-default values are chosen for \code{tau} or \code{lmax}, it is recommended to set \code{simu=1} which means that critical values are simulated based on the given data using M replications. 
#' Caution, for a time series of length \code{T=100} and \code{M=10,000} replications this takes approximately thirty minutes with increasing duration for higher T or M.  
#' It should be noted, however, that M smaller than 10,000 make the results unreliable.
#'
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param trend whether the time series exhibits a trend, \code{"none"} implies no trend and \code{"linear"} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature. Note that if \code{T*tau<11} the test statistic cannot be calculated.
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values are simulated. Default is \code{M=10000}.
#' @param lmax Maximum number of lagged differences to be included in the test regression. Default is \code{lmax=0}. Note that small sample critical values might differ for \code{lmax>0} so that simulation is recommended in this case.
#' @return Returns a matrix that consists of test statistic and critical values (corresponding to \code{alpha=0.1,0.05}) for testing against a change from I(1) to I(0), I(0) to I(1), and against a change in an unknown direction.
#' @seealso \code{\link{cusum_test}}, \code{\link{LBI_test}}, \code{\link{MR_test}}, \code{\link{ratio_test}}.
#' @author Janis Becker
#' @examples
#' set.seed(410)
#' 
#' # generate dummy-data
#' series <- c(rnorm(200), cumsum(rnorm(200)))
#' 
#' # test for a break in persistence
#' LKSN_test(series, trend="none")
#' @references
#' Leybourne, S., Kim, T., Smith, V., and Newbold, P. (2003): Tests for a change in persistence against the null of difference-stationarity. Econometrics Journal, 6, pp. 291-311.
#' @export
LKSN_test<-function(x,trend=c("none","linear"),tau=0.2,lmax=0,simu=0,M=10000)
{
  trend<-match.arg(trend,c("none","linear"))
  if(tau<=0 | tau>=0.5)
    stop("It must hold that 0<tau<0.5")
  if (any(is.na(x)))
    stop("x contains missing values")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric vector")
  T<-length(x)
  if ((T*tau)<11)
    stop("T*tau needs to be at least 11 to guarantee that the test statistic can be calculated")
  if(tau!=0.2 & simu==0)
    warning("Note that the critical values stated are not valid for a tau different from 0.2")
  if(lmax>0 & simu==0)
    warning("Note that the small sample critical values stated might be different for lmax different from 0")
  stat<-LKSN(x=x,trend=trend,tau=tau,lmax=lmax)
  t_stats<-c(min(stat$tstat1),min(stat$tstat2),min(stat$tstat1,stat$tstat2))
  if(simu==1){Crit<-CV(x=x,trend=trend,type="LKSN",M=M,tau=tau,lmax=lmax)}
  else{
    if(trend=="none") Crit<-getCV()$cv_LKSN_test[1:3,]
    if(trend=="linear") Crit<-getCV()$cv_LKSN_test[4:6,]
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

#' Unit root test by Elliot et al. (1996). For internal use only
#' @keywords internal
ers_test=function (y, trend, lag.max,T) 
{
  lag.max <- lag.max + 1
  nobs <- length(y)
  if (trend == "none") {
    ahat <- 1 - 25/T
    ya <- c(y[1], y[2:nobs] - ahat * y[1:(nobs - 1)])
    za1 <- c(1, rep(1 - ahat, nobs - 1))
    yd.reg <- summary(stats::lm(ya ~ -1 + za1))
    yd <- y - stats::coef(yd.reg)[1]
  }
  else if (trend == "linear") {
    ahat <- 1 - 25/T
    ya <- c(y[1], y[2:nobs] - ahat * y[1:(nobs - 1)])
    za1 <- c(1, rep(1 - ahat, nobs - 1))
    trd <- 1:nobs
    za2 <- c(1, trd[2:nobs] - ahat * trd[1:(nobs - 1)])
    yd.reg <- summary(stats::lm(ya ~ -1 + za1 + za2))
    yd <- y - stats::coef(yd.reg)[1] - stats::coef(yd.reg)[2] * trd
  }
  yd.l <- yd[1:(nobs - 1)]
  yd.diff <- diff(yd)
  if (lag.max > 1) {
    yd.dlags <- stats::embed(diff(yd), lag.max)[, -1]
    data.dfgls <- data.frame(cbind(yd.diff[-(1:(lag.max - 
                                                  1))], yd.l[-(1:(lag.max - 1))], yd.dlags))
    colnames(data.dfgls) <- c("yd.diff", "yd.lag", paste("yd.diff.lag", 
                                                         1:(lag.max - 1), sep = ""))
    dfgls.form <- stats::formula(paste("yd.diff ~ -1 + ", paste(colnames(data.dfgls)[-1], 
                                                         collapse = " + ")))
  }
  else if (lag.max <= 1) {
    data.dfgls <- data.frame(cbind(yd.diff, yd.l))
    colnames(data.dfgls) <- c("yd.diff", "yd.lag")
    dfgls.form <- stats::formula("yd.diff ~ -1 + yd.lag")
  }
  dfgls.reg <- summary(stats::lm(dfgls.form, data = data.dfgls))
  teststat <- stats::coef(dfgls.reg)[1, 3]
  return(teststat)
}

#' function to calculate sequence of LKSN test statistics. For internal use only
#' @keywords internal
LKSN<-function(x,trend,tau,lmax)
{
  T<-length(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  T1<-rep(NA,length(Ttau))
  T2<-rep(NA,length(Ttau))
  q<-1
  if(trend=="linear"){
    for(i in Ttau)
    {
      T1[q]<-ers_test(x[1:i],trend="linear",T=T,lag.max=lmax)
      T2[q]<-ers_test(rev(x)[1:i],trend="linear",T=T,lag.max=lmax)
      q<-q+1
    }
  }
  else{
    for(i in Ttau)
    {
      T1[q]<-ers_test(x[1:i],trend="none",T=T,lag.max=lmax)
      T2[q]<-ers_test(rev(x)[1:i],trend="none",T=T,lag.max=lmax)
      q<-q+1
    }
  }
  return(list(tstat1=T1,tstat2=T2))
}

