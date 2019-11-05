#' @title LM test against a change in persistence
#' @description This function performs a LM-type test for a change in persistence as suggested by Martins and Rodrigues (2014).
#' Under the null hypothesis the memory parameter d is constant over the sample. Under the alternative
#' an increase or a decrease of the memory parameter has occured over time.
#' @details
#' The critical values of the tests vary with sample size and memory parameter d. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Martins and Rodrigues (2014). These are, however, only valid for \code{tau=0.2} and \code{serial=FALSE}. 
#' In case that non-default values are chosen for \code{tau} or \code{serial}, it is recommended to set \code{simu=1} which means that critical values are simulated based on the given data using M replications. 
#' Caution, for a time series of length \code{T=750} and \code{M=10,000} replications this takes approximately twelve hours with increasing duration for higher T or M.  
#' It should be noted, however, that M smaller than 10,000 make the results unreliable.
#'
#' @param x the univariate numeric vector to be investigated. Missing values are not allowed.
#' @param trend whether the time series exhibits a trend, \code{"none"} implies no trend and \code{"linear"} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature. Note that if \code{T*tau<=2} the test statistic cannot be calculated.
#' @param statistic which type of test statistic should be used, \code{"squared"} for the squared t-statistic and \code{"standard"} for the standard t-test. Default is \code{statistic="squared"}.
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values should be simulated. Default is \code{M=10000}.
#' @param serial boolean, indicating whether to account for serial correlation of the errors. Default is \code{serial=FALSE} implying no correction for serial correlation.
#' @return Returns a matrix that consists of test statistic and critical values (corresponding to \code{alpha=0.1,0.05,0.01}) for testing against an increase in memory, against a decrease in memory, and against a change in an unknown direction.
#' @seealso \code{\link{cusum_test}}, \code{\link{LBI_test}}, \code{\link{LKSN_test}}, \code{\link{ratio_test}}.
#' @author Janis Becker
#' @examples
#' set.seed(410)
#' 
#' # generate dummy-data
#' series <- c(rnorm(200), cumsum(rnorm(200)))
#' 
#' # test for a break in persistence
#' MR_test(series, trend="none", statistic="squared")
#' @references
#' Martins, L.. and Rodrigues, P. (2014): Testing for persistence change in fractionally integrated models: An application to world inflation rates. Computational Statistics and Data Analysis, 76, pp. 502-522.
#' @export

MR_test<-function(x,trend=c("none","linear"),tau=0.2,statistic=c("squared","standard"),simu=0,M=10000,serial=c(FALSE,TRUE))
{
  statistic<-match.arg(statistic,c("squared","standard"))
  trend<-match.arg(trend,c("none","linear"))
  serial <- serial[1]
  if ((serial %in% c(FALSE, TRUE)) == FALSE) 
    stop("serial must be one of FALSE, TRUE. See details.")
  if(tau<=0 | tau>=0.5)
    stop("It must hold that 0<tau<0.5")
  if (any(is.na(x)))
    stop("x contains missing values")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric vector")
  T<-length(x)
  if ((T*tau)<=2)
    stop("increase T*tau to guarantee that the test statistic can be calculated")
  if(tau!=0.2 & simu==0)
    warning("Note that the critical values stated are not valid for a tau different from 0.2")
  #calculate test statistic
  stat<-MR(x=x,trend=trend,serial=serial,tau=tau)
  if(statistic=="standard") t_stats<-c(min(stat$tstat1),min(stat$tstat2),min(stat$tstat1,stat$tstat2))
  if(statistic=="squared") t_stats<-c(max(stat$tstat1^2),max(stat$tstat2^2),max(stat$tstat1^2,stat$tstat2^2))
  #simulate critical values
  if(simu==1){Crit<-CV(x=x,trend=trend,type="MR",M=M,statistic=statistic,d=stat$d0,tau=tau,serial=serial)}
  #extract critical values
  else{
    if(statistic=="standard") Crit<-getCV()$cv_MR_test
    if(statistic=="squared") Crit<-getCV()$cv_MR_test_squared
    Crit<-Crit[abs(as.numeric(rownames(Crit))-stat$d0)<.1,]
    ddif<-abs(as.numeric(rownames(Crit))-stat$d0)
    ddif<-ddif[c(1,4)]
    Crit[1:3,]<-(sum(ddif)-ddif[1])/(sum(ddif))*Crit[1:3,]+(sum(ddif)-ddif[2])/(sum(ddif))*Crit[4:6,]
    Crit<-Crit[1:3,]
    T<-stat$T
    if(T<100) Crit<-Crit[,1:3]
    if(T>750) Crit<-Crit[,10:12]
    if(T>99 & T<751){
      if(min(abs(as.numeric(colnames(Crit))-T))==0){Crit<-Crit[,rank(abs(as.numeric(colnames(Crit))-T))<3]}
      else{
        Tdif<-as.numeric(colnames(Crit))-T
        if(Tdif[which.min(abs(Tdif))]<0){
          Crit<-Crit[,(which.min(abs(Tdif))):(which.min(abs(Tdif))+5)]
          Tdif<-Tdif[(which.min(abs(Tdif))):(which.min(abs(Tdif))+5)]
        }
        else{
          Crit<-Crit[,(which.min(abs(Tdif))-3):(which.min(abs(Tdif))+2)]
          Tdif<-Tdif[(which.min(abs(Tdif))-3):(which.min(abs(Tdif))+2)]
        }
        Tdif<-abs(Tdif)[c(1,4)]
        Crit[,1:3]<-(sum(Tdif)-Tdif[1])/(sum(Tdif))*Crit[,1:3]+(sum(Tdif)-Tdif[2])/(sum(Tdif))*Crit[,4:6]
        Crit<-Crit[,1:3]
      }
    }
  }
  #arrange results
  result<-cbind(Crit,t_stats)
  if(statistic=="standard") for(j in 1:3) result[j,1:3]<-rev(result[j,1:3])
  colnames(result)<-c("90%","95%","99%","Teststatistic")
  rownames(result)<-c("Against increase in memory","Against decrease in memory","Against change in unknown direction")
  return(result)
}

#' function to calculate sequence of test statistics by Martins and Rodrigues (2014). For internal use only
#' @keywords internal
MR<-function(x,trend,serial,tau){
  T<-length(x)
  d<-fracdiff::fdGPH(x,.8)$d
  while(d>.5){
    x<-diff(x)
    d<-d-1
    T<-length(x)
  }
  d0<-stats::coef(longmemo::WhittleEst(x))-.5
  dif_series<-LongMemoryTS::fdiff(x,d=d0)
  pi<-LongMemoryTS::fdiff(rep(1,T),d=d0)
  mu<-stats::coef(stats::lm(dif_series~pi+0))
  x<-dif_series-pi*mu
  if(serial==TRUE){
    x<-stats::residuals(forecast::auto.arima(x,max.p=12*(T/100^(1/4)),max.q=0))
  }
  x_rev<-rev(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  T1<-rep(NA,length(Ttau))
  T2<-rep(NA,length(Ttau))
  q<-1
  index<-0
  for(i in Ttau)
  {
    xstar<-0
    for(l in 1:(i-1)) xstar[l]=sum(rev(x[1:l])/1:l) 
    model<-stats::lm(x[2:i]~xstar+0)
    if(serial==TRUE){
      T1[q]<-stats::coef(model)[1]/sqrt(sandwich::vcovHC(model)[1,1])
    }
    else{
      T1[q]<-sum(x[2:i]*xstar)/(sqrt(1/(i-2)*sum(stats::residuals(model)^2))*sqrt(sum(xstar^2)))
    }
    
    xstar<-0
    for(l in 1:(i-1)) xstar[l]=sum(rev(x_rev[1:l])/1:l) 
    model<-stats::lm(x_rev[2:i]~xstar+0)
    if(serial==TRUE){
      T2[q]<-stats::coef(model)[1]/sqrt(sandwich::vcovHC(model)[1,1])
    }
    else{
      T2[q]<-sum(x_rev[2:i]*xstar)/(sqrt(1/(i-2)*sum(stats::residuals(model)^2))*sqrt(sum(xstar^2)))
    }
    
    index[q]<-i
    q<-q+1
  }
  return(list(tstat1=T1,tstat2=T2,d0=d0,i=index,T=T))
}
