#' @title LM test against a change in persistence
#' @description This function performs a LM-type test for a change in persistence as suggested by Martins and Rodrigues (2014).
#' Under the null hypothesis the memory parameter d is constant over the sample. Under the alternative
#' an increase or a decrease of the memory parameter has occured over time.
#' @details
#' The critical values of the tests vary with sample size and memory parameter d. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Martins and Rodrigues (2014). These are based on \code{tau=0.2}.
#' If \code{simu=1}, the critical values are simulated based on the given data using M replications. Caution, for large M this might take a while,
#' small M, however, make the results unreliable.
#'
#' @param x the univariate numeric time series to be investigated.
#' @param trend whether the time series exhibits a trend, \code{none} implies no trend and \code{linear} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature.
#' @param statistic which type of test statistic should be used, \code{standard} for the standard t-test and \code{squared} for the squared test statistic.
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values should be simulated. Default is \code{M=10000}.
#' @param p order of the AR-model used to account for serial correlation of the errors. Default is \code{p=0}.
#' @param twostep boolean, indicating whether the two step procedure or the one step procedure should be used. Default is \code{twostep=FALSE}.
#' @return Returns a matrix that consists of test statistic and critical values (corresponding to \code{alpha=0.1,0.05,0.01}) for testing against an increase in memory, against a decrease in memory, and against a change in an unknown direction.
#' @author Janis Becker
#' @examples
#' library(fracdiff)
#' library(LongMemoryTS)
#' series<- c(rnorm(200),cumsum(rnorm(200)))
#' MR_test(series,trend="none",statistic="standard",twostep=FALSE)
#' @references
#' Martins, L.. and Rodrigues, P. (2014): Testing for persistence change in fractionally integrated models: An application to world inflation rates. Computational Statistics and Data Analysis, 76, pp. 502-522.
#' @export

MR_test<-function(x,trend=c("none","linear"),tau=0.2,statistic=c("standard","squared"),simu=0,M=10000,p=0,twostep=c(FALSE,TRUE))
{
  trend<-trend[1]
  statistic<-statistic[1]
  twostep<-twostep[1]
  if ((twostep %in% c(FALSE,TRUE)) == FALSE)
    stop("twostep must be one of FALSE, TRUE. See details.")
  if ((trend %in% c("none","linear")) == FALSE)
    stop("trend must be one of none, linear. See details.")
  if ((statistic %in% c("standard","squared")) == FALSE)
    stop("statistic must be one of standard, squared. See details.")
  if (any(is.na(x)))
    stop("missing values not allowed in time series")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric time series")
  T<-length(x)
  if ((T*tau)<=2+ as.numeric(p>1)*(p))
    stop("increase T*tau to guarantee that the test statistic can be calculated")
  #calculate test statistic
  stat<-MR(x=x,trend=trend,p=p,twostep=twostep,tau=tau)
  if(statistic=="standard") t_stats<-c(min(stat$tstat1),min(stat$tstat2),min(stat$tstat1,stat$tstat2))
  if(statistic=="squared") t_stats<-c(max(stat$tstat1^2),max(stat$tstat2^2),max(stat$tstat1^2,stat$tstat2^2))
  #simulate critical values
  if(simu==1){Crit<-CV(x=x,trend=trend,type="MR",M=M,statistic=statistic,d=stat$d0,tau=tau,p=p,twostep=twostep)}
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
  if(statistic=="standard") colnames(result)<-c("99%","95%","90%","Teststatistic")
  if(statistic=="squared") colnames(result)<-c("90%","95%","99%","Teststatistic")
  rownames(result)<-c("Against increase in memory","Against decrease in memory","Against change in unknown direction")
  return(result)
}

#' function to calculate sequence of test statistics by Martins and Rodrigues (2014). For internal use only
#' @keywords internal
MR<-function(x,trend,p,twostep,tau){
  T<-length(x)
  d<-fdGPH(x,.8)$d
  while(d>.5){
    x<-diff(x)
    d<-d-1
    T<-length(x)
  }
  d0<-coef(fracdiff(x,nar=p))[1]
  x<-fdiff(x,d0)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  if(trend=="none"){z<-rep(1,T)}else{z<-cbind(1,1:T)}
  zeta<-1
  for(i in 2:T){zeta[i]<-(i-1-d0)/i*zeta[i-1]}
  if(trend=="none"){z<-cumsum(zeta*z)}else{z<-apply(zeta*z,2,cumsum)}
  x<-residuals(lm(x~z))
  if(twostep==TRUE & p>0) x<-residuals(arima(x,order=c(p,0,0)));p<-0
  x_rev<-rev(x)
  T1<-rep(NA,length(Ttau))
  T2<-rep(NA,length(Ttau))
  q<-1
  index<-0
  for(i in Ttau)
  {
    xstar<-0
    for(l in 1:i){xstar[l]<-sum(rev(x[1:l])/1:l)}
    xstar<-xstar[-i]
    xx<-x[-1]
    if(p>0 |twostep==TRUE){
      if(p>0){
        xx_emb<-embed(as.matrix(xx),dimension = p+1)
        model<-lm(xx_emb[1:(i-1-p),1]~xstar[-p]+xx_emb[1:(i-1-p),-1]+0)
      }
      else{model<-lm(xx[1:(i-1)]~xstar+0)}
      T1[q]<-coef(model)[1]/sqrt(vcovHC(model)[1,1])
    }
    else{
      model<-lm(xx[1:(i-1)]~xstar+0)
      T1[q]<-coef(model)[1]/(sqrt(vcov(model)[1,1]))
    }
    xstar<-0
    for(l in 1:i){xstar[l]<-sum(rev(x_rev[1:l])/1:l)}
    xstar<-xstar[-i]
    xx<-x_rev[-1]
    if(p>0 |twostep==TRUE){
      if(p>0){
        xx_emb<-embed(as.matrix(xx),dimension = p+1)
        model<-lm(xx_emb[1:(i-1-p),1]~xstar[-p]+xx_emb[1:(i-1-p),-1]+0)
      }
      else{model<-lm(xx[1:(i-1)]~xstar+0)}
      T2[q]<-coef(model)[1]/sqrt(vcovHC(model)[1,1])
    }
    else{
      model<-lm(xx[1:(i-1)]~xstar+0)
      T2[q]<-coef(model)[1]/(sqrt(vcov(model)[1,1]))
    }
    index[q]<-i
    q<-q+1
  }
  return(list(tstat1=T1,tstat2=T2,d0=d0,i=index,T=T))
}
