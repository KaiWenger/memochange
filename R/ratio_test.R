#' @title Ratio-based test against a change in persistence
#' @description This function performs a ratio-based test against a change in persistence. Under the null hypothesis the time series is I(0) throughout and
#' under the alternative a change from either I(1) to I(0) or I(0) to I(1) has occured. 
#' @details
#' Busetti and Taylor (2004) (BT) introduced a test that is able to identify when time series
#' exhibit changes in persistence. Under the null
#' hypothesis, the series is constant I(0), i.e. stationary. Under the alternative the series exhibits a
#' break either from I(0) to I(1) or I(1) to I(0). As the test is oversized for weakly dependent time series,
#' Leybourne and Taylor (2004) (LT) standardized the test statistic by an estimate of the long run variance using m lags.
#' Another problem is that constant I(1) processes are neither covered under the null nor the alternative.
#' Here, the test often rejects the null although no change in persistence occured.
#' Harvey, Leybourne, and Taylor (2006) (HLT) introduced a modification where they multiply the test statistic by
#' a unit root test. This allows the test statistic to have the same critical values under both constant I(0) and constant I(1).
#' It should be noted, however, that only the critical values are identical, the distribution is highly irregular.
#'
#' The critical values of the tests vary with the sample size. If \code{simu=0}, the critical values provided
#' are based on linear interpolation of the critical values simulated by Harvey, Leybourne, and Taylor (2006). These are based on \code{tau=0.2}. If \code{simu=1},
#' the critical values are simulated based on the given setup using M replications. Caution, for large M this might take a while,
#' small M, however, make the results unreliable.
#'
#' @param x the univariate numeric time series to be investigated.
#' @param trend whether the time series exhibits a trend, \code{none} implies no trend and \code{linear} implies a linear trend.
#' @param tau the function tests in the interval \code{[T*tau,T*(1-tau)]} for a break in persistence with T being the length of the time series. It must hold that \code{0<tau<0.5}, default is \code{tau=0.2} as commonly used in the literature.
#' @param type which type of ratio test should be performed, \code{BT} for the standard ratio test by Busetti and Taylor (2004), \code{LT} for the modified ratio test by Leybourne and Taylor (2004), and \code{HLT} respectively \code{HLTmin} are the modified tests by Harvey, Leybourne, and Taylor (2006). See details.
#' @param statistic which type of test statistic should be used, \code{mean} corresponds to Hansen's (1991) mean score, \code{max} to Andrews' (1993) maximum statistic, and \code{exp} to Andrews and Ploberger's (1994) mean-exponential statistic
#' @param m Number of covariances used for the estimation of the long run variance if \code{LT} is considered. Default is \code{m=0}.
#' @param z Number of polynomials used if \code{HLT} or \code{HLTmin} are considered. Default is \code{z=9}.
#' @param simu whether critical values should be simulated or interpolated, \code{simu=1} means simulation, \code{simu=0} means interpolation based on critical values for \code{tau=0.2}. See details. Default is \code{simu=0}.
#' @param M number of replications in case critical values should be simulated. Default is \code{M=10000}.
#' @return Returns a matrix that consists of test statistic and critical values (corresponding to \code{alpha=0.1,0.05,0.01}) for testing against a change from I(1) to I(0), I(0) to I(1), and against a change in an unknown direction.
#' @author Janis Becker
#' @examples
#' series<- c(rnorm(200),cumsum(rnorm(200)))
#' ratio_test(series,trend="none",type="BT",statistic="mean")
#' @references
#' Busetti, F. and Taylor, R. (2004): Tests of stationarity against a change in persistence. Journal of Econometrics, 123, pp. 33-66.
#'
#' Leybourne, S. and Taylor, R. (2004): On tests for changes in persistence. Economics letters, 84, pp. 107-115.
#'
#' Harvey, D., Leybourne, S. and Taylor, R. (2006): Modified tests for a change in persistence. Journal of Econometrics, 134, pp. 441-469.
#' @export

ratio_test<-function(x,trend=c("none","linear"),tau=0.2,statistic=c("mean","max","exp"),type=c("BT","LT","HLT","HLTmin"),m=0,z=9,simu=0,M=10000)
{
  statistic<-statistic[1]
  trend<-trend[1]
  type<-type[1]
  if ((statistic %in% c("mean","max","exp")) == FALSE)
    stop("statistic must be one of mean, max, exp. See details.")
  if ((trend %in% c("none","linear")) == FALSE)
    stop("trend must be one of none, linear. See details.")
  if ((type %in% c("BT","LT","HLT","HLTmin")) == FALSE)
    stop("type must be one of BT, LT, HLT, HLTmin. See details.")
  if (any(is.na(x)))
    stop("missing values not allowed in time series")
  if (mode(x) %in% ("numeric") == FALSE | is.vector(x)==FALSE)
    stop("x must be a univariate numeric time series")
  T<-length(x)
  f<-as.numeric(trend=="linear")
  if ((T*tau)<=(1+f) & type%in% c("BT","HLT"))
    stop("increase T*tau to guarantee that the test statistic can be calculated")
  if ((T*tau)<=(1+f+as.numeric(m>3)*(m-3)) & type%in% c("LT"))
    stop("increase T*tau to guarantee that the test statistic can be calculated")
  if ((T*tau)<=(z+1) & type%in% c("HLTmin"))
    stop("increase T*tau to guarantee that the test statistic can be calculated")
  #calculate sequence of test statistics
  if(type=="LT") stat<-LT(x=x,trend=trend,m=m,tau=tau)
  if(type=="BT") stat<-BT(x=x,trend=trend,tau=tau)
  if(type=="HLT"){stat<-BT(x=x,trend=trend,tau=tau);WT<-HLT(x=x,trend=trend,z=z)}
  if(type=="HLTmin"){
    save<-HLTmin(x=x,trend=trend,z=z,tau=tau)
    stat<-save$tstat
    WT<-c(min(save$WT1),min(save$WT2),min(save$WT1,save$WT2))
  }
  #calculate function over test statistic depending on the statistic argument
  if(statistic=="mean"){
    t_stats<-c(mean(stat),mean(stat^(-1)),max(mean(stat),mean(stat^(-1))))
    if(type=="HLT"){
      if(trend=="none") t_stats<-t(exp(t(getCV()$b[1:3,c(1,4,7)])*-1*WT)*t_stats)
      if(trend=="linear") t_stats<-t(exp(t(getCV()$b[4:6,c(1,4,7)])*-1*WT)*t_stats)
    }
    if(type=="HLTmin"){
      if(trend=="none") t_stats<-t(exp(t(getCV()$b_min[1:3,c(1,4,7)])*-1*WT)*t_stats)
      if(trend=="linear") t_stats<-t(exp(t(getCV()$b_min[4:6,c(1,4,7)])*-1*WT)*t_stats)
    }
  }
  if(statistic=="max"){
    t_stats<-c(max(stat),max(stat^(-1)),max(stat,stat^(-1)))
    if(type=="HLT"){
      if(trend=="none") t_stats<-t(exp(t(getCV()$b[1:3,c(3,6,9)])*-1*WT)*t_stats)
      if(trend=="linear") t_stats<-t(exp(t(getCV()$b[4:6,c(3,6,9)])*-1*WT)*t_stats)
    }
    if(type=="HLTmin"){
      if(trend=="none") t_stats<-t(exp(t(getCV()$b_min[1:3,c(3,6,9)])*-1*WT)*t_stats)
      if(trend=="linear") t_stats<-t(exp(t(getCV()$b_min[4:6,c(3,6,9)])*-1*WT)*t_stats)
    }
  }
  if(statistic=="exp"){
    t_stats<-c(log(mean(exp(.5*stat))),log(mean(exp(.5*(stat^(-1))))),max(log(mean(exp(.5*stat))),log(mean(exp(.5*(stat^(-1)))))))
    if(type=="HLT"){
      if(trend=="none") t_stats<-t(exp(t(getCV()$b[1:3,c(2,5,8)])*-1*WT)*t_stats)
      if(trend=="linear") t_stats<-t(exp(t(getCV()$b[4:6,c(2,5,8)])*-1*WT)*t_stats)
    }
    if(type=="HLTmin"){
      if(trend=="none") t_stats<-t(exp(t(getCV()$b_min[1:3,c(2,5,8)])*-1*WT)*t_stats)
      if(trend=="linear") t_stats<-t(exp(t(getCV()$b_min[4:6,c(2,5,8)])*-1*WT)*t_stats)
    }
  }
  #generate critical values
  if(simu==1){Crit<-CV(x=x,statistic=statistic,trend=trend,type=type,m=m,M=M,tau=tau)}
  #extract interpolated critical values
  else{
    if(trend=="none" & statistic=="mean") Crit<-getCV()$cv_ratio_test[1:18,c(1,4,7)]
    if(trend=="none" & statistic=="exp") Crit<-getCV()$cv_ratio_test[1:18,c(2,5,8)]
    if(trend=="none" & statistic=="max") Crit<-getCV()$cv_ratio_test[1:18,c(3,6,9)]
    if(trend=="linear" & statistic=="mean") Crit<-getCV()$cv_ratio_test[19:36,c(1,4,7)]
    if(trend=="linear" & statistic=="exp") Crit<-getCV()$cv_ratio_test[19:36,c(2,5,8)]
    if(trend=="linear" & statistic=="max") Crit<-getCV()$cv_ratio_test[19:36,c(3,6,9)]
    T<-length(x)
    if(T<100) Crit<-Crit[1:3,]
    if(T>500) Crit<-Crit[16:18,]
    if(T>99 & T<501){
      if(min(abs(as.numeric(rownames(Crit))-T))==0){Crit<-Crit[rank(abs(as.numeric(rownames(Crit))-T))<3,]}
      else{
        Tdif<-as.numeric(rownames(Crit))-T
        if(Tdif[which.min(abs(Tdif))]<0){
          Crit<-Crit[(which.min(abs(Tdif))):(which.min(abs(Tdif))+5),]
          Tdif<-Tdif[(which.min(abs(Tdif))):(which.min(abs(Tdif))+5)]
        }
        else{
          Crit<-Crit[(which.min(abs(Tdif))-3):(which.min(abs(Tdif))+2),]
          Tdif<-Tdif[(which.min(abs(Tdif))-3):(which.min(abs(Tdif))+2)]
        }
        Tdif<-abs(Tdif)[c(1,4)]
        Crit[1:3,]<-(sum(Tdif)-Tdif[1])/(sum(Tdif))*Crit[1:3,]+(sum(Tdif)-Tdif[2])/(sum(Tdif))*Crit[4:6,]
        Crit<-Crit[1:3,]
      }
    }
    Crit<-t(Crit)
  }
  #arrange results
  if(type=="HLT" | type=="HLTmin"){
    result<-cbind(Crit,t(t_stats))
    colnames(result)<-c("90%","95%","99%","Teststatistic 90%", "Teststatistic 95%", "Teststatistic 99%")
  }
  else{
    result<-cbind(Crit,t_stats)
    colnames(result)<-c("90%","95%","99%","Teststatistic")
  }
  rownames(result)<-c("Against change from I(0) to I(1)","Against change from I(1) to I(0)","Against change in unknown direction")

  return(result)
}

#' function to calculate sequence of test statistics by Busetti and Taylor (2004). For internal use only.
#' @keywords internal
BT<-function(x,trend,tau)
{
  T<-length(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  tstat<-rep(NA,length(Ttau))
  q<-1
  if(trend=="none"){p<-0}else{p<-1}
  for(i in Ttau)
  {
    trend_1<-(1:i)^p
    trend_2<-((i+1):T)^p
    resi_1<-lm(x[1:i]~trend_1)$residuals
    resi_2<-lm(x[(i+1):T]~trend_2)$residuals
    tstat[q]<-(T-i)^(-2)*sum(cumsum(resi_2)^2)/(i^(-2)*sum(cumsum(resi_1)^2))
    q<-q+1
  }
  return(tstat)
}

#' function to calculate sequence of test statistics by Leybourne and Taylor (2004). For internal use only.
#' @keywords internal
LT<-function(x,trend,m,tau)
{
  T<-length(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  tstat<-rep(NA,length(Ttau))
  q<-1
  if(trend=="none"){p<-0}else{p<-1}
  for(i in Ttau)
  {
    trend_1<-(1:i)^p
    trend_2<-((i+1):T)^p
    resi_1<-lm(x[1:i]~trend_1)$residuals
    resi_2<-lm(x[(i+1):T]~trend_2)$residuals
    if(m>0){
      index <- 1:m
      cov_1<- sapply(index, function(x) t(resi_1[-c(1:x)]) %*%
                       resi_1[-c((length(resi_1) - x + 1):length(resi_1))])
      cov_2<- sapply(index, function(x) t(resi_2[-c(1:x)]) %*%
                       resi_2[-c((length(resi_2) - x + 1):length(resi_2))])
      bartlett <- 1 - index/(m + 1)
      var_1 <- mean(resi_1^2)+ 2/i * t(bartlett) %*% cov_1
      var_2 <- mean(resi_2^2)+ 2/(T-i) * t(bartlett) %*% cov_2
    }
    else{var_1 <- mean(resi_1^2);var_2 <- mean(resi_2^2)}
    tstat[q]<-var_1/var_2*(T-i)^(-2)*sum(cumsum(resi_2)^2)/(i^(-2)*sum(cumsum(resi_1)^2))
    q<-q+1
  }
  return(tstat)
}

#' function to calculate wald_test. For internal use only
#' @keywords internal
wald_test<-function(Sigma, b, Terms)
{
  w <- length(Terms)
  H0 <- rep(0, w)
  L <- matrix(rep(0, length(b) * w), ncol = length(b))
  for (i in 1:w) L[i, Terms[i]] <- 1
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))),
                            sep = ""), names(b))
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L),tol=1e-100)
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  structure(list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0,
                 L = L, result = res), class = "wald.test")
}

#' function to calculate unit root test to adjust test statistic as suggested by Harvey, Leybourne, and Taylor (2006).
#' @keywords internal
HLT<-function(x,trend,z)
{
  T<-length(x)
  tr<-matrix(nrow=T,ncol=z);for(f in 1:z){tr[,f]<-(1:T)^f}
  model<-lm(x~tr)
  if(trend=="none"){p<-0}else{p<-1}
  WT<-wald_test(vcov(model),coef(model),Terms=(p+2):(z+1))$result$chi2[1]
  return(WT/T)
}

#' function to calculate sequence of minimum test statistics by Harvey, Leybourne, and Taylor (2006). For internal use only
#' @keywords internal
HLTmin<-function(x,trend,z,tau)
{
  T<-length(x)
  Ttau<-(floor(T*tau)):(ceiling(T*(1-tau)))
  q<-1
  if(trend=="none"){p<-0}else{p<-1}
  tstat<-rep(NA,length(Ttau))
  WT1<-rep(NA,length(Ttau))
  WT2<-rep(NA,length(Ttau))
  for(i in Ttau)
  {
    trend_1<-(1:i)^p
    trend_2<-((i+1):T)^p
    resi_1<-lm(x[1:i]~trend_1)$residuals
    resi_2<-lm(x[(i+1):T]~trend_2)$residuals
    tstat[q]<-(T-i)^(-2)*sum(cumsum(resi_2)^2)/(i^(-2)*sum(cumsum(resi_1)^2))
    tr1<-matrix(nrow=i,ncol=z);for(f in 1:z){tr1[,f]<-(1:i)^f}
    tr2<-matrix(nrow=T-i,ncol=z);for(f in 1:z){tr2[,f]<-(1:(T-i))^f}
    model1<-lm(x[1:i]~tr1)
    model2<-lm(x[(i+1):T]~tr2)
    WT1[q]<-wald_test(vcov(model1),coef(model1),Terms=(p+2):(z+1))$result$chi2[1]
    WT2[q]<-wald_test(vcov(model2),coef(model2),Terms=(p+2):(z+1))$result$chi2[1]
    q<-q+1
  }
  return(list(tstat=tstat,WT1=WT1/T,WT2=WT2/T))
}
