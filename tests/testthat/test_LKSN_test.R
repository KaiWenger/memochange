context("LKSN_test")
set.seed(410)
x=rnorm(100)
expect_error(LKSN_test(x,trend="limear"))
expect_error(LKSN_test(x,tau=0))
expect_error(LKSN_test(x,tau=2))
expect_error(LKSN_test(x,m=1))
expect_warning(LKSN_test(x,simu=0,tau=0.15))
x=stats::ts(x)
expect_error(LKSN_test(x))
x=matrix(1:100,ncol=2)
expect_error(LKSN_test(x))

test_that("",{
  skip_on_cran()
  #size
  trend = c("none", "linear")
  statistic = c("mean", "max", "exp")
  for(a in 1:2){
    tr=trend[a]
    q=0
    for(i in 1:5){
      x=cumsum(rnorm(75))
      mod=LKSN_test(x,trend=tr)
      expect_lt(mod[1,2],mod[1,1]) #critical values should be in increasing order
      expect_equal(ncol(mod),3)
      expect_equal(min(mod[1,3],mod[2,3]),mod[3,3]) #max of test statistics should equal test statistic in both directions 
      q=q+sum(mod[,3]<mod[,2])
    }
    expect_lt(q,9)  #test should not reject H0 (which is true) in more than 8 of 10 cases at the 95 percent level
  }
  
  #test has power
  for(a in 1:2){
    tr=trend[a]
    q=0
    for(i in 1:5){
      x=pb_sim(100, 0.5, "none", d1=0, d2=1, mean=0, var=1)
      mod=LKSN_test(x,trend=tr)
      q=q+sum(mod[,3]<mod[,1])
    }
    expect_gt(q,2) #test should reject H0 at least in three of 10 cases at the 90 percent level
  }
})
