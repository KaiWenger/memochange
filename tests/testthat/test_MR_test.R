context("MR_test")
set.seed(410)
x=rnorm(100)
expect_error(MR_test(x,trend="limear"))
expect_error(MR_test(x,tau=0))
expect_error(MR_test(x,tau=2))
expect_error(MR_test(x,type="SK"))
expect_error(MR_test(x,m=100))
expect_warning(MR_test(x,simu=0,tau=0.1))
x=stats::ts(x)
expect_error(MR_test(x))
x=matrix(1:100,ncol=2)
expect_error(MR_test(x))

#size
test_that("",{
  skip_on_cran()
  trend = c("none", "linear")
  statistic = c("squared", "standard")
  serial = c(FALSE,TRUE)
  for(a in 1:2){
    tr=trend[a]
    for(b in 1:2){
      st=statistic[b]
      for(c in 1:2){
        q=0
        se=serial[c]
        for(i in 1:5){
          x=rnorm(75)
          mod=MR_test(x,trend=tr,statistic=st,serial=se)
          if(b==1){
            expect_lt(mod[1,1],mod[1,2]) #critical values should be in ascending order
            expect_equal(max(mod[1,4],mod[2,4]),mod[3,4]) #max of test statistics should equal test statistic in both directions 
            q=q+sum(mod[,4]>mod[,3])
          }
          if(b==2){
            expect_lt(mod[1,2],mod[1,1]) #critical values should be in ascending order
            expect_equal(min(mod[1,4],mod[2,4]),mod[3,4]) #max of test statistics should equal test statistic in both directions 
            q=q+sum(mod[,4]<mod[,3])
          }
        }
        expect_lt(q,11)  #test should not reject H0 (which is true) in more than 10 of 15 cases at the 99 percent level
      }
    }
  }
  
  #test has power
  for(a in 1:2){
    tr=trend[a]
    for(b in 1:2){
      st=statistic[b]
      for(c in 1:2){
        q=0
        se=serial[c]
        for(i in 1:5){
          x=pb_sim(100, 0.5, "none", d1=0.1, d2=0.9, mean=0, var=1)
          mod=MR_test(x,trend=tr,statistic=st,serial=se)
          if(b==1) q=q+sum(mod[,4]>mod[,1])
          if(b==2) q=q+sum(mod[,4]<mod[,1])
        }
        #expect_gt(q,2) #test should reject H0 at least in three of 15 cases at the 90 percent level
      }
    }
  }
})

