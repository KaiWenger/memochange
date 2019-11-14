context("cusum_test")
set.seed(410)
x=rnorm(100)
expect_equal(sum(cusum_test(x)<Inf),8)
expect_equal(sum(cusum_test(x)>-Inf),8)
expect_error(cusum_test(x,trend="limear"))
expect_error(cusum_test(x,tau=0))
expect_error(cusum_test(x,tau=2))
expect_error(cusum_test(x,type="BT"))
expect_error(cusum_test(x,m=100))
expect_warning(cusum_test(x,simu=0,tau=0.1))
x=stats::ts(x)
expect_error(cusum_test(x))
x=matrix(1:100,ncol=2)
expect_error(cusum_test(x))
test_that("",{
  skip_on_cran()
  #size
  trend = c("none", "linear")
  type = c("LKT", "SK")
  mm=c(0,1)
  for(a in 1:2){
    tr=trend[a]
    for(b in 1:2){
      ty=type[b]
      for(c in 1:2){
        m=mm[c]
        q=0
        for(i in 1:5){
          x=cumsum(rnorm(75))
          mod=cusum_test(x,trend=tr,type=ty,m=m)
          expect_gt(mod[1,1],mod[1,2]) #lower critical values should be in ascending order
          expect_lt(mod[2,1],mod[2,2]) #upper critical values in ascending order
          expect_equal(ncol(mod),4)
          expect_equal(mod[1,4],mod[2,4]) #test statistics should be equal 
          q=q+sum(mod[1,4]<mod[1,1])+sum(mod[2,4]>mod[2,1])
          }
        expect_lt(q,4) #test should not reject H0 (which is true) in more than 3 of 5 cases at the 99 percent level
      }
    } 
  }
  
  
  #test has power
  for(a in 1:2){
    tr=trend[a]
    for(b in 1:2){
      ty=type[b]
      q=0
      for(c in 1:2){
        m=mm[c]
        for(i in 1:5){
          x=pb_sim(100, 0.5, "none", d1=0, d2=1, mean=0, var=1)
          mod=cusum_test(x,trend=tr,type=ty,m=m)
          q=q+sum(mod[1,4]<mod[1,1])+sum(mod[2,4]>mod[2,1])
        }
        expect_gt(q,0)  #test should reject H0 in at least 1 of 5 cases at the 90 percent level
      }
    }
  }
})
