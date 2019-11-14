context("LBI_test")
set.seed(410)
x=rnorm(100)
expect_error(LBI_test(x,trend="limear"))
expect_error(LBI_test(x,tau=0))
expect_error(LBI_test(x,tau=2))
expect_error(LBI_test(x,m=1))
expect_warning(LBI_test(x,simu=0,tau=0.1))
x=stats::ts(x)
expect_error(LBI_test(x))
x=matrix(1:100,ncol=2)
expect_error(LBI_test(x))

test_that("",{
          skip_on_cran()
          x=rnorm(100)
          trend = c("none", "linear")
          statistic = c("mean", "max", "exp")
          for(a in 1:2){
            tr=trend[a]
            for(b in 1:3){
              st=statistic[b]
              q=0
              simu=LBI_test(x,trend=tr,statistic=st,simu=1,M=5000)
              true=LBI_test(x,trend=tr,statistic=st)
              test=true/simu
              test= test<0.9 | test>1.1
              expect_equal(sum(test),0) #expect that critical values given by Busetti and Taylor (2004) and simulated critical values should be close to identical
            }
          }
})