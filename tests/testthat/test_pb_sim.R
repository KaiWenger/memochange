context("pb_sim")

x=pb_sim(100,0.5,"linear",tp=1,d1=0,d2=0,mean=0,var=1)
index=1:length(x)
expect_gt(stats::coef(stats::lm(x~index))[2],0.5)

x=pb_sim(100,0.5,"none",tp=0,d1=0,d2=0,mean=2,var=1)
expect_gt(stats::coef(stats::lm(x~1))[1],0)

x=pb_sim(200,0.4,"none",tp=0,d1=0,d2=1,mean=0,var=1)
expect_equal(stats::na.omit(length(x)),200)
expect_lt(max(x),Inf)
expect_gt(min(x),-Inf)

x2=pb_sim(200,0.4,"none",tp=0,d1=0,d2=1,mean=0,var=10)
expect_gt(stats::var(x2),stats::var(x))

expect_true(is.vector(pb_sim(100,0.5,"none",d1=0.2,d2=0.8,mean=0,var=1)))
expect_error(pb_sim(100,0.5,"limear",tp=1,d1=0,d2=0,mean=0,var=1))

expect_warning(pb_sim(100,0.5,"linear",tp=0,d1=0,d2=0,mean=2,var=1))