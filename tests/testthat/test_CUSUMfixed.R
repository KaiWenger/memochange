context("CUSUMfixed")
set.seed(410)
T     <- 100
d     <- 0
x     <- fracdiff::fracdiff.sim(n=T, d=d)$series

expect_error(CUSUMfixed(x, d=d, procedure="something_else", bandw=0.1))
expect_error(CUSUMfixed(x, d=d, procedure="CUSUMfixedb_typeA", bandw=10))
expect_error(CUSUMfixed(x, d=d, procedure="CUSUMfixedm_typeA", bandw=0.1))
expect_error(CUSUMfixed(c(x,NA), d=d, procedure="CUSUMfixedb_typeA", bandw=0.1))
expect_error(CUSUMfixed(x, d=d, procedure="CUSUMfixedb_typeA", bandw=0.1,tau=2))
expect_warning(CUSUMfixed(x, d=d, procedure="CUSUMfixedb_typeA", bandw=0.1,tau=0.1))
x=stats::ts(x)
expect_error(CUSUMfixed(x, d=d, procedure="CUSUMfixedb_typeA", bandw=0.1))


# size
T         = 100
procedure = c("CUSUMfixedb_typeA", "CUSUMfixedb_typeB", "CUSUMfixedm_typeA", "CUSUMfixedm_typeB")
d_grid    = c(0.1,0.2)
for(a in 1:length(procedure)){
  proc     = procedure[a]
  if((proc == "CUSUMfixedb_typeA" | proc == "CUSUMfixedb_typeB")) bandw=0.1
  if((proc == "CUSUMfixedm_typeA" | proc == "CUSUMfixedm_typeB")) bandw=10
  for(b in 1:length(d_grid)){
    d     = d_grid[b]
    q     = 0
      for(i in 1:15){
        x     = fracdiff::fracdiff.sim(n=T, d=d)$series
        mod   = CUSUMfixed(x, d=d, procedure=proc, bandw=bandw)
        q     = q+sum(mod[4]>mod[3])
      }
      expect_lt(q,11)  #test should not reject H0 (which is true) in more than 10 of 15 cases at the 99 percent level
      
  }
}

# power
T         = 100
procedure = c("CUSUMfixedb_typeA", "CUSUMfixedb_typeB", "CUSUMfixedm_typeA", "CUSUMfixedm_typeB")
d_grid    = c(0.1,0.2)
for(a in 1:length(procedure)){
  proc     = procedure[a]
  if((proc == "CUSUMfixedb_typeA" | proc == "CUSUMfixedb_typeB")) bandw=0.1
  if((proc == "CUSUMfixedm_typeA" | proc == "CUSUMfixedm_typeB")) bandw=10
  for(b in 1:length(d_grid)){
    d     = d_grid[b]
    q     = 0
    for(i in 1:15)
    {
      x       = fracdiff::fracdiff.sim(n=T, d=d)$series
      changep = c(rep(0,T/2), rep(1,T/2))
      x       = x+changep
      mod     = CUSUMfixed(x, d=d, procedure=proc, bandw=bandw)
      q       = q+sum(mod[4]>mod[1])
    }
      expect_gt(q,2) #test should reject H0 at least in three of 15 cases at the 90 percent level
  }
}