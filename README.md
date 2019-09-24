
# memochange
## Testing for Structural Breaks under Long Memory and Testing for Changes in Persistence
<!--- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/memochange)](https://CRAN.R-project.org/package=memochange) -->


The `memochange` package implements the most prominent test procedures and break point estimators for persistent processes that exhibit structural breaks in mean or in persistence.
On the one hand, the package contains the most popular approaches for testing whether a time series exhibits a break in persistence from I(0) to I(1) or vice versa. In case the tests reject the null of constant persistence, various breakpoint estimators are available to detect the point of the break as well as the order of integration in the two regimes.
On the other hand, the package contains the most popular approaches to test for a change-in-mean in a long-memory time series. 
These include memory robust versions of the CUSUM, sup-Wald, and Wilcoxon type tests. The tests either utilize consistent estimates of the long-run variance or a self normalization approach in their test statistics.

`memochange` may be of interest to those who analyse macroeconomic and financial time series such as inflation rates, trading volume, interest rates, volatilities and so ons


You can install this R package from GitHub:


```r
install.packages("devtools")
library(devtools)
install_github("KaiWenger/memochange")
```

or directly from the CRAN repository:


```r
install.packages("memochange")
```


## Examples

In this section we will show how the implemented procedures can be used on a real data set as well as in Monte Carlo simulations to compare the performance of the tests.

### Tests for persistence change

As an empirical example for a process that might exhibit a break in persistence, we consider the price of crude oil.  First, we download the monthly price series from the FRED data base. For this purpose we need the `data.table` package.

```r
library(data.table)
oil=fread("https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=
line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&
txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=
yes&show_tooltip=yes&id=MCOILWTICO&scale=left&cosd=1986-01-01&coed=2019-08-01&line_color=
%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=
0&fml=a&fq=Monthly&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&transformation=lin&vintage_date=
2019-09-23&revision_date=2019-09-23&nd=1986-01-01")
```

To get a first visual impression, we plot the series.

```r
library(zoo)
library(xts)
oil=as.data.frame(oil)
oil$DATE=as.Date(oil$DATE)
oil_xts=xts(oil[,-1],order.by = oil$DATE)
plot.zoo(oil_xts,xlab="",ylab="Price",main="Crude Oil Price: West Texas Intermediate")
```
![](oil_plot.eps)

From the plot we observe that the series seems to be more variable in its second part from year 2000 onwards.  This is first evidence that a change in persistence has occurred. We can test this hypothesis using the functions `cusum_test`, `LBI_test`, `LKSN_test`, `MR_test`, and `ratio_test`. In this short example we use the ratio and MR test since these are the empirically most often applied one. The functionality of the other tests is similar. They all require a univariate numeric time series `x` as an input variable and yield a matrix of test statistic and critical values as an output variable.

```r
library(memochange)
x <- as.numeric(oil[,2])
```

As a starting point the default version of the ratio test is applied.

```r
ratio_test(x)

#>                                        90%    95%    99% Teststatistic
#> Against change from I(0) to I(1)    3.5148 4.6096 7.5536    225.943543
#> Against change from I(1) to I(0)    3.5588 4.6144 7.5304      1.170217
#> Against change in unknown direction 4.6144 5.7948 9.0840    225.943543

```

This yields a matrix that gives test statistic and critical values for the null of constant $I(0)$ against a change from $I(0)$ to $I(1)$ or vice versa. Furthermore, the statistics for a change in an unknown direction are included as well. This accounts for the fact that we perform two tests facing a multiple testing problem. The results suggest that a change from $I(0)$ to $I(1)$ has occurred somewhere in the series since the test statistic exceeds the critical value at the one percent level. In addition, this value is also significant when accounting for the multiple testing problem. Consequently, the default version of the ratio test suggests a break in persistence. 

We can modify this default version by choosing the arguments `trend`, `tau`, `statistic`, `type`, `m`, `z`, `simu`, and `M` (see the help page of the ratio test for details). 
The plot does not indicate a linear trend so that it seems unreasonable to change the trend argument. Also, the plot suggests that the break is rather in the middle of the series than at the beginning or the end so that changing the break fraction seems unnecessary as well. 
The type of test statistic calculated can be easily changed using the statistic argument. However, simulation results indicate mean, max, and exp statistics to deliver qualitatively similar results.

Something that is of more importance is the type of test performed. The default version considers the approach by Busetti and Taylor (2004). In case of a constant $I(1)$ process this test often spuriously identifies a break in persistence. Harvey, Leybourne and Taylor (2006) account for this issue by adjusting the test statistic such that its critical values are the same under constant $I(0)$ and constant $I(1)$.  We can calculate their test statistic by setting `type="HLT"`. For this purpose, we need to state the number of polynomials used in their test statistic. The default value is $9$ as suggested by Harvey, Leybourne and Taylor (2006).
Choosing another value is only sensible for very large data sets (number of obs. > $10000$) where the test statistic cannot be calculated due to computational singularity. 
In this case decreasing $z$ can allow the test statistic to be calculated. This, however, invalidates the critical values. 
Since our data set is rather small we can stick with the default of $z=9$.

```r
ratio_test(x,type="HLT")

#>                                        90%    95%    99% Teststatistic 90% Teststatistic 95% Teststatistic 99%
#> Against change from I(0) to I(1)    3.5148 4.6096 7.5536        58.9078204        43.4772689        25.3369507
#> Against change from I(1) to I(0)    3.5588 4.6144 7.5304         0.3085495         0.2290113         0.1290305
#> Against change in unknown direction 4.6144 5.7948 9.0840        44.2171379        34.1367566        20.0058559

```

Again the test results suggests that there is a break from $I(0)$ to $I(1)$. Consequently, it is not a constant $I(1)$ process that led to a spurious rejection of the test by Busetti and Taylor (2004).

Another test for a change in persistence is that by Martins and Rodrigues (2014). This is more general as it is not restricted to the $I(0)$/$I(1)$ framework, but can identify changes from $I(d_1)$ to $I(d_2)$ with $d_1 \neq d_2$ and $-1/2<d_1,d_2<2$.
The default version is applied by

```r
MR_test(x)

#>                                           99%       95%       90% Teststatistic
#> Against increase in memory          -2.503650 -1.925646 -1.638729     -3.374027
#> Against decrease in memory          -2.509242 -1.953697 -1.653875     -1.595476
#> Against change in unknown direction -2.726981 -2.210129 -1.937915     -3.374027

```

Again, the function returns a matrix consisting of test statistic an critical values. Here, the alternative of the test is an increase respectively a decrease in memory and the null is rejected in case the test statistic is lower than the critical values. Consequently, in line with the results of the ratio test, the approach by Martins and Rodrigues (2014) suggest that the series exhibits an increase in memory, i.e. that the memory of the series increases from $d_1$ to $d_2$ with $d_1<d_2$ at some point in time. Again, this also holds if we consider the critical values that account for the multiple testing problem.

Similar to the ratio test and all other tests against a change in persistence in the `memochange` package, the MR test also has the same arguments `trend`, `tau`, `simu`, and `M`. Furthermore, we can choose again the type of test statistic. This time we can decide whether to use the standard statistic or the squared statistic.

```r
MR_test(x,statistic="squared")

#>                                          90%      95%      99% Teststatistic
#> Against increase in memory          4.291558 5.426332 8.216863     12.114086
#> Against decrease in memory          4.012380 5.021871 7.531425      2.545542
#> Against change in unknown direction 5.062131 6.206093 8.976911     12.114086

```

As for the ratio test, changing the type of statistic has a rather small effect on the empirical performance of the test. The same is also true for the twostep argument which makes the function calculate the twostep version of the test if set to `TRUE`. Lastly, if we believe that the underlying process exhibits additional short run components, we can account for these by setting $p$ larger zero such as

```r
MR_test(x,statistic="squared",p=1)

#>                                          90%      95%      99% Teststatistic
#> Against increase in memory          4.129037 5.215535 7.924110      52.60968
#> Against decrease in memory          4.259275 5.388967 8.142119      19.69812
#> Against change in unknown direction 5.129512 6.280378 9.073346      52.60968

```

While the critical values and test statistic change, the conclusion remains the same.

All tests indicate that the oil price series exhibits an increase in memory over time. To correctly model and forecast the series, the exact location of the break is important.
This can be estimated by the `BP_estim` function. It is important for the function that the direction of the change is correctly specified. In our case, an increase in memory has occurred so that we set direction="01"

```r
BP_estim(x,direction="01")
```

This yields the location of the break (observation $151$), semiparametric estimates of the order of integration in the two regimes ($0.85$ and $1.03$) as well as the standard deviations of these estimates ($0.13$ and $0.15$).

```r
oil$DATE[151]

#> [1] "1998-07-01"
```

Consequently, the function indicates that there is a break in persistence in July, 1998. 
This means that from the beginning of the sample until June 1998 the series is integrated with an order of 0.85 and from July 1998 on the order of integration increased to 1.03.

As before, the function allows for various types of break point estimators. Instead of the default estimator of Busetti and Taylor (2004), one can also rely on the estimator of Leybourne, Kim, and Taylor (2006) by setting `type="LKT"`.
This estimator relies on estimates of the long-run variance. Therefore, it is also needed that $m$ is chosen, which determines how many covariances are used when estimating the long-run variance. Leybourne, Kim, and Taylor (2006) suggest $m=0$.

```r
BP_estim(x,direction="01",type="LKT",m=0)
```

This yields a similar result with the break point lying in the year 1998 and $d$ increasing from approximately $0.8$ to approximately $1$.

All other arguments of the function (`trend`, `tau`, `p`, `twostep`) were already discussed above except for `d_estim` and `d_bw`. These two arguments determine which estimator and bandwidth are used to estimate the order of integration in the two regimes. Concerning the estimator, the GPH (Geweke and Porter-Hudak (1983)) and the exact local Whittle estimator (Shimotsu and Phillips (2005)) can be selected. Although the exact local Whittle estimator has a lower variance, the GPH estimator is still often considered in empirical applications due to its simplicity. In our example the results of the two estimators are almost identical.

```r
BP_estim(x,direction="01",d_estim="GPH")
```

The `d_bw` argument determines how many frequencies are used for estimation. Larger values of $m$ imply a lower variance of the estimates, but also bias the estimator if the underlying process possesses short run dynamics.
Usually a value between $0.5$ and $0.8$ is considered.

```r
BP_estim(x,direction="01",d_bw=0.75)
BP_estim(x,direction="01",d_bw=0.65)
```

In our setup, it can be seen that increasing `d_bw` to $0.75$ does not severely change the estimated order of integration in the two regimes. Decreasing `d_bw`, however, leads to smaller estimates of $d$.\ 

### Tests for change-in-mean

As an example how to conduct the change-in-mean tests implemented in the `memochange` package,  we compare the performance (in terms of size and power) of two implemented tests via a Monte Carlo simulation. We choose the two sample sizes $n=[100,500]$ and memory parameters $d=[0.1,0.2]$. To estimate the memory parameter we apply the local Whittle estimator by Robinson (1995). The setup is very similar to the published paper by Wenger and Leschinski (2019) who introduce the two applied tests. The simulation can be extended by all other change-in-mean tests that are implemented in the `memochange` package, which is not been done here to save computing time. The whole following code should run just around a few seconds.

```r
library(memochange)
library(fracdiff)
library(longmemo)
library(LongMemoryTS)

test_func<-function(n,d)
{
  tseries     <- fracdiff.sim(n=n,d=d)$series
  
  changep     <- c(rep(0,n/2),rep(1,n/2))
  tseries2    <- tseries+changep

  d_est       <- local.W(tseries, m=floor(1+n^0.65))$d
  d_est2      <- local.W(tseries2, m=floor(1+n^0.65))$d
  
  typeAsize   <- CUSUMfixed(tseries,d=d_est,procedure="CUSUMfixedm_typeA",bandw=10)
  typeBsize   <- CUSUMfixed(tseries,d=d_est,procedure="CUSUMfixedm_typeB",bandw=10)
  
  typeApower  <- CUSUMfixed(tseries2,d=d_est2,procedure="CUSUMfixedm_typeA",bandw=10)
  typeBpower  <- CUSUMfixed(tseries2,d=d_est2,procedure="CUSUMfixedm_typeB",bandw=10)
  
  decAsize    <- typeAsize["Teststatistic"] > typeAsize["95%"]
  decBsize    <- typeBsize["Teststatistic"] > typeBsize["95%"]
  
  decApower   <- typeApower["Teststatistic"] > typeApower["95%"]
  decBpower   <- typeBpower["Teststatistic"] > typeBpower["95%"]
  
  return(list("decAsize"=unname(decAsize),"decBsize"=unname(decBsize),"decApower"=unname(decApower),"decBpower"=unname(decBpower)))
}
```

The code implements the function the Monte Carlo simulation is based on. First, a long-memory time series `tseries` of length $n$ with memory $d$ is simulated that is not subject to a change-in-mean. For the second time series `tseries2` a mean shift of size $1$ is added. Then the memory parameters are estimated by local Whittle approach and the fixed-$m$ CUSUM type A and B tests of Wenger and Leschinski (2019) are applied. Arguments that must be provided to the tests are the estimated long-memory parameters, the type of procedure to be applied, as well as the bandwidth that is used in fixed-$m$ estimation. We choose the values suggested in Wenger and Leschinski (2019).

In the next step the Monte Carlo simulation ($250$ replications) is done. To do so, the `MonteCarlo` package needs to be loaded. The results of the Monte Carlo simulation is viewed as a LaTeX table.

```r
library(MonteCarlo)

n_grid  <- c(100,500)
d_grid  <- c(0.1,0.2)

param_list=list("n"=n_grid, "d"=d_grid)
erg     <- MonteCarlo(func=test_func, nrep=250, param_list=param_list, ncpus=1)

rows    <- c("n")
cols    <- c("d")
MakeTable(output=erg, rows=rows, cols=cols, digits=2, include_meta = FALSE)
```

Both tests should hold their nominal size (first two tables), which is chosen as $5\%$. Furthermore, the power of the tests (third and fourth table) should increase as the sample size $n$ increases and decrease as the memory $d$ increases. Example results of a run of the Monte Carlo code are given in the next table.

| n/d  | 0.1    |  0.2  |  0.1 |  0.2    |
| :-   | :----: | :---: | :---:| :----:  | 
|      |        |       |      |         |
| size | type A |       |type B|         | 
|      |        |       |      |         |
| 100  |  0.052 |  0.040|0.052 | 0.056   | 
| 500  |  0.068 |  0.052|0.068 | 0.056   | 
|      |        |       |      |         |
|power | type A |       |type B|         | 
|      |        |       |      |         | 
| 100  | 0.688  | 0.440 | 0.640| 0.416   | 
| 500  | 1.000  | 0.900 | 1.000| 0.848   |  

## Contributing

You are welcome to help us improving our package. Please submit an issue if you encounter any bugs, have a feature request, or have  any other requests. You can also contact the corresponding author via email Kai.Wenger@gmx.de.

## References

Wenger, K. and Leschinski, C. (2019). Fixed-bandwidth CUSUM tests under long memory. Econometrics and Statistics, forthcoming.  <https://doi.org/10.1016/j.ecosta.2019.08.001>.
