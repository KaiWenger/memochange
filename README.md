
<!-- README.md is generated from README.Rmd. Please edit that file -->
memochange
==========

Testing for Structural Breaks under Long Memory and Testing for Changes in Persistence
--------------------------------------------------------------------------------------

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/memochange)](https://CRAN.R-project.org/package=memochange) [![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/memochange)](https://CRAN.R-project.org/package=memochange) [![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0) [![BuildStatus](https://api.travis-ci.org/KaiWenger/memochange.svg?branch=master)](https://travis-ci.com/KaiWenger/memochange)

The `memochange` package implements the most prominent test procedures and break point estimators for persistent processes that exhibit structural breaks in mean or in persistence. On the one hand, the package contains the most popular approaches for testing whether a time series exhibits a break in persistence from I(0) to I(1) or vice versa. In case the tests reject the null of constant persistence, various breakpoint estimators are available to detect the point of the break as well as the order of integration in the two regimes. On the other hand, the package contains the most popular approaches to test for a change-in-mean in a long-memory time series. These include memory robust versions of the CUSUM, sup-Wald, and Wilcoxon type tests. The tests either utilize consistent estimates of the long-run variance or a self normalization approach in their test statistics.

`memochange` may be of interest to those who analyse macroeconomic and financial time series such as inflation rates, trading volume, interest rates, volatilities and so on.

You can install this R package from GitHub:

``` r
install.packages("devtools")
library(devtools)
install_github("KaiWenger/memochange")
```

or directly from the CRAN repository:

``` r
install.packages("memochange")
```

Examples
--------

In this section we present two short examples that illustrate how the implemented procedures can be used on a real data. A more detailed presentation of the various tests and functions can be found in the vignette.

### Tests for persistence change

As an empirical example for a process that might exhibit a break in persistence, we consider the price of crude oil. First, we download the monthly price series from the FRED data base. For this purpose we need the `data.table` package.

``` r
oil=data.table::fread("https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=MCOILWTICO&scale=left&cosd=1986-01-01&coed=2019-08-01&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Monthly&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&transformation=lin&vintage_date=2019-09-23&revision_date=2019-09-23&nd=1986-01-01")
```

To get a first visual impression, we plot the series.

``` r
oil=as.data.frame(oil)
oil$DATE=zoo::as.Date(oil$DATE)
oil_xts=xts::xts(oil[,-1],order.by = oil$DATE)
zoo::plot.zoo(oil_xts,xlab="",ylab="Price",main="Crude Oil Price: West Texas Intermediate")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" style="display: block; margin: auto;" />

From the plot we observe that the series seems to be more variable in its second part from year 2000 onwards. This is first evidence that a change in persistence has occurred. We can test this hypothesis using the functions `cusum_test`, `LBI_test`, `LKSN_test`, `MR_test`, and `ratio_test`. In this short example we only consider the MR test as it is the most general one of the five implemented. The functionality of the other tests is similar. They all require a univariate numeric vector `x` as an input variable and yield a matrix of test statistic and critical values as an output variable.

``` r
library(memochange)
x <- as.numeric(oil[,2])
```

Applying the default version of the MR test by Martins and Rodrigues (2014) yields

``` r
MR_test(x)
#>                                           90%       95%       99%
#> Against increase in memory          -1.638729 -1.925646 -2.503650
#> Against decrease in memory          -1.653875 -1.953697 -2.509242
#> Against change in unknown direction -1.937915 -2.210129 -2.726981
#>                                     Teststatistic
#> Against increase in memory              -3.300664
#> Against decrease in memory              -1.629719
#> Against change in unknown direction     -3.300664
```

Here, test statistic and critical values for the null of constant persistence against an increase in persistence, a decrease in persistence, and a change in an unknown direction are displayed in a matrix. The latter accounts for the fact that we perform two tests facing a multiple testing problem. The results suggest that an increase in persistence has occurred somewhere in the series since the test statistic exceeds the critical value at the one percent level. In addition, this value is also significant when accounting for the multiple testing problem.

We can modify this default version by choosing the arguments `trend`, `tau`, `statistic`, `serial`, `simu`, and `M`. This amongst other things allows the test to be also applied for series who are suspected to exhibit linear trends or short run dynamics. Further details can be found in the vignette and on the help page of the MR test for details.

The test indicates that the oil price series exhibits an increase in memory over time. To correctly model and forecast the series, the exact location of the break is important. This can be estimated by the `BP_estim` function. It is important for the function that the direction of the change is correctly specified. In our case, an increase in memory has occurred so that we set direction="01".

``` r
BP_estim(x,direction="01")
#> $Breakpoint
#> [1] 151
#> 
#> $d_1
#> [1] 0.855238
#> 
#> $sd_1
#> [1] 0.129834
#> 
#> $d_2
#> [1] 1.034389
#> 
#> $sd_2
#> [1] 0.1468516
```

This yields a list stating the location of the break (observation 151), semiparametric estimates of the order of integration in the two regimes (0.86 and 1.03) as well as the standard deviations of these estimates (0.13 and 0.15).

``` r
oil$DATE[151]
#> [1] "1998-07-01"
```

Consequently, the function indicates that there is a break in persistence in July, 1998. This means that from the beginning of the sample until June 1998 the series is integrated with an order of 0.85 and from July 1998 on the order of integration increased to 1.03.

The function further allows for various types of break point and memory estimators. These are presented in the vignette.

### Tests for change-in-mean

As an example how to conduct the change-in-mean tests implemented in the `memochange` package, we compare the performance (in terms of size and power) of two implemented tests via a Monte Carlo simulation. We choose the two sample sizes n=\[100,500\] and memory parameters d=\[0.1,0.2\]. To estimate the memory parameter we apply the local Whittle estimator by Robinson (1995). The setup is very similar to the published paper by Wenger and Leschinski (2019) who introduce the two applied tests. The simulation can be extended by all other change-in-mean tests that are implemented in the `memochange` package, which is not been done here to save computing time. The whole following code should run just around a few seconds.

``` r
library(memochange)

test_func<-function(n,d)
{
  tseries     <- fracdiff::fracdiff.sim(n=n,d=d)$series
  
  changep     <- c(rep(0,n/2),rep(1,n/2))
  tseries2    <- tseries+changep

  d_est       <- LongMemoryTS::local.W(tseries, m=floor(1+n^0.65))$d
  d_est2      <- LongMemoryTS::local.W(tseries2, m=floor(1+n^0.65))$d
  
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

The code implements the function the Monte Carlo simulation is based on. First, a long-memory time series `tseries` of length n with memory d is simulated that is not subject to a change-in-mean. For the second time series `tseries2` a mean shift of size 1 is added. Then the memory parameters are estimated by local Whittle approach and the fixed-m CUSUM type A and B tests of Wenger and Leschinski (2019) are applied. Arguments that must be provided to the tests are the estimated long-memory parameters, the type of procedure to be applied, as well as the bandwidth that is used in fixed-m estimation. We choose the values suggested in Wenger and Leschinski (2019).

In the next step the Monte Carlo simulation (250 replications) is done. To do so, the `MonteCarlo` package needs to be loaded. The results of the Monte Carlo simulation is viewed as a LaTeX table.

``` r
n_grid  <- c(100,500)
d_grid  <- c(0.1,0.2)

param_list=list("n"=n_grid, "d"=d_grid)
erg     <- MonteCarlo::MonteCarlo(func=test_func, nrep=250, param_list=param_list, ncpus=1)

rows    <- c("n")
cols    <- c("d")
MonteCarlo::MakeTable(output=erg, rows=rows, cols=cols, digits=2, include_meta = FALSE)
```

Both tests should hold their nominal size (first two tables), which is chosen as 5%. Furthermore, the power of the tests (third and fourth table) should increase as the sample size n increases and decrease as the memory d increases. Example results of a run of the Monte Carlo code are given in the next table.

| n/d   |   0.1  |  0.2  |   0.1  |  0.2  |
|:------|:------:|:-----:|:------:|:-----:|
|       |        |       |        |       |
| size  | type A |       | type B |       |
|       |        |       |        |       |
| 100   |  0.052 | 0.040 |  0.052 | 0.056 |
| 500   |  0.068 | 0.052 |  0.068 | 0.056 |
|       |        |       |        |       |
| power | type A |       | type B |       |
|       |        |       |        |       |
| 100   |  0.688 | 0.440 |  0.640 | 0.416 |
| 500   |  1.000 | 0.900 |  1.000 | 0.848 |

Contributions
-------------

We welcome any and all contributions! To make the process as painless as possible for all involved, please see our [guide to contributing](contributing.md)

References
----------

Martins, Luis F, and Paulo MM Rodrigues. 2014. “Testing for Persistence Change in Fractionally Integrated Models: An Application to World Inflation Rates.” *Computational Statistics & Data Analysis* 76. Elsevier: 502–22. doi:[10.1016/j.csda.2012.07.021](https://doi.org/10.1016/j.csda.2012.07.021).
