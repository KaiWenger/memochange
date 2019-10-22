---
title: "memochange: An R package for estimation procedures and tests for persistent time series"
tags:
  - R
  - time series analysis
  - long memory
  - structural change
  - breaks in persistence

authors:
  - name: Kai Wenger
    orcid: 0000-0002-6566-0427
    affiliation: 1 
  - name: Janis Becker
    orcid: 0000-0001-5104-6330
    affiliation: 1

affiliations:
 - name: Institute of Statistics, Leibniz University Hannover
   index: 1

date: 01 September 2019
bibliography: paper.bib
biblio-style: "apalike"
output: html_document
---

# Summary
The degree of memory is an important determinant of the characteristics of a time series. For an $I(0)$, or short-memory, process (e.g., AR(1) or ARMA(1,1)), the impact of shocks is short-lived and dies out quickly. On the other hand, for an $I(1)$, or difference-stationary, process such as the random walk, shocks persist infinitely. Thus, any change in a variable will have an impact on all future realizations. For an $I(d)$, or long-memory, process with $0<d<1$, shocks neither die out quickly nor persist infinitely, but have a hyperbolically decaying impact. In this case, the current value of a variable depends on past shocks, but the less so the further these shocks are past.

There are plenty of procedures to determine the memory of a series (see @robinson1995gaussian, @shimotsu2010exact, among others). However, there is also the possibility that series exhibit a structural change in memory, often referred to as a change in persistence. Starting with @kim2000detection various procedures have been proposed to detect these changes and consistently estimate the change point. @busetti2004tests and @leybourne2004tests suggest approaches for testing the null of constant I(0) behaviour of the series against the alternative that a change from either $I(0)$ to $I(1)$ or $I(1)$ to $I(0)$ occurred. However, both approaches show serious distortions if neither the null nor the alternative is true, i.e. the series is constant $I(1)$. In this case the procedures by @leybourne2003tests and @leybourne2007cusum can be applied as they have the same alternative, but assume constant $I(1)$ behaviour under the null. Again, the procedures exhibit distortions when neither the null nor the alternative is true. To remedy this issue, @harvey2006modified suggest an approach that entails the same critical values for constant $I(0)$ and constant $I(1)$ behavior. Consequently, it accommodates both, constant $I(0)$ and constant $I(1)$ behavior under the null. 

While this earlier work focussed on the $I(0)$/$I(1)$ framework, more recent approaches are able to detect changes from $I(d_1)$ to $I(d_2)$ where $d_1$ and $d_2$ are allowed to be non-integers. @sibbertsen2009testing extend the approach of @leybourne2007cusum such that the testing procedure consistently detects changes from $0 \leq d_1<1/2$ to $1/2<d_2<3/2$ and vice versa. Under the null the test assumes constant $I(d)$ behavior with $0 \leq d <3/2$. The approach suggested by @martins2014testing is even able to identify changes from $-1/2<d_1<2$ to $-1/2<d_2<2$ with $d_1 \neq d_2$. Here, under the null the test assumes constant $I(d)$ behavior with $-1/2<d<2$. 


Structural changes cannot only occur in the persistence, but also in the mean of a persistent time series.  A large body of literature handles this topic for weakly dependent $I(0)$ time series. However, these tests for a mean shift are invalidated for $d>0$. The problem is that the limiting distributions of the standard tests are different under long memory. Therefore, they reject the null hypothesis of a constant mean with probability of one asymptotically.

In recent years, a lot of researchers modified standard testing procedures for a change-in-mean to account for $0 \leq d <0.5$. A review is given in @wenger2019change. The tests can be divided in two dimensions: First, whether they are adapted on the CUSUM, sup-Wald, or Wilcoxon test. Second, which type of long-run variance estimator they apply in their test statistics. Early approaches utilize consistent estimates of the long-run variance (@horvath1997effect, @wang2008change, @dehling2013non), while more recent contributions suggest self-normalized test statistics (@shao2011simple, @iacone2014fixed, @betken2016testing, @wenger2019fixed). This may be the reason since it is found in @wenger2019change that self-normalized tests are robust against size distortions.

# Statement of Need

memochange is an R (@RCoreTeam) package that contains implementations of the most prominent tests for changes in persistence and for change-in-mean in persistent time series. Stationary $I(0)$ and non-stationary $I(d)$, $d>0$ time series have very different properties. Therefore, our package can be used to check whether there is a break in persistence and/or a change-in-mean. This can avoid model misspecification and improve forecasting models. In macroeconomics and finance, there is a wide range of possible time series where the procedures can be applied to, e.g. inflation rates, volatilities, and trading volume (cf. @fleming2011long, @Xu14, @hassler2014detecting)

Recent scholarly publications that apply the approaches that the memochange package implements are, for example, @sibbertsen2014testing, @baillie2019long, and @wenger2019fixed.

The memochange package can be downloaded from CRAN or installed directly from the source code at the corresponding GitHub repository (which is the latest version of the package).

# Acknowledgements

Financial support of the Deutsche Forschungsgemeinschaft (DFG) is gratefully acknowledged. We would like to thank Christian Leschinski and Simon Wingert for their helpful comments and suggestions.

# References