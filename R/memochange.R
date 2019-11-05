#' memochange: Testing for Structural Breaks under Long Memory and Testing for Changes in Persistence
#'
#' Test procedures and break point estimators for persistent processes that exhibit structural breaks in mean or in persistence.\cr \cr 
#' On the one hand, the package contains the most popular approaches for testing whether a time series exhibits a break in persistence from I(0) to I(1) or vice versa, such as those of Busetti and Taylor (2004) and Leybourne, Kim, and Taylor (2007).
#' The approach by Martins and Rodrigues (2014), which allows to detect changes from I(d1) to I(d2) with d1 and d2 being non-integers, is included as well.\cr
#' In case the tests reject the null of constant persistence, various breakpoint estimators are available to detect the point of the break as well as the order of integration in the two regimes.\cr \cr
#' On the other hand, the package contains the most popular approaches to test for a change in mean in a long-memory time series, which were recently reviewed by Wenger, Leschinski, and Sibbertsen (2018). 
#' These include memory robust versions of the CUSUM, sup-Wald, and Wilcoxon type tests. The tests either utilize consistent estimates of the long-run variance or a self normalization approach in their test statistics.
#'
#' For details see the readme and vignettes in the corresponding GitHub repository (https://github.com/KaiWenger/memochange).
#'
#' @docType package
#' @seealso \code{\link{BP_estim}}, \code{\link{CUSUM_simple}},
#' \code{\link{cusum_test}}, \code{\link{CUSUMfixed}}, \code{\link{CUSUMLM}},
#' \code{\link{fixbsupw}}, \code{\link{LBI_test}}, \code{\link{LKSN_test}}, \code{\link{MR_test}},
#' \code{\link{pb_sim}}, \code{\link{ratio_test}}, \code{\link{snsupwald}},
#' \code{\link{snwilcoxon}}, \code{\link{wilcoxonLM}}
#'
#' @author Kai Wenger <Kai.Wenger@gmx.de>, Janis Becker
#' @name memochange
#' @references
#' Andrews, D. W. K. (1993): Tests for Parameter Instability and Structural Change With Unknown Change Point. Econometrica, 61, pp. 821-856.
#'
#' Betken, A. (2016): Testing for change-points in long-range dependent time series by means of a self-normalized wilcoxon test. Journal of Time Series Analysis, 37, pp. 785-908.
#'
#' Busetti, F. and Taylor, R. (2004): Tests of stationarity against a change in persistence. Journal of Econometrics, 123, pp. 33-66.
#'
#' Dehling, H. and Rooch, A. and Taqqu, M. S. (2012): Non-Parametric Change-Point Tests for Long-Range Dependent Data. Scandinavian Journal of Statistics, 40, pp. 153-173.
#'
#' Harvey, D., Leybourne, S. and Taylor, R. (2006): Modified tests for a change in persistence. Journal of Econometrics, 134, pp. 441-469.
#'
#' Horvath, L. and Kokoszka, P. (1997): The effect of long-range dependence on change-point estimators. Journal of Statistical Planung and Inference, 64, pp. 57-81.
#'
#' Hualde, J. and Iacone, F. (2017): Fixed bandwidth asymptotics for the studentized mean of fractionally integrated processes. Economics Letters, 150, pp. 39-43.
#'
#' Iacone, F. and Leybourne, S. J. and Taylor, R. A. M. (2014): A fixed-b Test for a Break in Level at an unknown Time under Fractional Integration. Journal of Time Series Analysis, 35, pp. 40-54.
#'
#' Leybourne, S., Kim, T., Smith, V., and Newbold, P. (2003): Tests for a change in persistence against the null of difference-stationarity. Econometrics Journal, 6, pp. 291-311.
#'
#' Leybourne, S. and Taylor, R. (2004): On tests for changes in persistence. Economics letters, 84, pp. 107-115.
#'
#' Leybourne, S., Kim, T., and Taylor, R. (2007): Cusum of squares-based tests for a change in persistence. Journal of Time Series Analysis, 28, pp. 408-433.
#'
#' Martins, L.. and Rodrigues, P. (2014): Testing for persistence change in fractionally integrated models: An application to world inflation rates Cusum of squares-based tests for a change in persistence. Computational Statistics and Data Analysis, 76, pp. 502-522.
#'
#' Shao, X. (2011): A simple test of changes in mean in the possible presence of long-range dependence. Journal of Time Series Analysis, 32, pp. 598-606.
#'
#' Sibbertsen, P. and Kruse, R. (2009): Testing for a break in persistence under long-range dependencies. Journal of Time Series Analysis, 30, pp. 263-285.
#'
#' Wang, L. (2008): Change-in-mean problem for long memory time series models with applications. Journal of Statistical Computation and Simulation, 78:7, pp. 653-668.
#'
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): Change-in-mean tests in long-memory time series: a review of recent developments. AStA Advances in Statistical Analysis, 103:2, pp. 237-256.
#'
#' Wenger, K. and Leschinski, C. and Sibbertsen, P. (2018): A simple test on structural change in long-memory time series. Economics Letters, 136, pp. 90-94.
NULL

