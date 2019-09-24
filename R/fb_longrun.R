#' function to estimate the fixed-b long-run variance. For internal use only.
#' @keywords internal
#'
fb_longrun <- function(tseries,m,type="Bartlett")
{
n               <- length(tseries)
cova            <- apply(matrix(1:(n-1), ncol=1), 1, function(i) return(sum(tseries[1:(n-i)]*tseries[(i+1):n])))/n
kern            <- as.matrix(1 - seq(1,(n-1),1)/(m+1))
kern[kern<0]    <- 0
sigm            <- (sum(tseries**2)/n)+2*sum(kern*cova)
return(sigm)
}


