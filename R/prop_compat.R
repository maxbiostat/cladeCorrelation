#' Proportion of clades which are pairwise compatible.
#'
#' @param n number of taxa.
#' @param trivial logical. If \code{TRUE}, trivial clades with 1 and
#'  n members are considered (default is \code{FALSE}).
#' @param log logical. \code{TRUE}, the log proportion is returned (default is \code{FALSE}).
#' 
#' @return (log) proportion of clades which are pairwise compatible.
#' @importFrom  matrixStats logSumExp
#' @export prop_compat
#'
#' @examples
#' prop_compat(5)
#' prop_compat(50)
#' prop_compat(500)
#' ns <- 4:100
#' plot(ns, prop_compat(ns), lwd = 3, xlab = "Number of taxa (n)",
#' ylab = expression(c(n)), main = "Proportion compatible", type = "l", cex.lab = 1.5)

prop_compat <- function(n, trivial = FALSE, log = FALSE){ 
  lA <- matrixStats::logSumExp(c(log_diff_exp(n*log(3), (n + 1) *log(2)), 0))
  lB <- 2*log_diff_exp(n*log(2), 0)
  lC <- 2*log_diff_exp(n*log(2), log(n+2))
  lcompat <- ifelse(trivial,
                     matrixStats::logSumExp(c(log(3) + lA, lB/2)),
                     log_diff_exp(matrixStats::logSumExp(c(log(3) + lA, lB/2, lC)), lB)
                    )
  lnumEntries <- ifelse(trivial, lB, lC)
  ans <- lcompat-lnumEntries
  if(!log) ans <- exp(ans)
  return(ans)
}
prop_compat <- Vectorize(prop_compat)