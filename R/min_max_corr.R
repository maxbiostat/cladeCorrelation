#' Minimum possible correlation between clades on \code{n} taxa.
#'
#' @param n number of taxa.
#'
#' @return mimimum value of the correlation matrix on \code{n} taxa.
#' @export min_corr
#'
#' @examples
#' ns <- round(seq(4, 1E3, length.out = 200))
#' plot(ns, min_corr(ns),  lwd = 3,
#' xlab = expression(n), ylab = expression(rho[min](n)),
#'  main = "Minimum correlation between clade indicators", 
#'  type = "l" , cex.lab = 1.5)
#'  abline(h = 0, lwd = 2, lty = 2)
min_corr <- function(n){
  return(
    -2/(3*n-5)
  )
}
min_corr <- Vectorize(min_corr)

#' Maximum possible correlation between clades on \code{n} taxa.
#'
#' @param n number of taxa.
#'
#' @return maximum value of the correlation matrix on \code{n} taxa.
#' @export max_corr
#'
#' @examples
#' ns <- round(seq(4, 1E3, length.out = 200))
#' plot(ns, max_corr(ns),  lwd = 3,
#' xlab = expression(n), ylab = expression(rho[max](n)),
#'  main = "Maximum correlation between clade indicators", 
#'  type = "l" , cex.lab = 1.5)
#'  abline(h = 0, lwd = 2, lty = 2)
max_corr <- function(n){
  kk <- ceiling(n/2)
  lpA <- pclade(kk, n, log = TRUE)
  lpB <- pclade(n-kk, n, log = TRUE)
  lpp <-  log(2)-log(n-1) - lchoose(n = n, k = kk) 
  log_rho <- log_diff_exp(lpp, lpA + lpB) - .5 * (lpA + log_diff_exp(0, lpA) + lpB + log_diff_exp(0, lpB))
  return(
    exp(log_rho)
  )
}
max_corr <- Vectorize(max_corr)