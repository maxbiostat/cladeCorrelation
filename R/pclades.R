
#' Auxiliary function for computing joint clade probability.
#'
#' @param a size of first clade.
#' @param b size of second clade.
#' @param n number of taxa.
#' @param log logical. If \code{TRUE}, returns the log (default is \code{FALSE}).
#'
#' @return (log) joint probability for cases 2 and 3.
#' @export Rn
Rn <- function(a, b, n, log = FALSE){
  ## This implementation automatically takes care of both cases: A \in B or B \in A
  mm <- min(a, b)
  MM <- max(a, b)
  ans <- ( log(4) + log(n) - (log(mm) + log1p(mm) + log1p(MM)) ) - lchoose(n = n, k = MM) - lchoose(n = MM, k = mm) 
  if(!log) ans <- exp(ans)  
  return(ans)
}

#' Auxiliary function for computing joint probabilities.
#'
#' @param a size of clade A.
#' @param b size of clade B.
#' @param n number of taxa.
#' @param log logical. If \code{TRUE}, returns the log joint probability (default is \code{FALSE}).
#'
#' @return (log) joint probability
#' @importFrom matrixStats logSumExp
#' @export Gn
#'
Gn <- function(a, b, n, log = FALSE){
  l1 <- log(n) - (log(a) + log(b) + log1p(a) + log1p(b))
  l2 <- matrixStats::logSumExp(c(
    log(a) + log1p(a),
    log(b) + log1p(b),
    log(a) + log(b)
  )) - {log(a) + log(b) + log1p(a) + log1p(b) + log1p(a + b)}
  l3 <- -{log(a+b) + log_diff_exp(2*log(a+b), 0)}
  ans <- matrixStats::logSumExp(c(log_diff_exp(l1, l2), l3))
  if(!log) ans <- exp(ans)
  return(ans)
}

#' Auxiliary function for computing joint probabilities 
#'
#' @param a size of clade A.
#' @param b size of clade B.
#' @param n number of taxa.
#' @param log logical. If \code{TRUE}, returns the log joint probability (default is \code{FALSE}).
#'
#' @return (log) joint probability.
#' @export rrn
#'
rrn <- function(a, b, n, log = FALSE){
  lK <- {log(4) + lfactorial(a) + lfactorial(b) + lfactorial(n - a - b)}-lfactorial(n-1)
  ans <- lK + Gn(a, b, n, log = TRUE)
  if(!log) ans <- exp(ans)
  return(ans)
}
#' Joint probability of two clades
#'
#' @param a  size of clade A.
#' @param b  size of clade B.
#' @param n number of taxa.
#' @param case see details
#' @param log logical if TRUE (default is FALSE) returns log probability
#'
#' @return (log) joint probability of clades A and B given their sizes.
#' @export pclades
#'
#' @examples
#' pclades(a = 2,  b = 3, n = 5,  case = 2)
#' @details case 1: A == B; 
#' 
#' case 2: A is contained in B;
#' 
#' case 3: B is contained in A;
#' 
#' case 4: A and B are disjoint and a + b = n;
#' 
#' case 5: A and B are disjoint and a + b < n;
#' 
#' case 6: A and B are incompatible.
pclades <- function(a, b, n, case = 1, log = FALSE){
  ans <- switch (case,
                 "1" = pclade(k = a, n = n, log = TRUE), ##  A == B
                 "2" = Rn(a = a, b = b, n = n, log = TRUE), ## A \in B
                 "3" = Rn(a = a, b = b, n = n, log = TRUE), ## B \in A
                 "4" = log(2)-log(n-1) -lchoose(n = n, k = a), ## intersect(A, B) = NULL and a + b = n 
                 "5" = rrn(a = a, b = b, n = n, log = TRUE), ## intersect(A, B) = NULL and a + b < n 
                 "6" =  -Inf ## A and B are incompatible
  )
  if(!log) ans <- exp(ans)
  return(ans)
}