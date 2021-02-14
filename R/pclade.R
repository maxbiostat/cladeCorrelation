#' Probability of a clade of size \code{k} on \code{n} taxa.
#'
#' @param k size of clade.
#' @param n number of taxa.
#' @param log logical if TRUE (default is FALSE) returns log probability.
#'
#' @return (log) probability of a clade of size \code{k} for \code{n} taxa.
#' @export pclade
#' @details Lemma 4.2 of Zhu, Degnan & Steel (2011).
#' @references Zhu, S., Degnan, J. H., & Steel, M. (2011).
#'  Clades, clans, and reciprocal monophyly under neutral evolutionary models.
#'  Theoretical population biology, 79(4), 220-227.
#' @examples
#' pclade(k = 2, n = 5)
pclade <- function(k, n, log = FALSE){
  if(k < 2 || k > n -1){
    ans <- 0
  }else{
    ans <- log(2) + log(n) - (log(k) + log(k + 1)) -lchoose(n, k)
    if(!log) ans <- exp(ans)
    return(ans)
  }
}
pclade <- Vectorize(pclade)

#' Probability of a clade under a PDA model on \code{n} taxa.
#'
#' @param clade a string containing a clade.
#' @param n number of taxa.
#' @param log logical. If \code{TRUE} returns the log probability (default is \code{FALSE}).
#'
#' @return (log) probability
#' @export prob_clade
#' @seealso \code{\link[cladeCorrelation]{pclade}}
#' @examples
#' all.c4 <- make_all_clades(4)
#' prob_clade(all.c4, n = 4)
prob_clade <- function(clade, n, log = FALSE){
  size <- length(get_clade_elements(clade))
  ans <- pclade(k = size, n = n, log = log)
  return(ans)
}
prob_clade <- Vectorize(prob_clade)